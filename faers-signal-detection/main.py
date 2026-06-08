"""
main.py
-------
FAERS Pharmacovigilance Signal Detection Pipeline

Usage:
    python main.py                             # snapshot analysis (metformin)
    python main.py --drug warfarin             # different drug
    python main.py --drug aspirin --top 30     # more AEs
    python main.py --drug metformin --trend    # + longitudinal analysis
    python main.py --drug metformin --trend --years 2019 2020 2021 2022 2023
"""

import os
import sys
import argparse
import logging
import pandas as pd

# ── add project root to path so `src` is importable ───────────────────────────
sys.path.insert(0, os.path.dirname(__file__))

from src.api_client          import FAERSClient
from src.signal_detection    import SignalDetector
from src.visualization       import SignalVisualizer
from src.trend_analysis      import TrendAnalyzer
from src.trend_visualization import TrendVisualizer

# ── logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="FAERS disproportionality analysis using the openFDA API"
    )
    parser.add_argument("--drug",    type=str, default="metformin",
                        help="Generic drug name to analyse (default: metformin)")
    parser.add_argument("--top",     type=int, default=25,
                        help="Number of top adverse events to evaluate (default: 25)")
    parser.add_argument("--api_key", type=str, default=None,
                        help="openFDA API key (optional — raises rate limit)")
    parser.add_argument("--out",     type=str, default="outputs",
                        help="Output directory (default: outputs/)")
    parser.add_argument("--trend",   action="store_true",
                        help="Run longitudinal trend analysis across years")
    parser.add_argument("--years",   type=int, nargs="+",
                        default=[2019, 2020, 2021, 2022, 2023],
                        help="Calendar years for trend analysis (default: 2019-2023)")
    return parser.parse_args()


def main():
    args = parse_args()
    drug_name = args.drug.lower()

    print(f"\n{'='*60}")
    print(f"  FAERS Signal Detection Pipeline")
    print(f"  Drug        : {drug_name.title()}")
    print(f"  Top AEs     : {args.top}")
    print(f"  Output dir  : {args.out}")
    print(f"{'='*60}\n")

    # ── initialise components ─────────────────────────────────────────────────
    client  = FAERSClient(api_key=args.api_key)
    detector = SignalDetector(prr_threshold=2.0, chi2_threshold=4.0, min_count=3)
    viz      = SignalVisualizer(output_dir=os.path.join(args.out, "figures"))

    # ── Step 1: background totals ─────────────────────────────────────────────
    logger.info("Fetching background totals from openFDA…")
    n_total  = client.get_total_reports()
    logger.info(f"  Total FAERS reports       : {n_total:,}")

    drug_total = client.get_drug_total(drug_name)
    if drug_total == 0:
        print(f"\nERROR: No reports found for drug '{drug_name}'. "
              f"Check spelling or try a different generic name.\n")
        sys.exit(1)
    logger.info(f"  Reports for {drug_name:<15} : {drug_total:,}")

    # ── Step 2: top adverse events ────────────────────────────────────────────
    logger.info(f"Fetching top {args.top} adverse events…")
    top_aes = client.get_top_adverse_events(drug_name, limit=args.top)
    if not top_aes:
        print(f"ERROR: Could not retrieve adverse events for '{drug_name}'.\n")
        sys.exit(1)
    logger.info(f"  Retrieved {len(top_aes)} adverse events")

    # ── Step 3: background count per AE ──────────────────────────────────────
    logger.info("Fetching AE background counts (one call per AE)…")
    ae_background = {}
    for i, ae in enumerate(top_aes, 1):
        ae_name = ae["term"]
        ae_background[ae_name] = client.get_ae_background_total(ae_name)
        if i % 5 == 0:
            logger.info(f"  {i}/{len(top_aes)} AEs processed")

    # ── Step 4: compute signals ───────────────────────────────────────────────
    logger.info("Computing PRR, ROR, chi-squared…")
    df = detector.run_batch(drug_name, top_aes, drug_total, n_total, ae_background)

    # ── Step 5: save results ──────────────────────────────────────────────────
    os.makedirs(args.out, exist_ok=True)
    csv_path = os.path.join(args.out, f"{drug_name}_signals.csv")
    df.to_csv(csv_path, index=False)
    logger.info(f"Results saved → {csv_path}")

    # ── Step 6: print summary ─────────────────────────────────────────────────
    detector.print_summary(df)

    # ── Step 7: generate visualisations ──────────────────────────────────────
    logger.info("Generating visualisations…")
    viz.plot_prr_bar(df, drug_name)
    viz.plot_signal_scatter(df, drug_name)
    viz.plot_signal_heatmap(df, drug_name)
    viz.plot_summary_dashboard(df, drug_name)

    print(f"\nAll outputs saved to: {args.out}/\n")

    # ── Step 8: longitudinal trend analysis (optional) ────────────────────────
    if args.trend:
        print(f"\n{'='*60}")
        print(f"  Longitudinal Trend Analysis  ({args.years[0]}–{args.years[-1]})")
        print(f"{'='*60}\n")

        # Use top 5 detected signals for trend tracking
        top_signal_aes = (
            df[df["signal"]].nlargest(5, "PRR")["ae_name"].tolist()
            if df["signal"].any()
            else df.nlargest(5, "PRR")["ae_name"].tolist()
        )
        if not top_signal_aes:
            print("No signals detected — skipping trend analysis.\n")
        else:
            analyzer  = TrendAnalyzer(client, detector)
            trend_viz = TrendVisualizer(
                output_dir=os.path.join(args.out, "figures")
            )

            results = analyzer.run(drug_name, args.years, top_signal_aes)

            # Save CSVs
            results["drug_totals"].to_csv(
                os.path.join(args.out, f"{drug_name}_yearly_volume.csv"),
                index=False,
            )
            results["prr_trends"].to_csv(
                os.path.join(args.out, f"{drug_name}_prr_trends.csv"),
                index=False,
            )

            # Generate trend plots
            logger.info("Generating longitudinal visualisations…")
            trend_viz.plot_report_volume(results["drug_totals"], drug_name)
            trend_viz.plot_ae_trends(results["ae_counts"], drug_name)
            trend_viz.plot_prr_trends(results["prr_trends"], drug_name)
            trend_viz.plot_prr_heatmap(results["prr_trends"], drug_name)
            trend_viz.plot_longitudinal_dashboard(
                results["drug_totals"],
                results["ae_counts"],
                results["prr_trends"],
                drug_name,
            )

            print(f"\nLongitudinal outputs saved to: {args.out}/\n")


if __name__ == "__main__":
    main()
