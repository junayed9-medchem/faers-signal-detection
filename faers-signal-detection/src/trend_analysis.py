"""
trend_analysis.py
-----------------
Longitudinal analysis of FAERS adverse event reporting trends.

Tracks how drug safety reporting patterns change year-over-year —
mirroring FDA's approach to RWE submission monitoring (FY 2020-2024).

Analyses:
  1. Annual report volume for the drug
  2. Top AE composition changes over time
  3. PRR trend for key safety signals (strengthening vs weakening)
  4. Year-over-year reporting rate change
"""

import logging
import pandas as pd
import numpy as np
from typing import List, Dict

logger = logging.getLogger(__name__)


class TrendAnalyzer:
    """
    Builds a longitudinal view of FAERS drug safety reporting by
    querying the openFDA API for each year in a specified range.

    Each year's data forms one row in the output DataFrames, enabling
    trend analysis across the observation window.
    """

    def __init__(self, client, detector):
        """
        Parameters
        ----------
        client   : FAERSClient  — rate-limited openFDA API wrapper
        detector : SignalDetector — PRR/ROR computation engine
        """
        self.client   = client
        self.detector = detector

    # ── internal helpers ───────────────────────────────────────────────────────

    def _date_filter(self, year: int) -> str:
        """Return an openFDA date filter for a full calendar year."""
        return f"receivedate:[{year}0101+TO+{year}1231]"

    def _count(self, search: str) -> int:
        """Return total report count for a search expression, or 0."""
        return self.client._total(search)

    # ── public methods ─────────────────────────────────────────────────────────

    def get_yearly_drug_totals(
        self, drug_name: str, years: List[int]
    ) -> pd.DataFrame:
        """
        Total FAERS reports mentioning *drug_name* for each year.

        Returns
        -------
        DataFrame with columns: year, total_reports, yoy_change_pct
        """
        logger.info("  Fetching yearly report volume for %s…", drug_name)
        rows = []
        for year in years:
            search = (
                f'patient.drug.openfda.generic_name:"{drug_name}"'
                f"+AND+{self._date_filter(year)}"
            )
            count = self._count(search)
            logger.info("    %d : %s reports", year, f"{count:,}")
            rows.append({"year": year, "total_reports": count})

        df = pd.DataFrame(rows)
        df["yoy_change_pct"] = df["total_reports"].pct_change() * 100
        return df

    def get_yearly_ae_counts(
        self, drug_name: str, ae_names: List[str], years: List[int]
    ) -> pd.DataFrame:
        """
        Report count per (year, AE) for a given drug.

        Returns
        -------
        DataFrame with columns: year, ae_name, count
        """
        logger.info("  Fetching yearly AE counts (%d AEs × %d years)…",
                    len(ae_names), len(years))
        rows = []
        for year in years:
            for ae in ae_names:
                ae_safe = ae.replace('"', '\\"')
                search = (
                    f'patient.drug.openfda.generic_name:"{drug_name}"'
                    f'+AND+patient.reaction.reactionmeddrapt.exact:"{ae_safe}"'
                    f"+AND+{self._date_filter(year)}"
                )
                count = self._count(search)
                rows.append({"year": year, "ae_name": ae, "count": count})

        return pd.DataFrame(rows)

    def compute_prr_trends(
        self, drug_name: str, ae_names: List[str], years: List[int]
    ) -> pd.DataFrame:
        """
        Compute PRR and ROR for each (year, AE) combination.

        A rising PRR trend suggests a strengthening safety signal;
        a flat or falling trend may indicate reporting bias stabilising.

        Returns
        -------
        DataFrame with columns:
            year, ae_name, count, PRR, ROR, chi2, signal
        """
        logger.info("  Computing PRR trends…")

        # ── yearly drug totals (a + b) ────────────────────────────────────────
        drug_totals_df = self.get_yearly_drug_totals(drug_name, years)
        drug_totals = dict(zip(drug_totals_df["year"],
                               drug_totals_df["total_reports"]))

        # ── yearly background totals (N) ─────────────────────────────────────
        logger.info("  Fetching yearly background totals…")
        bg_totals: Dict[int, int] = {}
        for year in years:
            bg_totals[year] = self._count(self._date_filter(year))
            logger.info("    %d background: %s", year,
                        f"{bg_totals[year]:,}")

        # ── drug+AE counts per year (a) ───────────────────────────────────────
        ae_counts_df = self.get_yearly_ae_counts(drug_name, ae_names, years)

        # ── AE background per year (a + c) ────────────────────────────────────
        logger.info("  Fetching AE background counts by year…")
        ae_bg_rows = []
        for year in years:
            for ae in ae_names:
                ae_safe = ae.replace('"', '\\"')
                search = (
                    f'patient.reaction.reactionmeddrapt.exact:"{ae_safe}"'
                    f"+AND+{self._date_filter(year)}"
                )
                bg = self._count(search)
                ae_bg_rows.append({"year": year, "ae_name": ae, "ae_bg": bg})
        ae_bg_df = pd.DataFrame(ae_bg_rows)

        # ── merge and compute PRR ─────────────────────────────────────────────
        merged = ae_counts_df.merge(ae_bg_df, on=["year", "ae_name"])

        results = []
        for _, row in merged.iterrows():
            year = int(row["year"])
            signal = self.detector.compute_signal(
                drug_name=drug_name,
                ae_name=row["ae_name"],
                a=int(row["count"]),
                drug_total=drug_totals.get(year, 0),
                ae_total=int(row["ae_bg"]),
                n_total=bg_totals.get(year, 0),
            )
            signal["year"] = year
            results.append(signal)

        df = pd.DataFrame(results)
        return df.sort_values(["ae_name", "year"]).reset_index(drop=True)

    # ── full pipeline ──────────────────────────────────────────────────────────

    def run(
        self,
        drug_name: str,
        years: List[int],
        top_ae_names: List[str],
    ) -> dict:
        """
        Run the full longitudinal analysis pipeline.

        Returns
        -------
        dict with keys:
            'drug_totals' : yearly report volume DataFrame
            'ae_counts'   : yearly AE count DataFrame
            'prr_trends'  : yearly PRR / ROR DataFrame
        """
        print(f"\nLongitudinal analysis: {drug_name.title()}  "
              f"({years[0]}–{years[-1]})\n")

        drug_totals = self.get_yearly_drug_totals(drug_name, years)
        ae_counts   = self.get_yearly_ae_counts(drug_name, top_ae_names, years)
        prr_trends  = self.compute_prr_trends(drug_name, top_ae_names, years)

        # ── print summary ─────────────────────────────────────────────────────
        print("\nYearly Report Volume:")
        print(drug_totals.to_string(index=False))

        print("\nPRR Trend for Key Signals:")
        pivot = prr_trends.pivot_table(
            index="ae_name", columns="year", values="PRR", aggfunc="first"
        )
        print(pivot.round(2).to_string())

        return {
            "drug_totals": drug_totals,
            "ae_counts":   ae_counts,
            "prr_trends":  prr_trends,
        }
