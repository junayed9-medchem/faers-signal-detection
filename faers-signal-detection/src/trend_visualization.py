"""
trend_visualization.py
----------------------
Visualisations for longitudinal FAERS trend analysis.

Plots:
  1. Annual report volume line chart (with YoY % change)
  2. Top AE composition over time (multi-line)
  3. PRR trend lines for key signals
  4. PRR heatmap: year × AE
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

DPI          = 150
GRID_ALPHA   = 0.25
BACKGROUND   = "#FAFAFA"
SIGNAL_COLOR = "#C0392B"
ACCENT       = "#2980B9"


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor(BACKGROUND)
    ax.grid(alpha=GRID_ALPHA)


class TrendVisualizer:

    def __init__(self, output_dir: str = "outputs/figures"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": DPI})

    # ── 1. Report volume trend ─────────────────────────────────────────────────

    def plot_report_volume(self, drug_totals: pd.DataFrame, drug_name: str) -> None:
        """
        Dual-axis chart: absolute report count (bar) and
        year-over-year % change (line).
        """
        fig, ax1 = plt.subplots(figsize=(10, 5))
        ax2 = ax1.twinx()

        years  = drug_totals["year"]
        counts = drug_totals["total_reports"]
        yoy    = drug_totals["yoy_change_pct"]

        ax1.bar(years, counts, color=ACCENT, alpha=0.75,
                edgecolor="white", label="Annual reports")
        ax1.set_ylabel("Total FAERS Reports", fontsize=11, color=ACCENT)
        ax1.tick_params(axis="y", labelcolor=ACCENT)

        ax2.plot(years, yoy, color=SIGNAL_COLOR, marker="o",
                 linewidth=2, markersize=6, label="YoY % change")
        ax2.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
        ax2.set_ylabel("Year-over-Year Change (%)", fontsize=11,
                       color=SIGNAL_COLOR)
        ax2.tick_params(axis="y", labelcolor=SIGNAL_COLOR)

        ax1.set_xlabel("Calendar Year", fontsize=11)
        ax1.set_xticks(years)
        fig.suptitle(
            f"Annual FAERS Report Volume — {drug_name.title()}\n"
            f"Post-Marketing Surveillance Trend",
            fontsize=13, fontweight="bold", y=1.01,
        )

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2,
                   loc="upper left", fontsize=10)
        _style(ax1)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_report_volume.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 2. AE composition trend ────────────────────────────────────────────────

    def plot_ae_trends(self, ae_counts: pd.DataFrame, drug_name: str) -> None:
        """
        Multi-line chart of report counts per AE per year.
        Shows how the AE reporting landscape evolves over time.
        """
        ae_names = ae_counts["ae_name"].unique()
        colors   = cm.tab10(np.linspace(0, 0.9, len(ae_names)))

        fig, ax = plt.subplots(figsize=(11, 6))

        for ae, color in zip(ae_names, colors):
            sub = ae_counts[ae_counts["ae_name"] == ae].sort_values("year")
            ax.plot(sub["year"], sub["count"],
                    marker="o", linewidth=2, markersize=6,
                    label=ae.title(), color=color, alpha=0.85)

        ax.set_xlabel("Calendar Year", fontsize=11)
        ax.set_ylabel("Co-reported Cases with Drug", fontsize=11)
        ax.set_title(
            f"Adverse Event Reporting Trends — {drug_name.title()}\n"
            f"Annual co-report counts per MedDRA preferred term",
            fontsize=13, fontweight="bold", pad=12,
        )
        ax.legend(loc="upper left", fontsize=9, framealpha=0.9,
                  ncol=2 if len(ae_names) > 5 else 1)
        ax.set_xticks(ae_counts["year"].unique())
        _style(ax)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_ae_trends.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 3. PRR trend lines ────────────────────────────────────────────────────

    def plot_prr_trends(self, prr_df: pd.DataFrame, drug_name: str) -> None:
        """
        PRR trend line per AE across years.
        Signals (PRR ≥ 2) shown in red, others in blue.
        """
        ae_names = prr_df["ae_name"].unique()
        colors   = cm.tab10(np.linspace(0, 0.9, len(ae_names)))

        fig, ax = plt.subplots(figsize=(11, 6))

        for ae, color in zip(ae_names, colors):
            sub = prr_df[prr_df["ae_name"] == ae].sort_values("year")
            sub_clean = sub.dropna(subset=["PRR"])
            if sub_clean.empty:
                continue
            is_signal = sub_clean["signal"].any()
            lw = 2.5 if is_signal else 1.5
            ls = "-"  if is_signal else "--"
            ax.plot(sub_clean["year"], sub_clean["PRR"],
                    marker="o", linewidth=lw, linestyle=ls,
                    markersize=6, label=ae.title(),
                    color=SIGNAL_COLOR if is_signal else color,
                    alpha=0.85)

        # Signal threshold
        ax.axhline(2.0, color="black", linestyle=":", linewidth=1.2,
                   alpha=0.6, label="Signal threshold (PRR = 2)")

        ax.set_xlabel("Calendar Year", fontsize=11)
        ax.set_ylabel("Proportional Reporting Ratio (PRR)", fontsize=11)
        ax.set_title(
            f"PRR Trend Over Time — {drug_name.title()}\n"
            f"Solid red lines = confirmed signals  ·  dashed = below threshold",
            fontsize=13, fontweight="bold", pad=12,
        )
        ax.legend(loc="upper left", fontsize=9, framealpha=0.9,
                  ncol=2 if len(ae_names) > 5 else 1)
        ax.set_xticks(prr_df["year"].unique())
        _style(ax)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_prr_trends.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 4. PRR heatmap ─────────────────────────────────────────────────────────

    def plot_prr_heatmap(self, prr_df: pd.DataFrame, drug_name: str) -> None:
        """
        Heatmap: rows = AEs, columns = years, cells = PRR value.
        Colour intensity shows signal strength across time.
        """
        pivot = prr_df.pivot_table(
            index="ae_name", columns="year", values="PRR", aggfunc="first"
        ).fillna(0)

        fig, ax = plt.subplots(figsize=(10, max(5, len(pivot) * 0.55)))
        im = ax.imshow(pivot.values, aspect="auto",
                       cmap="YlOrRd", vmin=0,
                       vmax=min(pivot.values.max(), 20))

        # Axis labels
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns.astype(int), fontsize=11)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index, fontsize=9)

        # Cell text + signal star
        signal_lookup = prr_df.groupby(["ae_name", "year"])["signal"].first()
        for i, ae in enumerate(pivot.index):
            for j, yr in enumerate(pivot.columns):
                val = pivot.iloc[i, j]
                is_sig = signal_lookup.get((ae, yr), False)
                star = " ★" if is_sig else ""
                text_color = "white" if val > pivot.values.max() * 0.6 else "black"
                ax.text(j, i, f"{val:.1f}{star}",
                        ha="center", va="center",
                        fontsize=8, color=text_color, fontweight="bold")

        plt.colorbar(im, ax=ax, fraction=0.025,
                     label="PRR  (★ = confirmed signal)")
        ax.set_title(
            f"PRR Heatmap — {drug_name.title()}  (Year × Adverse Event)\n"
            f"★ = PRR≥2, χ²≥4, n≥3",
            fontsize=12, fontweight="bold", pad=14,
        )
        ax.spines[:].set_visible(False)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_prr_heatmap.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 5. Full longitudinal dashboard ────────────────────────────────────────

    def plot_longitudinal_dashboard(
        self,
        drug_totals: pd.DataFrame,
        ae_counts: pd.DataFrame,
        prr_df: pd.DataFrame,
        drug_name: str,
    ) -> None:
        """
        2×2 summary dashboard combining all longitudinal views.
        Ideal for README screenshot and presentations.
        """
        fig, axes = plt.subplots(2, 2, figsize=(16, 11))
        fig.suptitle(
            f"Longitudinal FAERS Analysis — {drug_name.title()}\n"
            f"Post-Marketing Safety Surveillance  ·  "
            f"{int(drug_totals['year'].min())}–{int(drug_totals['year'].max())}",
            fontsize=14, fontweight="bold", y=1.0,
        )

        years  = drug_totals["year"]
        counts = drug_totals["total_reports"]

        # (0,0) Report volume bars
        axes[0, 0].bar(years, counts, color=ACCENT, alpha=0.8, edgecolor="white")
        axes[0, 0].set_title("Annual Report Volume", fontsize=11, fontweight="bold")
        axes[0, 0].set_xlabel("Year"); axes[0, 0].set_ylabel("Reports")
        axes[0, 0].set_xticks(years)
        _style(axes[0, 0])

        # (0,1) YoY % change
        yoy = drug_totals["yoy_change_pct"].fillna(0)
        bar_colors = [SIGNAL_COLOR if v > 0 else ACCENT for v in yoy]
        axes[0, 1].bar(years, yoy, color=bar_colors, alpha=0.8, edgecolor="white")
        axes[0, 1].axhline(0, color="black", linewidth=0.8)
        axes[0, 1].set_title("Year-over-Year Change (%)", fontsize=11, fontweight="bold")
        axes[0, 1].set_xlabel("Year"); axes[0, 1].set_ylabel("% Change")
        axes[0, 1].set_xticks(years)
        _style(axes[0, 1])

        # (1,0) AE composition trend (top 5)
        ae_names = (
            ae_counts.groupby("ae_name")["count"].sum()
            .nlargest(5).index.tolist()
        )
        colors = cm.tab10(np.linspace(0, 0.9, len(ae_names)))
        for ae, color in zip(ae_names, colors):
            sub = ae_counts[ae_counts["ae_name"] == ae].sort_values("year")
            axes[1, 0].plot(sub["year"], sub["count"],
                            marker="o", linewidth=2, color=color,
                            label=ae.title(), alpha=0.85)
        axes[1, 0].set_title("Top AE Trends", fontsize=11, fontweight="bold")
        axes[1, 0].set_xlabel("Year"); axes[1, 0].set_ylabel("Co-reports")
        axes[1, 0].legend(fontsize=7, ncol=2)
        axes[1, 0].set_xticks(years)
        _style(axes[1, 0])

        # (1,1) PRR trends
        ae_names_prr = (
            prr_df.groupby("ae_name")["PRR"].mean()
            .nlargest(5).index.tolist()
        )
        colors2 = cm.tab10(np.linspace(0, 0.9, len(ae_names_prr)))
        for ae, color in zip(ae_names_prr, colors2):
            sub = prr_df[prr_df["ae_name"] == ae].sort_values("year").dropna(subset=["PRR"])
            axes[1, 1].plot(sub["year"], sub["PRR"],
                            marker="o", linewidth=2, color=color,
                            label=ae.title(), alpha=0.85)
        axes[1, 1].axhline(2.0, color="black", linestyle=":",
                            linewidth=1, alpha=0.6, label="PRR=2")
        axes[1, 1].set_title("PRR Trends", fontsize=11, fontweight="bold")
        axes[1, 1].set_xlabel("Year"); axes[1, 1].set_ylabel("PRR")
        axes[1, 1].legend(fontsize=7, ncol=2)
        axes[1, 1].set_xticks(years)
        _style(axes[1, 1])

        fig.tight_layout()
        path = os.path.join(
            self.output_dir, f"{drug_name}_longitudinal_dashboard.png"
        )
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")
