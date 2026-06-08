"""
visualization.py
----------------
Publication-ready plots for pharmacovigilance signal detection results.

Plots generated:
  1. PRR horizontal bar chart with 95% CI
  2. PRR vs ROR scatter plot
  3. Signal metrics heatmap
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── colour palette ─────────────────────────────────────────────────────────────
SIGNAL_COLOR    = "#C0392B"   # red   — confirmed signal
NO_SIGNAL_COLOR = "#2980B9"   # blue  — below threshold
BACKGROUND      = "#FAFAFA"
GRID_ALPHA      = 0.25
DPI             = 150


def _style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor(BACKGROUND)
    ax.grid(alpha=GRID_ALPHA)


class SignalVisualizer:

    def __init__(self, output_dir: str = "outputs/figures"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": DPI})

    # ── 1. PRR bar chart ───────────────────────────────────────────────────────

    def plot_prr_bar(self, df: pd.DataFrame, drug_name: str) -> None:
        """Horizontal bar chart of PRR ± 95% CI for the top 20 AEs."""
        plot_df = df.dropna(subset=["PRR"]).nlargest(20, "PRR").copy()
        colors = [SIGNAL_COLOR if s else NO_SIGNAL_COLOR for s in plot_df["signal"]]

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.barh(plot_df["ae_name"], plot_df["PRR"],
                color=colors, alpha=0.85, edgecolor="white", height=0.65)

        # 95% CI error bars
        xerr_lo = (plot_df["PRR"] - plot_df["PRR_lower"]).clip(lower=0)
        xerr_hi = (plot_df["PRR_upper"] - plot_df["PRR"]).clip(lower=0)
        ax.errorbar(
            plot_df["PRR"], plot_df["ae_name"],
            xerr=[xerr_lo, xerr_hi],
            fmt="none", color="black", capsize=3, linewidth=0.8, alpha=0.55,
        )

        # Signal threshold line
        ax.axvline(2.0, color="black", linestyle="--", linewidth=1.2, alpha=0.6,
                   label="Signal threshold (PRR = 2)")

        # Legend
        patches = [
            mpatches.Patch(color=SIGNAL_COLOR,    label="Signal  (PRR≥2, χ²≥4, n≥3)"),
            mpatches.Patch(color=NO_SIGNAL_COLOR, label="Below threshold"),
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=10, framealpha=0.9)

        ax.set_xlabel("Proportional Reporting Ratio (PRR)", fontsize=12)
        ax.set_title(
            f"Pharmacovigilance Signal Detection — {drug_name.title()}\n"
            f"FAERS Disproportionality Analysis · PRR ± 95% CI",
            fontsize=13, fontweight="bold", pad=14,
        )
        _style_ax(ax)
        ax.set_xlim(left=0)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_prr_bar.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 2. PRR vs ROR scatter ──────────────────────────────────────────────────

    def plot_signal_scatter(self, df: pd.DataFrame, drug_name: str) -> None:
        """Scatter plot: PRR (x) vs ROR (y), bubble size ∝ report count."""
        plot_df = df.dropna(subset=["PRR", "ROR"]).copy()
        colors = [SIGNAL_COLOR if s else NO_SIGNAL_COLOR for s in plot_df["signal"]]
        sizes  = np.clip(plot_df["count"] / plot_df["count"].max() * 500, 25, 500)

        fig, ax = plt.subplots(figsize=(10, 7))
        ax.scatter(plot_df["PRR"], plot_df["ROR"],
                   c=colors, s=sizes, alpha=0.75, edgecolors="white", linewidths=0.5)

        ax.axvline(2.0, color="grey", linestyle="--", linewidth=1, alpha=0.6)
        ax.axhline(2.0, color="grey", linestyle=":",  linewidth=1, alpha=0.6)

        # Annotate top signals
        for _, row in plot_df[plot_df["signal"]].nlargest(8, "PRR").iterrows():
            ax.annotate(
                row["ae_name"], xy=(row["PRR"], row["ROR"]),
                xytext=(5, 4), textcoords="offset points",
                fontsize=7, alpha=0.85,
            )

        patches = [
            mpatches.Patch(color=SIGNAL_COLOR,    label="Signal detected"),
            mpatches.Patch(color=NO_SIGNAL_COLOR, label="Below threshold"),
        ]
        ax.legend(handles=patches, fontsize=10, framealpha=0.9)
        ax.set_xlabel("Proportional Reporting Ratio (PRR)", fontsize=12)
        ax.set_ylabel("Reporting Odds Ratio (ROR)", fontsize=12)
        ax.set_title(
            f"PRR vs ROR — {drug_name.title()}  (bubble size ∝ report count)",
            fontsize=13, fontweight="bold", pad=14,
        )
        _style_ax(ax)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_scatter.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 3. Heatmap ─────────────────────────────────────────────────────────────

    def plot_signal_heatmap(self, df: pd.DataFrame, drug_name: str) -> None:
        """Colour-coded table of PRR, ROR, and chi-squared for top 15 AEs."""
        plot_df = df.dropna(subset=["PRR"]).nlargest(15, "PRR").copy()
        metrics = ["PRR", "ROR", "chi2"]
        data    = plot_df[metrics].values.astype(float)

        # Normalise per column for colour intensity
        col_min = data.min(axis=0)
        col_max = data.max(axis=0)
        data_norm = (data - col_min) / (col_max - col_min + 1e-8)

        fig, ax = plt.subplots(figsize=(9, 9))
        ax.imshow(data_norm, aspect="auto", cmap="YlOrRd", alpha=0.85, vmin=0, vmax=1)

        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(["PRR", "ROR", "Chi²"], fontsize=12, fontweight="bold")
        ax.set_yticks(range(len(plot_df)))
        ax.set_yticklabels(plot_df["ae_name"], fontsize=9)

        # Cell text
        for i in range(len(plot_df)):
            for j, col in enumerate(metrics):
                val = plot_df.iloc[i][col]
                star = " ★" if (col == "PRR" and plot_df.iloc[i]["signal"]) else ""
                text_color = "white" if data_norm[i, j] > 0.65 else "black"
                ax.text(j, i, f"{val:.1f}{star}", ha="center", va="center",
                        fontsize=8, color=text_color, fontweight="bold")

        ax.set_title(
            f"Signal Metrics Heatmap — {drug_name.title()}\n"
            f"★ = confirmed signal  ·  colour intensity ∝ metric value",
            fontsize=12, fontweight="bold", pad=14,
        )
        ax.spines[:].set_visible(False)

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_heatmap.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")

    # ── 4. Summary report ──────────────────────────────────────────────────────

    def plot_summary_dashboard(self, df: pd.DataFrame, drug_name: str) -> None:
        """2×2 dashboard combining key views in one figure."""
        df_clean = df.dropna(subset=["PRR"])
        top20    = df_clean.nlargest(20, "PRR")
        colors   = [SIGNAL_COLOR if s else NO_SIGNAL_COLOR for s in top20["signal"]]

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(
            f"FAERS Pharmacovigilance Signal Detection — {drug_name.title()}",
            fontsize=15, fontweight="bold", y=0.98,
        )

        # (0,0) PRR bar
        axes[0, 0].barh(top20["ae_name"], top20["PRR"], color=colors, alpha=0.85, edgecolor="white")
        axes[0, 0].axvline(2.0, color="black", linestyle="--", linewidth=1)
        axes[0, 0].set_xlabel("PRR")
        axes[0, 0].set_title("PRR — Top 20 Adverse Events")
        _style_ax(axes[0, 0])

        # (0,1) Signal count donut
        n_signal   = int(df["signal"].sum())
        n_nosignal = len(df) - n_signal
        axes[0, 1].pie(
            [n_signal, n_nosignal],
            labels=[f"Signals ({n_signal})", f"Below threshold ({n_nosignal})"],
            colors=[SIGNAL_COLOR, NO_SIGNAL_COLOR],
            autopct="%1.0f%%", startangle=90, textprops={"fontsize": 11},
            wedgeprops={"edgecolor": "white", "linewidth": 2},
        )
        axes[0, 1].set_title("Signal vs Non-Signal AEs")

        # (1,0) PRR vs count scatter
        axes[1, 0].scatter(df_clean["count"], df_clean["PRR"],
                           c=[SIGNAL_COLOR if s else NO_SIGNAL_COLOR for s in df_clean["signal"]],
                           s=50, alpha=0.65, edgecolors="white")
        axes[1, 0].axhline(2.0, color="grey", linestyle="--", linewidth=1)
        axes[1, 0].set_xlabel("Report Count (n)")
        axes[1, 0].set_ylabel("PRR")
        axes[1, 0].set_title("PRR vs Report Volume")
        _style_ax(axes[1, 0])

        # (1,1) ROR vs PRR
        axes[1, 1].scatter(df_clean["PRR"], df_clean["ROR"],
                           c=[SIGNAL_COLOR if s else NO_SIGNAL_COLOR for s in df_clean["signal"]],
                           s=50, alpha=0.65, edgecolors="white")
        axes[1, 1].axvline(2.0, color="grey", linestyle="--", linewidth=1)
        axes[1, 1].axhline(2.0, color="grey", linestyle=":",  linewidth=1)
        axes[1, 1].set_xlabel("PRR")
        axes[1, 1].set_ylabel("ROR")
        axes[1, 1].set_title("PRR vs ROR Concordance")
        _style_ax(axes[1, 1])

        fig.tight_layout()
        path = os.path.join(self.output_dir, f"{drug_name}_dashboard.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {path}")
