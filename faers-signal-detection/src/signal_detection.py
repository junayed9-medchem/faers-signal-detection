"""
signal_detection.py
-------------------
Disproportionality analysis methods for pharmacovigilance signal detection.

Methods implemented:
  - Proportional Reporting Ratio (PRR)        Evans et al., 2001
  - Reporting Odds Ratio (ROR)                van Puijenbroek et al., 2002
  - Chi-squared test (Pearson)

Signal threshold (Evans 2001):
  PRR >= 2  AND  chi² >= 4  AND  n >= 3
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import List


class SignalDetector:
    """
    Computes PRR, ROR, and chi-squared for each drug–adverse event pair
    using the 2x2 contingency table:

                    AE present    AE absent   |  Total
        Drug          a             b          |  a + b
        No drug       c             d          |  c + d
        ─────────────────────────────────────────────────
        Total         a+c           b+d        |  N

    Where:
        a  = reports with BOTH drug X and adverse event Y
        b  = reports with drug X but NOT adverse event Y
        c  = reports with adverse event Y but NOT drug X
        d  = reports with NEITHER drug X nor adverse event Y
        N  = total reports in database
    """

    def __init__(
        self,
        prr_threshold: float = 2.0,
        chi2_threshold: float = 4.0,
        min_count: int = 3,
    ):
        self.prr_threshold = prr_threshold
        self.chi2_threshold = chi2_threshold
        self.min_count = min_count

    # ── core computation ───────────────────────────────────────────────────────

    def compute_signal(
        self,
        drug_name: str,
        ae_name: str,
        a: int,
        drug_total: int,
        ae_total: int,
        n_total: int,
    ) -> dict:
        """
        Compute disproportionality metrics for a single drug–AE pair.

        Parameters
        ----------
        drug_name  : str  drug identifier
        ae_name    : str  adverse event (MedDRA PT)
        a          : int  reports with drug AND ae
        drug_total : int  a + b  (total drug reports)
        ae_total   : int  a + c  (total AE reports across all drugs)
        n_total    : int  N      (total FAERS reports)
        """
        a = int(a)
        ab = int(drug_total)   # a + b
        ac = int(ae_total)     # a + c
        N = int(n_total)

        b = ab - a             # drug present, AE absent
        c = ac - a             # drug absent,  AE present
        d = N - a - b - c      # drug absent,  AE absent

        null_result = {
            "drug": drug_name, "ae_name": ae_name, "count": a,
            "PRR": np.nan, "PRR_lower": np.nan, "PRR_upper": np.nan,
            "ROR": np.nan, "ROR_lower": np.nan, "ROR_upper": np.nan,
            "chi2": np.nan, "p_value": np.nan, "signal": False,
        }

        if any(x <= 0 for x in [a, b, c, d]):
            return null_result

        # ── Proportional Reporting Ratio ───────────────────────────────────────
        # PRR = [a / (a+b)] / [c / (c+d)]
        prr = (a / ab) / (c / (c + d))
        se_log_prr = np.sqrt(1/a - 1/ab + 1/c - 1/(c + d))
        prr_lo = np.exp(np.log(prr) - 1.96 * se_log_prr)
        prr_hi = np.exp(np.log(prr) + 1.96 * se_log_prr)

        # ── Reporting Odds Ratio ───────────────────────────────────────────────
        # ROR = (a * d) / (b * c)
        ror = (a * d) / (b * c)
        se_log_ror = np.sqrt(1/a + 1/b + 1/c + 1/d)
        ror_lo = np.exp(np.log(ror) - 1.96 * se_log_ror)
        ror_hi = np.exp(np.log(ror) + 1.96 * se_log_ror)

        # ── Chi-squared ────────────────────────────────────────────────────────
        contingency = np.array([[a, b], [c, d]])
        chi2, p_val, _, _ = stats.chi2_contingency(contingency, correction=False)

        # ── Signal flag ────────────────────────────────────────────────────────
        is_signal = (
            prr >= self.prr_threshold
            and chi2 >= self.chi2_threshold
            and a >= self.min_count
        )

        return {
            "drug": drug_name,
            "ae_name": ae_name,
            "count": a,
            "PRR": round(prr, 3),
            "PRR_lower": round(prr_lo, 3),
            "PRR_upper": round(prr_hi, 3),
            "ROR": round(ror, 3),
            "ROR_lower": round(ror_lo, 3),
            "ROR_upper": round(ror_hi, 3),
            "chi2": round(chi2, 3),
            "p_value": round(p_val, 6),
            "signal": is_signal,
        }

    # ── batch ──────────────────────────────────────────────────────────────────

    def run_batch(
        self,
        drug_name: str,
        top_aes: list,
        drug_total: int,
        n_total: int,
        ae_background: dict,
    ) -> pd.DataFrame:
        """
        Compute signals for a list of adverse events.

        Parameters
        ----------
        top_aes       : list of {term, count} dicts from openFDA
        ae_background : dict mapping ae_term -> background total
        """
        rows = []
        for ae in top_aes:
            ae_name = ae["term"]
            a = ae["count"]
            ae_total = ae_background.get(ae_name, 0)
            row = self.compute_signal(drug_name, ae_name, a, drug_total, ae_total, n_total)
            rows.append(row)

        df = pd.DataFrame(rows)
        return df.sort_values("PRR", ascending=False).reset_index(drop=True)

    # ── reporting ──────────────────────────────────────────────────────────────

    def print_summary(self, df: pd.DataFrame) -> None:
        total = len(df)
        n_signals = int(df["signal"].sum())
        print(f"\n{'=' * 55}")
        print(f"  Signal Detection Summary")
        print(f"{'=' * 55}")
        print(f"  AEs analysed    : {total}")
        print(f"  Signals detected: {n_signals}  ({100 * n_signals / total:.1f}%)")
        print(f"  Threshold       : PRR≥{self.prr_threshold}, χ²≥{self.chi2_threshold}, n≥{self.min_count}")
        print(f"{'=' * 55}")
        if n_signals > 0:
            cols = ["ae_name", "count", "PRR", "ROR", "chi2"]
            print(df[df["signal"]][cols].head(10).to_string(index=False))
        print()
