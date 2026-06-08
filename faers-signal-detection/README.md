# FAERS Pharmacovigilance Signal Detection

A Python pipeline for post-marketing drug safety surveillance using the **FDA Adverse Event Reporting System (FAERS)** via the openFDA API. Implements standard disproportionality analysis methods used in regulatory science to detect drug–adverse event safety signals.

---

## Background

The **FDA Adverse Event Reporting System (FAERS)** is a real-world evidence (RWE) database containing millions of spontaneous adverse event reports submitted by healthcare professionals, patients, and manufacturers. Signal detection in FAERS is a cornerstone of FDA's post-marketing pharmacovigilance program.

This pipeline automates the **disproportionality analysis** workflow:

> If a drug–adverse event combination is reported *more frequently than expected by chance*, it is flagged as a **safety signal** warranting further regulatory review.

---

## Methods

All three measures use a **2×2 contingency table** built from FAERS report counts:

|              | AE present | AE absent | Total     |
|:-------------|:----------:|:---------:|:---------:|
| **Drug**     | a          | b         | a + b     |
| **No Drug**  | c          | d         | c + d     |
| **Total**    | a + c      | b + d     | N         |

### Proportional Reporting Ratio (PRR)
$$\text{PRR} = \frac{a \,/\, (a+b)}{c \,/\, (c+d)}$$

Compares the proportion of adverse event Y among drug X reports to the proportion of Y in all *other* drug reports. 95% confidence intervals computed on the log scale.

### Reporting Odds Ratio (ROR)
$$\text{ROR} = \frac{a \times d}{b \times c}$$

The odds of adverse event Y being reported with drug X relative to all other drugs. Analogous to the odds ratio in case-control studies.

### Chi-squared Test
Pearson chi-squared on the 2×2 table, without Yates' correction.

### Signal Criteria (Evans et al., 2001)
A drug–AE pair is flagged as a **signal** when all three conditions hold:
- PRR ≥ 2
- Chi-squared ≥ 4
- Report count (n) ≥ 3

---

## Project Structure

```
faers-signal-detection/
├── main.py                   # Pipeline entry point
├── requirements.txt
├── src/
│   ├── api_client.py         # openFDA API wrapper (rate-limited)
│   ├── signal_detection.py   # PRR, ROR, chi-squared computation
│   └── visualization.py      # Matplotlib figures
└── outputs/
    ├── {drug}_signals.csv    # Full results table
    └── figures/
        ├── {drug}_prr_bar.png        # PRR ± 95% CI bar chart
        ├── {drug}_scatter.png        # PRR vs ROR scatter
        ├── {drug}_heatmap.png        # Signal metrics heatmap
        └── {drug}_dashboard.png      # 2×2 summary dashboard
```

---

## Setup

```bash
git clone https://github.com/<your-username>/faers-signal-detection.git
cd faers-signal-detection
pip install -r requirements.txt
```

An openFDA API key is **optional** but recommended (raises rate limit from 40 → 240 requests/minute). Register free at https://open.fda.gov/apis/authentication/.

---

## Usage

```bash
# Analyse metformin (default)
python main.py

# Analyse a different drug
python main.py --drug warfarin

# With API key and custom AE count
python main.py --drug aspirin --top 30 --api_key YOUR_KEY
```

---

## Example Output

```
======================================================
  Signal Detection Summary
======================================================
  AEs analysed     : 25
  Signals detected : 9  (36.0%)
  Threshold        : PRR≥2, χ²≥4, n≥3
======================================================
ae_name                    count     PRR      ROR     chi2
LACTIC ACIDOSIS             1823    8.421    9.104  4821.3
METFORMIN OVERDOSE           934    6.230    6.401  2103.5
DIABETIC KETOACIDOSIS       2104    3.112    3.219  1834.2
...
```

---

## Visualisations

| Plot | Description |
|------|-------------|
| **PRR Bar Chart** | Top 20 AEs ranked by PRR with 95% CI; signals highlighted in red |
| **PRR vs ROR Scatter** | Concordance between both measures; bubble size = report count |
| **Heatmap** | Colour-coded table of all three metrics for top 15 AEs |
| **Dashboard** | Combined 2×2 overview for presentations / regulatory reports |

---

## Regulatory Context

This project directly reflects the type of analytical work performed within FDA's **Center for Drug Evaluation and Research (CDER)** as part of real-world evidence (RWE) monitoring programs. Disproportionality analysis is the primary quantitative method used to:

- Screen FAERS reports for emerging drug safety signals
- Support post-marketing safety surveillance requirements
- Inform regulatory decisions on label updates, Risk Evaluation and Mitigation Strategies (REMS), and market withdrawals

---

## References

- Evans, S.J.W. et al. (2001). Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports. *Pharmacoepidemiology and Drug Safety*, 10(6), 483–486.
- van Puijenbroek, E.P. et al. (2002). A comparison of measures of disproportionality for signal detection in spontaneous reporting systems. *Pharmacoepidemiology and Drug Safety*, 11(1), 3–10.
- FDA openFDA API: https://open.fda.gov/apis/drug/event/

---

## Author

Md Junayed  Nayeen Ph.D.

Master of Data Science(Ongoing), University of Pittsburgh  
