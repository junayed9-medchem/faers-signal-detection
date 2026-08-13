# FAERS Pharmacovigilance Signal Detection

**Post-marketing drug safety signal detection using the FDA Adverse Event Reporting System (FAERS), via the openFDA API.**

Md Junayed Nayeen, Ph.D.; Master of Data Science, University of Pittsburgh

---

## Background

The **FDA Adverse Event Reporting System (FAERS)** is the FDA's giant database of adverse event reports. Every time a doctor, patient, or drug company reports a side effect, it lands here. This project aims to spot which drug–side-effect pairs show up more often than expected by random chance. When that happens, it is noted as safety signal: a red flag worth investigating further, not proof of harm. Signal detection in FAERS is a cornerstone of the FDA's post-marketing pharmacovigilance program.

This notebook implements a disproportionality analysis workflow. The standard method regulatory agencies like the FDA use to scan spontaneous adverse event databases for potential drug safety signals without requiring a controlled clinical trial. Rather than asking whether an adverse event was ever reported with a drug (which is almost always true, given baseline disease rates in the population), the workflow asks whether it is reported disproportionately more often with that drug than expected from its background rate across all other drugs. This is done by building a (2×2) contingency table for each drug–event pair and computing PRR, ROR, and chi-squared, then flagging pairs that clear the Evans (2001) signal threshold. All data used in this analysis is pulled live via the openFDA Drug Adverse Event API, a free, public interface to the FDA Adverse Event Reporting System (FAERS) — the same underlying spontaneous-reporting database used in FDA's own post-marketing surveillance. This approach was chosen because it mirrors the exact methodology used in CDER's real-world evidence (RWE) monitoring programs, and is fully reproducible by anyone using only free, publicly available data.

> If a drug–adverse event combination is reported *more frequently than expected by chance*, it is flagged as a **safety signal** warranting further regulatory review.

The following cells are **in order**  and each one builds on the previous. A brief explanation is included in every step for comprehension.
