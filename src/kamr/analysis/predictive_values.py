"""Translate sensitivity and specificity into predictive values across plausible prevalences.

Sensitivity and specificity do not depend on how common resistance is; positive and negative
predictive value do. Our cohort is a surveillance collection in which 24% of colistin-tested
genomes are resistant, which is far above most clinical settings, so the rule's discrimination
says less about its clinical usefulness than its predictive values at a realistic prevalence do.

This module exists so the predictive values quoted in the limitations are derived by code from the
reported sensitivity and specificity rather than asserted.

Outputs: reports/review/predictive_values.md
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

# Illustrative prevalences: a low-endemicity clinical setting, a moderate one, a high-burden
# setting, and our own cohort. These bracket reported ranges; they are not estimates for any
# particular institution.
PREVALENCES = (0.02, 0.05, 0.10, 0.24)


def predictive_values(
    sensitivity: float, specificity: float, prevalence: float
) -> tuple[float, float]:
    """Return (PPV, NPV) by Bayes' theorem. Raises on inputs outside [0, 1]."""
    for name, v in (
        ("sensitivity", sensitivity),
        ("specificity", specificity),
        ("prevalence", prevalence),
    ):
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1], got {v}")
    tp = sensitivity * prevalence
    fp = (1.0 - specificity) * (1.0 - prevalence)
    tn = specificity * (1.0 - prevalence)
    fn = (1.0 - sensitivity) * prevalence
    ppv = tp / (tp + fp) if tp + fp else float("nan")
    npv = tn / (tn + fn) if tn + fn else float("nan")
    return ppv, npv


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sensitivity", type=float, default=0.567, help="internal headline rule")
    p.add_argument("--specificity", type=float, default=0.952)
    p.add_argument("--label", default="AMRFinderPlus + mgrB disrupted (internal, §3.3)")
    p.add_argument("--out", type=Path, default=Path("reports/review/predictive_values.md"))
    a = p.parse_args()

    lines = [
        "# Predictive values across plausible resistance prevalences",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Rule: {a.label}",
        f"- Sensitivity {a.sensitivity:.1%}, specificity {a.specificity:.1%}",
        "",
        "Sensitivity and specificity are properties of the rule; predictive values also depend on",
        "how common resistance is in the population tested. Prevalences below are illustrative.",
        "",
        "| Colistin resistance prevalence | PPV | NPV |",
        "|---|---|---|",
    ]
    for prev in PREVALENCES:
        ppv, npv = predictive_values(a.sensitivity, a.specificity, prev)
        note = " (this cohort)" if abs(prev - 0.24) < 1e-9 else ""
        lines.append(f"| {prev:.0%}{note} | {ppv:.1%} | {npv:.1%} |")
    text = "\n".join(lines) + "\n"
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
