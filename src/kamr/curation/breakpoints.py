"""Load a pinned breakpoint table and interpret MIC measurements."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

MIC_UNITS = {"mg/l", "ug/ml", "mcg/ml"}


@dataclass(frozen=True)
class Breakpoint:
    antibiotic: str
    s_le: float
    r_gt: float
    applies: bool
    version: str


def load_breakpoints(path: Path) -> tuple[dict[str, Breakpoint], bool]:
    """Return breakpoints keyed by antibiotic, and whether every row is verified."""
    table = pd.read_csv(path, dtype=str).fillna("")
    bps = {}
    for row in table.itertuples(index=False):
        name = row.antibiotic.strip().lower()
        bps[name] = Breakpoint(
            antibiotic=name,
            s_le=float(row.s_le_mg_per_l),
            r_gt=float(row.r_gt_mg_per_l),
            applies=row.applies_to_kpneumoniae.strip().lower() == "yes",
            version=f"{row.standard} {row.version}",
        )
    verified = bool((table["verified"].str.strip().str.lower() == "yes").all())
    return bps, verified


def normalize_unit(unit: object) -> str:
    text = "" if unit is None or (isinstance(unit, float) and pd.isna(unit)) else str(unit)
    return text.strip().lower().replace("µ", "u").replace("μ", "u").replace(" ", "")


def interpret_mic(sign: object, value: object, bp: Breakpoint) -> str | None:
    """Return 'S', 'I' or 'R' for an MIC, or None if censoring makes the call undecidable.

    EUCAST convention: S if MIC <= s_le, R if MIC > r_gt, otherwise I
    ("susceptible, increased exposure").
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    s = "" if sign is None or (isinstance(sign, float) and pd.isna(sign)) else str(sign).strip()
    if s in ("", "=", "=="):
        if v <= bp.s_le:
            return "S"
        if v > bp.r_gt:
            return "R"
        return "I"
    if s in ("<=", "<"):
        return "S" if v <= bp.s_le else None
    if s == ">":
        return "R" if v >= bp.r_gt else None
    if s == ">=":
        return "R" if v > bp.r_gt else None
    return None
