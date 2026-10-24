"""Rule-based resistance predictions per genome and drug, from AMRFinderPlus or Kleborate output.

Each rule returns a table with one row per genome and drug: predicted ("R" or "S") and the
evidence (markers that triggered the call). Rules are deliberately simple and documented, so the
concordance analysis measures what known mechanisms explain, not a tuned model.

AMRFinderPlus rules (Type == AMR only). For every drug, the chromosomal oqxA/oqxB efflux genes
(intrinsic in K. pneumoniae) are ignored. Efflux-regulator variants (oqxR, ramR, acrR, marR,
soxR/soxS, ramA) are ignored by the core rule and included only by the "+regulators" rule,
because many are lineage polymorphisms rather than resistance mechanisms.
  meropenem      acquired gene with Subclass CARBAPENEM (carbapenemases; porin mutations excluded)
  ciprofloxacin  quinolone target mutation (gyrA/parC/parE) or acquired quinolone gene
  gentamicin     Subclass GENTAMICIN, or a 16S rRNA methyltransferase (armA, rmt*, npmA)
  colistin       any element with Class or Subclass COLISTIN (mcr genes, mgrB/pmrAB/phoPQ variants)
  tigecycline    Subclass TIGECYCLINE (core: acquired tet(X), tmexCD-toprJ)

Kleborate rules (kpsc AMR columns):
  meropenem      Bla_Carb_acquired
  ciprofloxacin  Flq_acquired or Flq_mutations
  gentamicin     AGly_acquired genes aac(3), ant(2''), aadB, armA, rmt*, npmA
  colistin       Col_acquired or Col_mutations
  tigecycline    Tgc_acquired or Tgc_mutations
"""

from __future__ import annotations

import re

import pandas as pd

DRUGS = ("meropenem", "ciprofloxacin", "gentamicin", "colistin", "tigecycline")
RMTASE = re.compile(r"^(armA|rmt|npmA)", re.IGNORECASE)
INTRINSIC = re.compile(r"^oqx[AB](\d|$)", re.IGNORECASE)
REGULATOR = re.compile(r"^(oqxR|ramR|ramA|acrR|marR|soxR|soxS)(_|$)", re.IGNORECASE)
KLEB_GENT = re.compile(r"^(aac\(3\)|aacC|ant\(2''\)|aadB|armA|rmt|npmA)", re.IGNORECASE)
KLEB_COLUMNS = {
    "meropenem": ("Bla_Carb_acquired",),
    "ciprofloxacin": ("Flq_acquired", "Flq_mutations"),
    "gentamicin": ("AGly_acquired",),
    "colistin": ("Col_acquired", "Col_mutations"),
    "tigecycline": ("Tgc_acquired", "Tgc_mutations"),
}


def _col(df: pd.DataFrame, *names: str) -> pd.Series:
    for name in names:
        if name in df.columns:
            return df[name].fillna("").astype(str)
    return pd.Series("", index=df.index)


def amrfinder_markers(amr: pd.DataFrame, include_regulators: bool = False) -> pd.DataFrame:
    """Return long table genome_id, drug, marker for every AMRFinderPlus hit a rule uses."""
    symbol = _col(amr, "Element symbol", "Gene symbol")
    etype = _col(amr, "Type", "Element type").str.upper()
    subtype = _col(amr, "Subtype", "Element subtype").str.upper()
    cls = _col(amr, "Class").str.upper()
    sub = _col(amr, "Subclass").str.upper()
    is_amr = (etype == "AMR") & ~symbol.str.match(INTRINSIC)
    if not include_regulators:
        is_amr &= ~symbol.str.match(REGULATOR)
    point = subtype.str.startswith("POINT")
    masks = {
        "meropenem": is_amr & ~point & sub.str.contains("CARBAPENEM"),
        "ciprofloxacin": is_amr & (cls.str.contains("QUINOLONE") | sub.str.contains("QUINOLONE")),
        "gentamicin": is_amr & (sub.str.contains("GENTAMICIN") | symbol.str.match(RMTASE)),
        "colistin": is_amr & (cls.str.contains("COLISTIN") | sub.str.contains("COLISTIN")),
        "tigecycline": is_amr & sub.str.contains("TIGECYCLINE"),
    }
    frames = [
        pd.DataFrame({"genome_id": amr.loc[m, "genome_id"], "drug": d, "marker": symbol[m]})
        for d, m in masks.items()
    ]
    return pd.concat(frames, ignore_index=True)


def kleborate_markers(kl: pd.DataFrame) -> pd.DataFrame:
    """Return long table genome_id, drug, marker from Kleborate AMR columns."""
    rows = []
    for drug, suffixes in KLEB_COLUMNS.items():
        cols = [c for c in kl.columns if any(c.endswith(f"__amr__{s}") for s in suffixes)]
        for col in cols:
            for gid, value in zip(kl["genome_id"], kl[col].fillna("-").astype(str), strict=True):
                if value.strip() in ("", "-"):
                    continue
                for token in value.split(";"):
                    token = token.strip()
                    if not token:
                        continue
                    if drug == "gentamicin" and not KLEB_GENT.match(token):
                        continue
                    rows.append((gid, drug, token))
    return pd.DataFrame(rows, columns=["genome_id", "drug", "marker"])


def predictions(markers: pd.DataFrame, genomes: list[str], tool: str) -> pd.DataFrame:
    """One row per genome x drug: predicted R if the rule found any marker, else S."""
    grid = pd.MultiIndex.from_product([sorted(genomes), DRUGS], names=["genome_id", "drug"])
    evidence = markers.groupby(["genome_id", "drug"])["marker"].apply(
        lambda s: ";".join(sorted(set(s)))
    )
    out = pd.DataFrame(index=grid).join(evidence.rename("evidence")).reset_index()
    out["evidence"] = out["evidence"].fillna("")
    out["predicted"] = out["evidence"].map(lambda e: "R" if e else "S")
    out["tool"] = tool
    return out
