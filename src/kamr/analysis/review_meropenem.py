"""Sprint 2 review: why meropenem rules disagree with laboratory phenotype.

Question 1 (false positives): of genomes carrying a carbapenemase, which enzyme families still
test susceptible, and how close are their MICs to the EUCAST breakpoint (S <= 2 mg/L)?
Question 2 (false negatives): resistant genomes without a carbapenemase - do they carry an ESBL or
AmpC plus outer-membrane porin (OmpK35/OmpK36) defects, the classic non-carbapenemase route?

Inputs:  results/cohort/amrfinder.parquet, data/curated/cohort_labels.parquet,
         data/curated/ast_records.parquet, data/curated/genome_qc.parquet,
         reports/cohort/qc_pass_genomes.tsv
Outputs: reports/review/meropenem.md, reports/review/meropenem_carbapenemase.csv,
         reports/review/meropenem_genomes.tsv
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
from pathlib import Path

import pandas as pd

from kamr.analysis.rules import INTRINSIC

OXA48_LIKE = re.compile(
    r"^blaOXA-(48|162|163|181|204|232|244|245|370|405|436|484|505|515|517|519|535|546)(\D|$)"
)
FAMILIES = (
    ("KPC", re.compile(r"^blaKPC")),
    ("NDM", re.compile(r"^blaNDM")),
    ("OXA-48-like", OXA48_LIKE),
    ("VIM", re.compile(r"^blaVIM")),
    ("IMP", re.compile(r"^blaIMP")),
)
ESBL_AMPC = re.compile(
    r"^bla(CTX-M|DHA|CMY|ACT|FOX|MOX|ACC|SHV-(2|5|12|27|31|38)|TEM-(10|12|26|52))"
)
PORIN = re.compile(r"^ompK3[56]", re.IGNORECASE)
MIC_BINS = [0, 0.125, 0.25, 0.5, 1, 2, 8, 1e9]
MIC_LABELS = ["<=0.125", "0.25", "0.5", "1", "2", "4-8", ">8"]


def carbapenemase_family(symbol: str) -> str:
    for name, pattern in FAMILIES:
        if pattern.match(symbol):
            return name
    return "other"


def genome_mechanisms(amr: pd.DataFrame) -> pd.DataFrame:
    """Per genome: carbapenemase families and alleles, ESBL/AmpC genes, porin defects."""
    sym = amr["Element symbol"].fillna("").astype(str)
    etype = amr["Type"].fillna("").astype(str).str.upper()
    sub = amr["Subclass"].fillna("").astype(str).str.upper()
    subtype = amr["Subtype"].fillna("").astype(str).str.upper()
    is_amr = (etype == "AMR") & ~sym.str.match(INTRINSIC)
    carb = is_amr & ~subtype.str.startswith("POINT") & sub.str.contains("CARBAPENEM")
    rows = amr.assign(symbol=sym)

    def collect(mask: pd.Series, col: str, fn=None) -> pd.Series:
        vals = rows.loc[mask, ["genome_id", "symbol"]].copy()
        if fn is not None:
            vals["symbol"] = vals["symbol"].map(fn)
        return (
            vals.groupby("genome_id")["symbol"]
            .apply(lambda s: ";".join(sorted(set(s))))
            .rename(col)
        )

    out = pd.concat(
        [
            collect(carb, "carbapenemase_alleles"),
            collect(carb, "carbapenemase_family", carbapenemase_family),
            collect(is_amr & sym.str.match(ESBL_AMPC), "esbl_ampc"),
            collect(sym.str.match(PORIN) & subtype.str.startswith("POINT"), "porin_defects"),
        ],
        axis=1,
    )
    return out.fillna("")


def meropenem_mic(records: pd.DataFrame) -> pd.DataFrame:
    """Highest re-interpretable meropenem MIC per genome, with its sign, for laboratory records."""
    mero = records[(records["antibiotic"] == "meropenem") & records["eucast_call"].notna()].copy()
    mero["mic"] = pd.to_numeric(mero["measurement_value"], errors="coerce")
    mero = mero.dropna(subset=["mic"]).sort_values("mic")
    top = mero.groupby("genome_id").tail(1).set_index("genome_id")
    sign = top["measurement_sign"].fillna("").astype(str).str.strip()
    sign = sign.where(~sign.isin(["=", "=="]), "")
    return pd.DataFrame(
        {
            "mic": top["mic"],
            "mic_text": sign.where(sign.str.len() > 0, "") + top["mic"].map(lambda v: f"{v:g}"),
            "mic_bin": pd.cut(top["mic"], MIC_BINS, labels=MIC_LABELS, right=True),
            "method": top["laboratory_typing_method"].fillna("not reported"),
            "reported_standard": top["testing_standard"].fillna("not reported"),
        }
    )


def build_table(amr, labels, records, qc, genomes) -> pd.DataFrame:
    lab = labels[(labels["antibiotic"] == "meropenem") & labels["label"].isin(["R", "S"])]
    lab = lab[lab["genome_id"].isin(genomes)].set_index("genome_id")[["label"]]
    mech = genome_mechanisms(amr[amr["genome_id"].isin(genomes)])
    meta = qc.set_index("genome_id")[
        ["st", "bioproject_accession", "collection_year", "isolation_country"]
    ]
    table = (
        lab.join(mech, how="left").join(meropenem_mic(records), how="left").join(meta, how="left")
    )
    for col in ("carbapenemase_alleles", "carbapenemase_family", "esbl_ampc", "porin_defects"):
        table[col] = table[col].fillna("")
    table["has_carbapenemase"] = table["carbapenemase_family"] != ""
    table["family_group"] = table["carbapenemase_family"].map(
        lambda f: "none" if not f else (f if ";" not in f else "multiple")
    )
    return table.reset_index()


def family_summary(table: pd.DataFrame) -> pd.DataFrame:
    carb = table[table["has_carbapenemase"]]
    g = carb.groupby("family_group")["label"].value_counts().unstack(fill_value=0)
    for col in ("R", "S"):
        if col not in g:
            g[col] = 0
    g["genomes"] = g["R"] + g["S"]
    g["S_share"] = (g["S"] / g["genomes"]).round(3)
    return g[["genomes", "R", "S", "S_share"]].sort_values("genomes", ascending=False)


def md_table(df: pd.DataFrame, index_name: str) -> list[str]:
    """Markdown table; float columns are shown as percentages, integer columns as counts."""
    cols = [index_name, *map(str, df.columns)]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    is_share = {c: pd.api.types.is_float_dtype(df[c]) for c in df.columns}
    for idx, row in df.iterrows():
        cells = [
            f"{row[c]:.1%}"
            if is_share[c]
            else f"{int(row[c]):,}"
            if pd.api.types.is_number(row[c])
            else str(row[c])
            for c in df.columns
        ]
        lines.append("| " + " | ".join([str(idx), *cells]) + " |")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinder", type=Path, default=Path("results/cohort/amrfinder.parquet"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--genomes", type=Path, default=Path("reports/cohort/qc_pass_genomes.tsv"))
    parser.add_argument("--out", type=Path, default=Path("reports/review"))
    args = parser.parse_args()

    genomes = set(pd.read_csv(args.genomes, sep="\t", dtype=str)["genome_id"])
    table = build_table(
        pd.read_parquet(args.amrfinder),
        pd.read_parquet(args.curated / "cohort_labels.parquet"),
        pd.read_parquet(args.curated / "ast_records.parquet"),
        pd.read_parquet(args.curated / "genome_qc.parquet"),
        genomes,
    )
    fams = family_summary(table)

    fp = table[table["has_carbapenemase"] & (table["label"] == "S")]
    tp = table[table["has_carbapenemase"] & (table["label"] == "R")]
    fn = table[~table["has_carbapenemase"] & (table["label"] == "R")]
    tn = table[~table["has_carbapenemase"] & (table["label"] == "S")]

    mic_fp = pd.crosstab(fp["family_group"], fp["mic_bin"].astype(str)).reindex(
        columns=[c for c in MIC_LABELS if c in set(fp["mic_bin"].astype(str))], fill_value=0
    )
    near = int(fp["mic_bin"].astype(str).isin(["1", "2"]).sum())

    def route(df: pd.DataFrame) -> pd.Series:
        esbl = df["esbl_ampc"] != ""
        porin = df["porin_defects"] != ""
        return pd.Series(
            {
                "ESBL/AmpC + porin defect": int((esbl & porin).sum()),
                "ESBL/AmpC only": int((esbl & ~porin).sum()),
                "porin defect only": int((~esbl & porin).sum()),
                "neither": int((~esbl & ~porin).sum()),
            }
        )

    routes = pd.DataFrame({"R, no carbapenemase": route(fn), "S, no carbapenemase": route(tn)})

    lines = [
        "# Meropenem discordance review",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes with a meropenem S/R label: {len(table):,} "
        f"(R {int((table.label == 'R').sum()):,}, S {int((table.label == 'S').sum()):,})",
        "- Carbapenemase calls: AMRFinderPlus acquired genes with Subclass CARBAPENEM",
        "- EUCAST v16.1 meropenem: S <= 2 mg/L, R > 8 mg/L",
        "",
        f"Carbapenemase present: {len(tp) + len(fp):,} genomes, of which {len(fp):,} test S "
        "(rule false positives).",
        f"No carbapenemase: {len(fn) + len(tn):,} genomes, of which {len(fn):,} test R "
        "(rule false negatives).",
        "",
        "## 1. Carbapenemase family vs phenotype",
        "",
        *md_table(fams, "Family"),
        "",
        "## 2. MIC of susceptible carbapenemase carriers (mg/L)",
        "",
        f"{near:,} of {len(fp):,} false positives sit at 1-2 mg/L, within two doubling dilutions "
        "of the breakpoint.",
        "",
        *md_table(mic_fp, "Family"),
        "",
        "## 3. Mechanism routes without a carbapenemase",
        "",
        "Porin defects = AMRFinderPlus point/disruption calls in ompK35 or ompK36.",
        "",
        *md_table(routes, "Route"),
        "",
        "## 4. Where the false positives come from",
        "",
        "Top BioProjects:",
        "",
        *md_table(
            fp["bioproject_accession"].fillna("missing").value_counts().head(8).to_frame("genomes"),
            "BioProject",
        ),
        "",
        "Top STs:",
        "",
        *md_table(fp["st"].fillna("unknown").value_counts().head(8).to_frame("genomes"), "ST"),
        "",
        "Per-genome detail: `reports/review/meropenem_genomes.tsv`.",
    ]
    args.out.mkdir(parents=True, exist_ok=True)
    fams.to_csv(args.out / "meropenem_carbapenemase.csv")
    keep = [
        "genome_id",
        "label",
        "carbapenemase_alleles",
        "carbapenemase_family",
        "mic_text",
        "method",
        "reported_standard",
        "esbl_ampc",
        "porin_defects",
        "st",
        "bioproject_accession",
        "collection_year",
        "isolation_country",
    ]
    discordant = table[table["has_carbapenemase"] != (table["label"] == "R")]
    discordant[keep].sort_values(["label", "carbapenemase_family", "genome_id"]).to_csv(
        args.out / "meropenem_genomes.tsv", sep="\t", index=False
    )
    (args.out / "meropenem.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
