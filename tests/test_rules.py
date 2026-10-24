import pandas as pd

from kamr.analysis.concordance import evaluation_labels, metrics, wilson
from kamr.analysis.rules import amrfinder_markers, kleborate_markers, predictions


def amr_row(gid, symbol, cls, sub, subtype="AMR", etype="AMR"):
    return {
        "genome_id": gid,
        "Element symbol": symbol,
        "Type": etype,
        "Subtype": subtype,
        "Class": cls,
        "Subclass": sub,
    }


AMR = pd.DataFrame(
    [
        amr_row("1", "blaKPC-2", "BETA-LACTAM", "CARBAPENEM"),
        amr_row("1", "ompK35_G134AfsTer11", "BETA-LACTAM", "CARBAPENEM", subtype="POINT"),
        amr_row("1", "gyrA_S83I", "QUINOLONE", "QUINOLONE", subtype="POINT"),
        amr_row("1", "oqxA", "PHENICOL/QUINOLONE", "PHENICOL/QUINOLONE"),
        amr_row("1", "rmtB1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE"),
        amr_row("2", "ompK36_D135DGD", "BETA-LACTAM", "CARBAPENEM", subtype="POINT"),
        amr_row("2", "oqxB", "PHENICOL/QUINOLONE", "PHENICOL/QUINOLONE"),
        amr_row("2", "aac(3)-IId", "AMINOGLYCOSIDE", "GENTAMICIN"),
        amr_row("2", "mgrB_Q30STOP", "COLISTIN", "COLISTIN", subtype="POINT_DISRUPT"),
        amr_row("2", "fieF", "EFFLUX", "EFFLUX", etype="STRESS"),
    ]
)


def pred(markers, genomes, drug, gid):
    p = predictions(markers, genomes, "t")
    return p[(p.drug == drug) & (p.genome_id == gid)].iloc[0]


def test_amrfinder_rules():
    m = amrfinder_markers(AMR)
    g = ["1", "2", "3"]
    assert pred(m, g, "meropenem", "1").predicted == "R"
    assert pred(m, g, "meropenem", "1").evidence == "blaKPC-2"
    assert pred(m, g, "meropenem", "2").predicted == "S"  # porin mutation alone
    assert pred(m, g, "ciprofloxacin", "1").evidence == "gyrA_S83I"  # oqxA ignored
    assert pred(m, g, "ciprofloxacin", "2").predicted == "S"  # oqxB only
    assert pred(m, g, "gentamicin", "1").evidence == "rmtB1"
    assert pred(m, g, "gentamicin", "2").evidence == "aac(3)-IId"
    assert pred(m, g, "colistin", "2").evidence == "mgrB_Q30STOP"
    assert pred(m, g, "colistin", "3").predicted == "S"


def test_kleborate_rules():
    kl = pd.DataFrame(
        {
            "genome_id": ["1", "2"],
            "klebsiella_pneumo_complex__amr__Bla_Carb_acquired": ["KPC-2", "-"],
            "klebsiella_pneumo_complex__amr__Flq_mutations": ["GyrA:p.Ser83Ile", "-"],
            "klebsiella_pneumo_complex__amr__Flq_acquired": ["-", "QnrB1"],
            "klebsiella_pneumo_complex__amr__AGly_acquired": ["strA.v1^;aadA2^", "aac(3)-IId;strB"],
            "klebsiella_pneumo_complex__amr__Col_mutations": ["pmrb:p.Gly20fs", "-"],
        }
    )
    m = kleborate_markers(kl)
    g = ["1", "2"]
    assert pred(m, g, "meropenem", "1").predicted == "R"
    assert pred(m, g, "ciprofloxacin", "2").evidence == "QnrB1"
    assert pred(m, g, "gentamicin", "1").predicted == "S"  # streptomycin genes only
    assert pred(m, g, "gentamicin", "2").evidence == "aac(3)-IId"
    assert pred(m, g, "colistin", "1").predicted == "R"
    assert pred(m, g, "tigecycline", "1").predicted == "S"  # column absent


def test_metrics_and_wilson():
    df = pd.DataFrame({"label": list("RRRSS"), "predicted": list("RRSSR")})
    m = metrics(df)
    assert (m["TP"], m["FN"], m["FP"], m["TN"]) == (2, 1, 1, 1)
    assert abs(m["sensitivity"] - 2 / 3) < 1e-9
    lo, hi = wilson(2, 3)
    assert 0 < lo < 2 / 3 < hi < 1


def test_colistin_primary_set_is_broth_microdilution_only():
    labels = pd.DataFrame(
        {"genome_id": ["1", "2", "3"], "antibiotic": "colistin", "label": ["R", "S", "R"]}
    )
    records = pd.DataFrame(
        {
            "genome_id": ["1", "2", "3"],
            "antibiotic": "colistin",
            "laboratory_typing_method": ["Broth dilution", "Disk diffusion", "Broth dilution"],
            "record_label": ["R", "S", "R"],
        }
    )
    ev = evaluation_labels(labels, records, {"1", "2"})
    assert set(ev[ev.drug == "colistin"].genome_id) == {"1"}
    assert set(ev[ev.drug == "colistin_all"].genome_id) == {"1", "2"}
