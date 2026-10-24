import pandas as pd

from kamr.analysis.review_gentamicin import aminoglycoside_genes, enrichment, highest_mic


def test_highest_mic_bins():
    rec = pd.DataFrame(
        {
            "genome_id": ["1", "1", "2"],
            "antibiotic": "gentamicin",
            "eucast_call": ["R", "R", "S"],
            "measurement_value": ["4", "8", "1"],
            "laboratory_typing_method": ["Broth dilution", "Agar dilution", None],
        }
    )
    m = highest_mic(rec, "gentamicin")
    assert m.loc["1", "mic"] == 8 and m.loc["1", "mic_bin"] == "8"
    assert m.loc["1", "method"] == "Agar dilution"
    assert m.loc["2", "mic_bin"] == "1" and m.loc["2", "method"] == "not reported"


def test_aminoglycoside_genes_and_enrichment():
    amr = pd.DataFrame(
        {
            "genome_id": ["1", "1", "2", "3", "3"],
            "Element symbol": ["aac(6')-Ib", "oqxA", "aac(6')-Ib", "aph(3')-Ia", "blaKPC-2"],
            "Type": "AMR",
            "Class": [
                "AMINOGLYCOSIDE",
                "AMINOGLYCOSIDE",
                "AMINOGLYCOSIDE",
                "AMINOGLYCOSIDE",
                "BETA-LACTAM",
            ],
        }
    )
    genes = aminoglycoside_genes(amr)
    assert set(genes["gene"]) == {"aac(6')-Ib", "aph(3')-Ia"}
    e = enrichment(genes, fn={"1", "2"}, tn={"3"})
    assert e.loc["aac(6')-Ib", "in_FN"] == 2
    assert e.loc["aac(6')-Ib", "FN_share"] == 1.0
    assert e.loc["aph(3')-Ia", "in_TN"] == 1
