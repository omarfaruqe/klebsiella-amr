import pandas as pd

from kamr.external.validate import evaluate, table


def tables():
    g = [f"GCA_{i}.1" for i in range(10)]
    # a rule that catches the first three resistant genomes and nothing else
    return g, {"rule A": pd.Series(["R"] * 3 + ["S"] * 7, index=g)}


def labels(g):
    return pd.DataFrame(
        {
            "genome_id": g,
            "label": ["R"] * 5 + ["S"] * 5,
            "label_source": ["mic_eucast"] * 8 + ["submitted"] * 2,
        }
    )


def test_evaluate_reports_sensitivity_specificity_and_counts():
    g, tab = tables()
    res = evaluate(labels(g), tab)
    row = res.iloc[0]
    assert row["n"] == 10
    assert row["sensitivity"] == 3 / 5  # 3 of 5 resistant caught
    assert row["specificity"] == 1.0  # no susceptible called resistant
    assert row["FN"] == 2
    assert row["FP"] == 0


def test_evaluate_drops_genomes_the_rule_does_not_cover():
    # a genome absent from the rule's series must be excluded, not counted as susceptible,
    # which would silently inflate specificity
    g, tab = tables()
    tab["rule A"] = tab["rule A"].drop(g[0])
    res = evaluate(labels(g), tab)
    assert res.iloc[0]["n"] == 9


def test_evaluate_respects_a_label_stratum():
    g, tab = tables()
    lab = labels(g)
    res = evaluate(lab[lab["label_source"] == "mic_eucast"], tab)
    assert res.iloc[0]["n"] == 8


def test_table_renders_intervals_per_rule():
    g, tab = tables()
    lines = table(evaluate(labels(g), tab))
    assert lines[0].startswith("| Rule |")
    assert len([ln for ln in lines if ln.startswith("| rule A")]) == 1


def test_composition_separates_determinant_from_mgrb_and_unexplained():
    from kamr.external.validate import composition

    g = [f"GCA_{i}.1" for i in range(4)]
    lab = pd.DataFrame({"genome_id": g, "label": ["R"] * 4})
    # g0 has mcr only, g1 mgrB only, g2 both, g3 neither
    amr = pd.DataFrame(
        {
            "genome_id": [g[0], g[2]],
            "Element symbol": ["mcr-1.1", "mcr-9"],
            "Type": ["AMR", "AMR"],
            "Subtype": ["AMR", "AMR"],
            "Class": ["COLISTIN", "COLISTIN"],
            "Subclass": ["COLISTIN", "COLISTIN"],
        }
    )
    status = pd.DataFrame(
        {
            "genome_id": g,
            "mgrb_status": ["intact", "absent", "absent", "intact"],
            "mgrb_disrupted": [False, True, True, False],
        }
    )
    counts, genes = composition(lab, amr, status)
    assert counts["acquired determinant"] == 1
    assert counts["mgrB disrupted"] == 1
    assert counts["both"] == 1
    assert counts["unexplained"] == 1
    assert set(genes.index) == {"mcr-1.1", "mcr-9"}


def test_composition_ignores_susceptible_isolates():
    from kamr.external.validate import composition

    g = ["GCA_0.1", "GCA_1.1"]
    lab = pd.DataFrame({"genome_id": g, "label": ["R", "S"]})
    amr = pd.DataFrame(
        {
            "genome_id": [g[1]],
            "Element symbol": ["mcr-1.1"],
            "Type": ["AMR"],
            "Subtype": ["AMR"],
            "Class": ["COLISTIN"],
            "Subclass": ["COLISTIN"],
        }
    )
    status = pd.DataFrame(
        {"genome_id": g, "mgrb_status": ["intact"] * 2, "mgrb_disrupted": [False, False]}
    )
    counts, genes = composition(lab, amr, status)
    assert counts.sum() == 1
    assert counts["unexplained"] == 1
    assert genes.empty
