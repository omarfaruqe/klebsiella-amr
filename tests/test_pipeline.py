import pandas as pd

from kamr.ingestion.assemblies import check_fasta, fasta_stats
from kamr.pipeline.collate import read_amrfinder, read_kleborate

FASTA = ">c1\nACGT\nACGT\n>c2\nAAAA\n"


def test_fasta_stats():
    assert fasta_stats(FASTA) == (2, 12)


def test_check_fasta():
    assert check_fasta(FASTA, 12) == "ok"
    assert check_fasta(FASTA, 12.4) == "ok"
    assert check_fasta(FASTA, 20) == "length_mismatch"
    assert check_fasta(FASTA, None) == "ok"
    assert check_fasta("<html>error</html>", 12) == "not_fasta"
    assert check_fasta(">c1\n", 12) == "empty"


def test_read_amrfinder(tmp_path):
    d = tmp_path / "amrfinder"
    d.mkdir()
    pd.DataFrame({"Name": ["573.1"], "Element symbol": ["blaKPC-2"]}).to_csv(
        d / "573.1.tsv", sep="\t", index=False
    )
    pd.DataFrame({"Name": ["573.2"] * 2, "Element symbol": ["fosA", "blaSHV-11"]}).to_csv(
        d / "573.2.tsv", sep="\t", index=False
    )
    (d / "573.3.tsv.tmp").write_text("partial")
    df = read_amrfinder(d)
    assert list(df.columns) == ["genome_id", "Element symbol"]
    assert df.groupby("genome_id").size().to_dict() == {"573.1": 1, "573.2": 2}


def test_read_kleborate_uses_full_table_and_done_marker(tmp_path):
    root = tmp_path / "kleborate"
    done = root / "573.1"
    done.mkdir(parents=True)
    (done / ".done").touch()
    pd.DataFrame({"strain": ["x"], "ST": ["ST11"]}).to_csv(
        done / "a_summary.txt", sep="\t", index=False
    )
    pd.DataFrame({"strain": ["x"], "klebsiella_pneumo_complex__mlst__ST": ["ST11"]}).to_csv(
        done / "b_full.txt", sep="\t", index=False
    )
    unfinished = root / "573.2"
    unfinished.mkdir()
    pd.DataFrame({"klebsiella_pneumo_complex__mlst__ST": ["ST258"]}).to_csv(
        unfinished / "full.txt", sep="\t", index=False
    )
    (root / ".573.3.tmp").mkdir()
    df = read_kleborate(root)
    assert df["genome_id"].tolist() == ["573.1"]
    assert "klebsiella_pneumo_complex__mlst__ST" in df.columns
