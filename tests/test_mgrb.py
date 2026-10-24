import pandas as pd

from kamr.analysis.mgrb import QLEN, classify

FULL = "MKKLRWVLLIVIIAGCLLLWTQMLNVMCDQDVQFFSGICTINKFIPW"


def hits(*rows):
    cols = [
        "sseqid",
        "qstart",
        "qend",
        "sstart",
        "send",
        "pident",
        "length",
        "evalue",
        "sseq",
        "slen",
    ]
    return pd.DataFrame(rows, columns=cols)


def test_absent_and_intact():
    assert classify(hits())[0] == "absent"
    assert classify(hits(("c1", 1, 47, 1000, 1141, 100.0, 47, 1e-30, FULL, 50000)))[0] == "intact"
    assert classify(hits(("c1", 1, 47, 1000, 1141, 97.9, 47, 1e-30, FULL, 50000)))[0] == "missense"


def test_premature_stop():
    sseq = FULL[:20] + "*" + FULL[21:]
    assert classify(hits(("c1", 1, 47, 1000, 1141, 97.0, 47, 1e-25, sseq, 50000)))[0] == (
        "premature_stop"
    )


def test_contig_break_and_partial():
    edge = hits(("c7", 1, 30, 9910, 9999, 100.0, 30, 1e-15, FULL[:30], 10000))
    assert classify(edge)[0] == "contig_break"
    mid = hits(("c7", 1, 30, 5000, 5089, 100.0, 30, 1e-15, FULL[:30], 10000))
    assert classify(mid)[0] == "partial"


def test_two_parts():
    same_far = hits(
        ("c1", 1, 25, 1000, 1074, 100.0, 25, 1e-12, FULL[:25], 50000),
        ("c1", 24, QLEN, 2400, 2471, 100.0, 24, 1e-10, FULL[23:], 50000),
    )
    assert classify(same_far)[0] == "insertion_assembled"
    same_near = hits(
        ("c1", 1, 25, 1000, 1074, 100.0, 25, 1e-12, FULL[:25], 50000),
        ("c1", 24, QLEN, 1080, 1151, 100.0, 24, 1e-10, FULL[23:], 50000),
    )
    assert classify(same_near)[0] == "frameshift_or_indel"
    split = hits(
        ("c1", 1, 25, 9926, 10000, 100.0, 25, 1e-12, FULL[:25], 10000),
        ("c2", 24, QLEN, 1, 72, 100.0, 24, 1e-10, FULL[23:], 8000),
    )
    assert classify(split)[0] == "split_contigs"


def test_short_spurious_hits_ignored():
    assert classify(hits(("c1", 5, 12, 100, 123, 90.0, 8, 1e-6, "LRWVLLIV", 50000)))[0] == "absent"
