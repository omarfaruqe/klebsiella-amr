import pandas as pd

from kamr.analysis.mgrb_flanks import (
    best_fragments,
    distance_to_contig_end,
    flank_sequence,
    read_fasta,
)


def frag(lo, hi, slen, length=40):
    return pd.Series(
        {
            "sseqid": "c1",
            "sstart": lo,
            "send": hi,
            "lo": lo,
            "hi": hi,
            "slen": slen,
            "length": length,
        }
    )


def test_read_fasta(tmp_path):
    p = tmp_path / "a.fna"
    p.write_text(">c1 desc here\nACGT\nACGT\n>c2\nTTTT\n")
    assert read_fasta(p) == {"c1": "ACGTACGT", "c2": "TTTT"}


def test_best_fragments_normalises_orientation():
    hits = pd.DataFrame(
        [
            {"sseqid": "c1", "sstart": 900, "send": 800, "length": 40, "slen": 5000},
            {"sseqid": "c1", "sstart": 10, "send": 20, "length": 5, "slen": 5000},
        ]
    )
    f = best_fragments(hits)
    assert len(f) == 1 and f.iloc[0]["lo"] == 800 and f.iloc[0]["hi"] == 900


def test_distance_to_contig_end_picks_nearer_edge():
    assert distance_to_contig_end(frag(51, 150, 5000)) == 50
    assert distance_to_contig_end(frag(4800, 4950, 5000)) == 50
    assert distance_to_contig_end(frag(2000, 2100, 5000)) == 1999


def test_flank_sequence_takes_the_nearer_end():
    seq = "A" * 100 + "G" * 50 + "C" * 2000
    side, fl = flank_sequence(seq, frag(101, 150, len(seq)), flank_bp=60)
    assert side == "left" and fl == "A" * 60
    side, fl = flank_sequence(seq, frag(101, 150, 200), flank_bp=60)
    assert side == "right" and set(fl) == {"C"}
