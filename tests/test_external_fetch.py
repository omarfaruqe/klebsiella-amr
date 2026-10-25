import hashlib
import io
import zipfile

import pytest

from kamr.external.fetch_assemblies import download_one, expected_md5, extract_fasta

ACC = "GCA_000000001.1"
GENOME = ">contig1\n" + "ACGT" * 1_250_000 + "\n"  # 5 Mb, inside the accepted window


def archive(acc=ACC, body=GENOME, md5=None, member=None, extra=()):
    member = member or f"ncbi_dataset/data/{acc}/{acc}_ASM1v1_genomic.fna"
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr(member, body)
        digest = md5 if md5 is not None else hashlib.md5(body.encode()).hexdigest()
        if md5 != "omit":
            z.writestr("md5sum.txt", f"{digest}  {member}\n")
        for name in extra:
            z.writestr(name, body)
    return buf.getvalue()


def test_extract_verifies_md5():
    text, status = extract_fasta(archive(), ACC)
    assert text.startswith(">contig1")
    assert status == "md5 verified"


def test_extract_rejects_a_corrupted_payload():
    # a truncated transfer must fail loudly, not yield a short genome that QC would later
    # misread as a fragmented assembly
    with pytest.raises(ValueError, match="md5 mismatch"):
        extract_fasta(archive(md5="0" * 32), ACC)


def test_extract_tolerates_a_missing_md5_listing():
    _, status = extract_fasta(archive(md5="omit"), ACC)
    assert status == "md5 not listed"


def test_extract_rejects_an_ambiguous_archive():
    extra = (f"ncbi_dataset/data/{ACC}/{ACC}_ASM2v2_genomic.fna",)
    with pytest.raises(ValueError, match="expected one genomic FASTA"):
        extract_fasta(archive(extra=extra), ACC)


def test_extract_ignores_another_accessions_directory():
    other = "ncbi_dataset/data/GCA_999.9/GCA_999.9_x_genomic.fna"
    with pytest.raises(ValueError, match="expected one genomic FASTA"):
        extract_fasta(archive(member=other), ACC)


def test_expected_md5_matches_despite_leading_dot_slash():
    buf = io.BytesIO()
    member = f"ncbi_dataset/data/{ACC}/{ACC}_ASM1v1_genomic.fna"
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr(member, GENOME)
        z.writestr("md5sum.txt", f"abc123  ./{member}\n")
    with zipfile.ZipFile(io.BytesIO(buf.getvalue())) as z:
        assert expected_md5(z, member) == "abc123"


def test_download_one_is_idempotent(tmp_path):
    # an existing file is reported without a network call, so an interrupted run can restart
    (tmp_path / f"{ACC}.fna").write_text(GENOME)
    res = download_one(ACC, tmp_path)
    assert res["status"] == "already present"
    assert res["length"] == 5_000_000
    assert res["sha256"] == hashlib.sha256(GENOME.encode()).hexdigest()


def test_download_one_rejects_an_implausible_length(tmp_path, monkeypatch):
    small = ">c\n" + "ACGT" * 10 + "\n"

    class Resp:
        status_code = 200
        content = archive(body=small)

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "kamr.external.fetch_assemblies.session",
        lambda: type("S", (), {"get": lambda *a, **k: Resp()})(),
    )
    res = download_one(ACC, tmp_path)
    assert "outside" in res["status"]
    assert not (tmp_path / f"{ACC}.fna").exists()
