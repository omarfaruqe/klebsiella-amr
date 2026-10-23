#!/usr/bin/env bash
# Verify the AMR toolchain on this machine against reference genome HS11286.
# Run inside the kamr environment. Prints versions, key calls and a fingerprint;
# matching fingerprints on two machines mean identical normalized results.
set -euo pipefail

ACC="GCF_000240185.1"
WORK="${KAMR_VERIFY_DIR:-$HOME/kamr-scratch/verify}"
THREADS="${THREADS:-4}"
mkdir -p "$WORK" && cd "$WORK"

echo "== host: $(hostname)  os: $(uname -sm)  cores: $(getconf _NPROCESSORS_ONLN)"
echo "== python: $(command -v python)"
[[ "$(command -v python)" == *"/envs/kamr/"* ]] || { echo "ERROR: activate the kamr env first"; exit 1; }

echo "== versions"
AMR_VER=$(amrfinder --version)
AMR_DB=$(amrfinder -l 2>&1 | sed -n 's/^Database version: //p' | head -1)
KLEB_VER=$(kleborate --version 2>&1 | tail -1)
KAPT_VER=$(python -c "import importlib.metadata as m; print(m.version('kaptive'))")
python -c "from kaptive.database import load_database" || { echo "ERROR: kaptive.database import failed"; exit 1; }
echo "amrfinder $AMR_VER  db $AMR_DB"
echo "kleborate $KLEB_VER  kaptive $KAPT_VER"
echo "mash $(mash --version)  seqkit $(seqkit version | awk '{print $2}')  blastn $(blastn -version | head -1 | awk '{print $2}')"

if [[ ! -s "ref/$ACC.fna" ]]; then
  echo "== downloading $ACC"
  rm -rf ref ncbi_dataset ref.zip && mkdir -p ref
  datasets download genome accession "$ACC" --include genome --filename ref.zip >/dev/null
  unzip -o -q ref.zip
  mv ncbi_dataset/data/"$ACC"/*.fna "ref/$ACC.fna"
  rm -rf ncbi_dataset ref.zip README.md md5sum.txt
fi
echo "== assembly"
seqkit stats -T "ref/$ACC.fna" | cut -f4,5
echo "fasta sha256: $(python -c "import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest()[:16])" "ref/$ACC.fna")"

echo "== amrfinder (threads=$THREADS)"
start=$(date +%s)
amrfinder -n "ref/$ACC.fna" --organism Klebsiella_pneumoniae --plus --threads "$THREADS" -o amrfinder.tsv -q
echo "amrfinder seconds: $(( $(date +%s) - start ))"

echo "== kleborate"
start=$(date +%s)
rm -rf kleborate && kleborate -a "ref/$ACC.fna" -o kleborate -p kpsc >/dev/null
echo "kleborate seconds: $(( $(date +%s) - start ))"

echo "== key calls"
python - << 'PY'
import glob, hashlib
import pandas as pd

amr = pd.read_csv("amrfinder.tsv", sep="\t")
col = "Element symbol" if "Element symbol" in amr.columns else "Gene symbol"
amr_calls = sorted(amr[col].astype(str))

kl = pd.read_csv(glob.glob("kleborate/*.txt")[0], sep="\t").iloc[0]
want = ["mlst__ST", "Bla_Carb_acquired", "Bla_ESBL_acquired", "Col_mutations", "AGly_acquired",
        "Flq_mutations", "Omp_mutations", "resistance_score__resistance_score",
        "kaptive__K_locus", "kaptive__O_locus"]
kl_calls = {}
for key in want:
    hits = [c for c in kl.index if c.endswith(key)]
    kl_calls[key] = str(kl[hits[0]]) if hits else "MISSING"

expected = {"blaKPC-2", "gyrA_S83I", "parC_S80I", "aac(3)-IId"}
missing = expected - set(amr_calls)
print(f"amrfinder elements ({len(amr_calls)}): {';'.join(amr_calls)}")
for k, v in kl_calls.items():
    print(f"kleborate {k}: {v}")

checks = {
    "amrfinder expected markers": not missing,
    "kleborate ST11": kl_calls["mlst__ST"] == "ST11",
    "kleborate KPC-2": "KPC-2" in kl_calls["Bla_Carb_acquired"],
    "kleborate KL103": kl_calls["kaptive__K_locus"] == "KL103",
}
for name, ok in checks.items():
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
if missing:
    print("missing:", sorted(missing))

blob = "|".join(amr_calls) + "#" + "|".join(f"{k}={v}" for k, v in kl_calls.items())
print("FINGERPRINT", hashlib.sha256(blob.encode()).hexdigest()[:16])
PY
