#!/usr/bin/env bash
# Build a submission-ready PDF from the manuscript Markdown.
#
# Requires pandoc and xelatex. Figure paths in the manuscript are relative to docs/preprint, so
# the resource path is set to the repository root.
#
# The Markdown is preprocessed so the PDF does not duplicate what the metadata block provides:
#   - the leading H1 title and byline are removed (metadata.yaml supplies them);
#   - image alt text is stripped, because pandoc turns it into a second caption above the
#     manuscript's own "**Figure N.**" caption. Alt text is kept in the Markdown for GitHub.
#
# Usage: docs/preprint/build_pdf.sh [manuscript.pdf] [supplementary.pdf]
set -euo pipefail
cd "$(dirname "$0")/../.."

OUT="${1:-docs/preprint/manuscript.pdf}"
for tool in pandoc xelatex; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found" >&2; exit 1; }
done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
python3 - "$TMP/body.md" <<'PY'
import re
import sys

src = open("docs/preprint/manuscript.md").read()
lines = src.split("\n")
# drop the leading H1 and the byline line that follows it; metadata.yaml carries both
out, dropped_title, dropped_author = [], False, False
for line in lines:
    if not dropped_title and line.startswith("# "):
        dropped_title = True
        continue
    if dropped_title and not dropped_author and line.strip() == "Md Omar Faruqe":
        dropped_author = True
        continue
    out.append(line)
body = "\n".join(out)
# strip alt text so pandoc does not emit its own figure caption
body = re.sub(r"!\[[^\]]*\]\(", "![](", body)
open(sys.argv[1], "w").write(body)
PY

pandoc "$TMP/body.md" \
  --from=markdown+pipe_tables \
  --pdf-engine=xelatex \
  --resource-path=.:docs/preprint \
  --metadata-file=docs/preprint/metadata.yaml \
  -V geometry:"a4paper,margin=2.4cm" \
  -V mainfont="DejaVu Serif" \
  -V sansfont="DejaVu Sans" \
  -V monofont="DejaVu Sans Mono" \
  -V fontsize=10pt \
  -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
  -o "$OUT"

echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"

# Supplementary material builds as its own PDF, so it can be uploaded separately.
SUPP_SRC="docs/preprint/supplementary.md"
SUPP_OUT="${2:-docs/preprint/supplementary.pdf}"
if [ -f "$SUPP_SRC" ]; then
  pandoc "$SUPP_SRC" \
    --from=markdown+pipe_tables \
    --pdf-engine=xelatex \
    --resource-path=.:docs/preprint \
    -V geometry:"a4paper,margin=2.4cm" \
    -V mainfont="DejaVu Serif" \
    -V sansfont="DejaVu Sans" \
    -V monofont="DejaVu Sans Mono" \
    -V fontsize=10pt \
    -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
    -o "$SUPP_OUT"
  echo "wrote $SUPP_OUT ($(du -h "$SUPP_OUT" | cut -f1))"
fi
