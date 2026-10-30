#!/usr/bin/env bash
# Build a Word manuscript for journals that require .docx (Microbial Genomics does).
#
# Pandoc cannot add line numbers to a .docx, and Microbial Genomics requires them, so enable
# them once in Word after opening: Layout > Line Numbers > Continuous.
#
# Figure images are embedded. Separate high-resolution figure files are required at revision and
# live in reports/figures/.
#
# Usage: docs/preprint/build_docx.sh [output.docx] [supplementary.docx]
set -euo pipefail
cd "$(dirname "$0")/../.."

OUT="${1:-docs/preprint/manuscript.docx}"
command -v pandoc >/dev/null || { echo "ERROR: pandoc not found" >&2; exit 1; }

pandoc docs/preprint/manuscript.md \
  --from=markdown+pipe_tables \
  --resource-path=.:docs/preprint \
  --reference-doc=docs/preprint/reference.docx 2>/dev/null \
  -o "$OUT" \
  || pandoc docs/preprint/manuscript.md \
       --from=markdown+pipe_tables \
       --resource-path=.:docs/preprint \
       -o "$OUT"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"

SUPP_OUT="${2:-docs/preprint/supplementary.docx}"
if [ -f docs/preprint/supplementary.md ]; then
  pandoc docs/preprint/supplementary.md \
    --from=markdown+pipe_tables \
    --resource-path=.:docs/preprint \
    -o "$SUPP_OUT"
  echo "wrote $SUPP_OUT ($(du -h "$SUPP_OUT" | cut -f1))"
fi
