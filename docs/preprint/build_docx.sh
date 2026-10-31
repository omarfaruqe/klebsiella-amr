#!/usr/bin/env bash
# Build Word manuscripts for journals that require .docx (Microbial Genomics does).
#
# reference.docx sets the typeface. Microbial Genomics asks for "a standard typeface such as
# Times New Roman or Arial"; pandoc's default is Calibri, so reference.docx carries a theme
# with Times New Roman for both major and minor fonts.
#
# Pandoc cannot add line numbers to a .docx, and Microbial Genomics requires them for review,
# so enable them once in Word after opening: Layout > Line Numbers > Continuous.
#
# Usage: docs/preprint/build_docx.sh [manuscript.docx] [supplementary.docx]
set -euo pipefail
cd "$(dirname "$0")/../.."

command -v pandoc >/dev/null || { echo "ERROR: pandoc not found" >&2; exit 1; }

REF="docs/preprint/reference.docx"
REF_ARG=()
if [ -f "$REF" ]; then
  REF_ARG=(--reference-doc="$REF")
else
  echo "note: $REF absent; Word output will use pandoc's default typeface" >&2
fi

build () {
  local src="$1" out="$2"
  pandoc "$src" \
    --from=markdown+pipe_tables+autolink_bare_uris \
    --resource-path=.:docs/preprint \
    "${REF_ARG[@]}" \
    -o "$out"
  echo "wrote $out ($(du -h "$out" | cut -f1))"
}

build docs/preprint/manuscript.md "${1:-docs/preprint/manuscript.docx}"
[ -f docs/preprint/supplementary.md ] && \
  build docs/preprint/supplementary.md "${2:-docs/preprint/supplementary.docx}"
