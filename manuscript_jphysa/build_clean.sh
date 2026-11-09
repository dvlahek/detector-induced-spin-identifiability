#!/usr/bin/env bash
# Clean, reproducible PDF compilation for the SPIN / Journal of Physics A manuscript.
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${1:-${SRC_DIR}/../submission_pdfs}"
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

for tool in pdflatex biber pdftotext; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required executable not available: $tool" >&2
    exit 2
  fi
done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$SRC_DIR"/*.tex "$SRC_DIR"/*.bib "$SRC_DIR"/*.cls "$TMP"/
cp -R "$SRC_DIR/figures" "$TMP/"
cd "$TMP"

for doc in main supplementary; do
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" > "$OUT_DIR/${doc}_pdflatex_first.log" 2>&1
  biber "$doc" > "$OUT_DIR/${doc}_biber.log" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" > "$OUT_DIR/${doc}_pdflatex_second.log" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" > "$OUT_DIR/${doc}_pdflatex_final.log" 2>&1
  pdftotext -layout "$doc.pdf" "$doc.txt"

  if grep -Eq 'CMS2017ChargedTGC|CMS2026EFT|Lasserre2001|Parrilo2003|5c48258f69758818e3ed524cb4d7c1e9d5697811' "$doc.txt"; then
    echo "Incorrect or unresolved bibliography in $doc.pdf" >&2
    exit 3
  fi
  if grep -Eq '^WARN -|^ERROR -' "$OUT_DIR/${doc}_biber.log"; then
    echo "Biber reported a warning/error for $doc" >&2
    exit 4
  fi
  cp "$doc.pdf" "$OUT_DIR/${doc}_biber.pdf"
done

if ! grep -Fq '93c70ed482443730c1abc6c2891732fb1bd5bd4f' "$TMP/main.txt"; then
  echo 'The PDF does not cite the frozen alpha=0 scientific commit' >&2
  exit 5
fi
if ! grep -Fq '10.1016/j.physletb.2017.06.009' "$TMP/main.txt"; then
  echo 'The CMS 2017 DOI is missing from the PDF' >&2
  exit 6
fi
if ! grep -Fq 'Lasserre' "$TMP/supplementary.txt" || ! grep -Fq 'Parrilo' "$TMP/supplementary.txt"; then
  echo 'Expected SOS references are missing from the supplementary PDF' >&2
  exit 7
fi

echo "Clean builds and bibliography checks passed. PDF files are in $OUT_DIR"
