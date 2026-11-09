# Manuscript sources

The files `main.tex`, `supplementary.tex`, `references.bib`, `iopjournal.cls` and `figures/` build the article and its Supplementary Material. The manuscript is prepared for journal review and has not been accepted or published.

The numerical data and software are archived at the [fixed scientific commit](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/93c70ed482443730c1abc6c2891732fb1bd5bd4f). All manuscript results use the primary $\alpha=0$ detector model.

## Clean PDF build

Use a LaTeX distribution containing `pdflatex` and `biber`. Start from a clean checkout or remove `.bbl` and `.bcf` files. In this directory, for each document:

```bash
pdflatex -interaction=nonstopmode main.tex
biber main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Repeat the same sequence with `supplementary`. Confirm that the CMS citations and the software reference are numbered, not shown as BibLaTeX citation keys. The paper's numerical source citation points to the fixed commit above.
