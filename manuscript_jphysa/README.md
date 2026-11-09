# Manuscript sources

The files `main.tex`, `supplementary.tex`, `references.bib`, `iopjournal.cls` and `figures/` build the article and its Supplementary Material. The manuscript is prepared for journal review and has not been accepted or published.

The numerical data and software are archived at the [fixed scientific commit](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/93c70ed482443730c1abc6c2891732fb1bd5bd4f). All manuscript results use the primary $\alpha=0$ detector model.

## Clean PDF build

With a LaTeX distribution that provides `pdflatex`, `biber` and `pdftotext`, run from the repository root:

```bash
bash manuscript_jphysa/build_clean.sh
```

The script compiles both documents in a fresh temporary directory using the complete `pdflatex → biber → pdflatex → pdflatex` sequence. It writes `main_biber.pdf` and `supplementary_biber.pdf` to `submission_pdfs/` and stops if citation keys remain unresolved or the manuscript points to the outdated filtered-charge code snapshot.

To select another output directory, pass its path as the first argument. The sources and frozen numerical snapshot are never modified by this build.
