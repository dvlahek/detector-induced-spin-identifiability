# J. Phys. A manuscript source — unfiltered charge (v1.9.0)

This is the submitted-candidate LaTeX source for the main article and Supplementary Material. The manuscript is **not yet published or peer-reviewed**.

**Frozen numerical code and data:** [`f711505b03b1a01ad4d4be8a50ce4d46dbda2080`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/f711505b03b1a01ad4d4be8a50ce4d46dbda2080). This immutable commit contains the primary unfiltered (`alpha=0`) results; it must not be confused with the earlier filtered-charge v1.8.0 snapshot.

Build from this directory:

```bash
pdflatex -interaction=nonstopmode main.tex
biber main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary.tex
biber supplementary
pdflatex -interaction=nonstopmode supplementary.tex
pdflatex -interaction=nonstopmode supplementary.tex
```

The IOP class is included for standalone compilation. The typesetting placeholders are suppressed for review PDF builds. The cover letter is separate from the public research-code repository.
