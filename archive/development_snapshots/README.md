# Development snapshot archive

This directory preserves superseded files from the development releases
`v1.0.0`, `v1.1.0`, `v1.4.0`, and `v1.5.0`.

- `workflows/` contains retired release-creation workflows and is outside
  `.github/workflows/`, so the legacy automation cannot run.
- `release_notes/` contains the corresponding historical notes.
- `theory_numerics/` contains superseded scripts and versioned result trees.
- `provenance/` retains the original internal workflow/artifact labels; the
  active provenance file keeps the scientific identifiers but removes obsolete
  submission-branding and transient workflow-run labels.

The historical Git tags and GitHub releases remain unchanged. Publication 1.0
uses only the unversioned active scripts and `results/publication_v1_0/` in the
repository root.
