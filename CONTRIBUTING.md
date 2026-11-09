# Contributing and reproducibility reports

This repository is a research artifact accompanying a manuscript, not a general-purpose spin-analysis library. Correctness, reproducibility, portability and documentation contributions are welcome.

## Report an issue

Open a [GitHub issue](https://github.com/dvlahek/detector-induced-spin-identifiability/issues) with:

1. The tested **commit SHA** (frozen manuscript snapshot or current `main`).
2. OS, Python version and dependencies (`python -m pip freeze`).
3. Minimal command, actual output and expected output.
4. Whether the issue affects an analytic proof, numerical evaluation, constructed witness, manuscript text or an auxiliary collider artifact.

Do not post personal credentials or nonpublic event-level data.

## Propose a change

Open a focused pull request against `main`. Identify any affected manuscript equation or table and distinguish scientific changes from presentation-only edits. At minimum run:

```bash
python theory_numerics/verify_publication_tables.py
python theory_numerics/charge_kernel_identity_check.py
```

See the [GitHub Actions reproduction workflow](.github/workflows/reproduce.yml) for the full validation chain. Do **not** silently rewrite immutable cited results or historical archives; explain a correction and create a clearly identified new snapshot when needed.

The auxiliary collider folders are not part of the current Journal of Physics A theoretical claims.
