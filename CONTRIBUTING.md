# Reproducibility reports and contributions

The repository is a reproducible research artifact associated with a mathematical-physics manuscript. It is not a general-purpose spin-determination package.

If a calculation fails, open a GitHub issue with the following details.

1. Exact commit SHA and the command you ran.
2. Operating system, Python version and `python -m pip freeze`.
3. Expected and actual output, including the first relevant traceback.
4. The affected manuscript equation, figure or evidence table.

A minimal initial check is:

```bash
python -m pip install -r requirements.txt
python reproduce.py --mode smoke
```

Use `python reproduce.py --mode full` for the full numerical regression suite. Pull requests should distinguish a change to the scientific model from a correction to implementation or presentation. Never overwrite a frozen cited result without documenting a new reproducibility snapshot.

Do not include private data or credentials in issues. See [the workflow](.github/workflows/reproduce.yml) and the [reproduction guide](docs/REPRODUCIBILITY.md).
