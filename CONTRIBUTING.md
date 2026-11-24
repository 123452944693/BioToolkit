# Contributing

Contributions should improve scientific correctness, robustness, documentation or usability. Open an issue describing the problem, expected behavior, and a minimal synthetic example before substantial changes.

1. Fork the repository and create a feature branch.
2. Install the package in editable mode (`python -m pip install -e .`).
3. Add a test covering the change and run `python -m unittest discover -s tests -v`.
4. Submit a pull request describing the change, test results and biological assumptions. Disclose relevant AI assistance and verify generated code yourself.

Do not include patient data, unpublished sequence data, access tokens or third-party code without permission. Contributions should be compatible with the MIT license. Keep changes focused and avoid trivial or repetitive submissions.

## Initial roadmap

- Compare quality summaries and ORF coordinates against established tools on synthetic fixtures.
- Add optional alternative translation tables with authoritative reference fixtures.
- Document benchmarking against established libraries on synthetic sequences.
- Improve FASTA error messages with record and line context.

These are planned improvements, not implemented features.
