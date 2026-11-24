# BioToolkit — DNA sequence utilities

A small, dependency-free Python toolkit for students and researchers learning sequence analysis. BioToolkit emphasizes explicit biological assumptions, readable algorithms and reproducible command-line outputs.

Version 0.1.0 | Python 3.10+ | MIT license

## Features

- Parse multi-record, multiline DNA FASTA files.
- Calculate nucleotide counts and GC percentage, with explicit treatment of ambiguous bases.
- Compute reverse complements with IUPAC ambiguity codes.
- Translate DNA using the standard nuclear genetic code (NCBI table 1), including ambiguous codons.
- Export tab-separated sequence statistics and FASTA-formatted transformations.

## Installation

Clone or download the repository, open a terminal in the project directory, and run:

```sh
python -m pip install .
```

No third-party runtime packages are required. Python 3.10 or newer is needed. The project is not published on PyPI.

## Command-line examples

```sh
biotoolkit stats examples/sample.fasta
biotoolkit revcomp examples/sample.fasta
biotoolkit translate examples/sample.fasta --frame 0 --partial trim
```

Use `-` instead of a file path to read from standard input. Redirect output to a file if needed, for example `biotoolkit stats examples/sample.fasta > stats.tsv`. On invalid input, the command exits with status 1 and prints an error to stderr. Because parsing is streaming, a malformed later record may leave partial output; discard that output on failure.

## Python API

```python
from biotoolkit import gc_content, reverse_complement, translate, read_fasta

print(gc_content("ACGTNN"))       # 50.0
print(reverse_complement("ATGC")) # GCAT
print(translate("ATGGCNTAA"))     # MA*

with open("examples/sample.fasta", encoding="utf-8") as handle:
    for record in read_fasta(handle):
        print(record.id, record.sequence)
```

## Scientific conventions

| Operation | Convention |
| --- | --- |
| Alphabet | DNA IUPAC symbols A C G T R Y S W K M B D H V N; no RNA U or gaps |
| GC percentage | 100 × (G + C) / (A + C + G + T), excluding ambiguous positions; undefined as `None` / `NA` when the denominator is zero |
| Reverse complement | IUPAC complements, output in 5′ to 3′ orientation |
| Translation | Standard genetic code, no start codon inference; stops shown as `*`; ambiguous codons resolve only when all possibilities encode the same symbol, otherwise `X` |
| Reading frame | Zero-based offset 0, 1 or 2; reverse-strand translation requires reverse complement first |
| Partial codons | Error by default; optional trimming of incomplete trailing codons |
| FASTA | Empty records and missing identifiers rejected; duplicate identifiers retained |

The parser processes one record at a time but stores that record's full sequence in memory. This tool is educational and intended for small-scale analyses; it is not validated for clinical diagnosis, variant interpretation or production-scale genome processing. It does not infer open reading frames, handle splicing or support alternative genetic codes.

## Tests

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
```

Tests cover all 64 standard codons, ambiguity handling, invalid inputs, FASTA parsing and command-line behavior. A GitHub Actions workflow runs these tests on Python 3.10–3.13 when configured in a repository with Actions enabled.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Meaningful bug reports, reproducible test cases and scientifically justified improvements are welcome.

## References

- NCBI, [The Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — standard translation table 1.
- IUBMB, [Nomenclature for incompletely specified bases in nucleic acid sequences](https://iubmb.qmul.ac.uk/misc/naseq.html) — ambiguity codes.
- NCBI, [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) — background on FASTA sequence representation.

## License

This independent toolkit is released under the MIT license. See [LICENSE](LICENSE). It contains no Helix or SmartLib source code and uses only synthetic examples.
