# BioToolkit — DNA sequence utilities

A small, dependency-free Python toolkit for students and researchers learning sequence analysis. BioToolkit emphasizes explicit biological assumptions, readable algorithms and reproducible command-line outputs.

Version 0.2.0 | Python 3.10+ | MIT license

## Features

- Parse multi-record, multiline DNA FASTA files.
- Calculate nucleotide counts and GC percentage, with explicit treatment of ambiguous bases.
- Compute reverse complements with IUPAC ambiguity codes.
- Translate DNA using the standard nuclear genetic code (NCBI table 1), including ambiguous codons.
- Export tab-separated sequence statistics and FASTA-formatted transformations.

- Analyze four-line FASTQ reads (Phred+33), including gzip-compressed input.
- Find complete candidate ORFs across both strands.
- Perform small global DNA alignments via a Python API.

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
biotoolkit fastq-stats examples/sample.fastq
biotoolkit orfs examples/sample.fasta --min-aa 1
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

The parser processes one record at a time but stores that record's full sequence in memory. This tool is educational and intended for small-scale analyses; it is not validated for clinical diagnosis, variant interpretation or production-scale genome processing. ORF candidates are not validated gene predictions. It does not handle splicing or support alternative genetic codes.

### New analysis APIs

```python
from biotoolkit import find_orfs, global_align
print(find_orfs("ATGGCCTAA"))
print(global_align("ACGT", "AGT"))
```

FASTQ accepts unwrapped four-line records and Phred+33 scores (ASCII 33–126). Mean quality is the arithmetic mean of Phred scores, not mean error probability. Q20/Q30 percentages count bases meeting or exceeding those thresholds. FASTA and FASTQ CLI commands support `.gz` files. ORFs require ATG and an exact standard-code stop, include nested starts, report 0-based half-open forward-reference coordinates including the stop, and exclude the stop from the returned protein. Ambiguous stops are not used to terminate candidates. Alignment uses Needleman–Wunsch with linear gaps, literal IUPAC comparison and a one-million-cell limit; it is not a production-scale aligner.

## Tests

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
```

Tests cover all 64 standard codons, ambiguity handling, invalid inputs, FASTA parsing command-line behavior, FASTQ quality, gzip input, ORF coordinates and global alignment. A GitHub Actions workflow runs these tests on Python 3.10–3.13 when configured in a repository with Actions enabled.

## Tutorial

Try the [hands-on sequence analysis tutorial](docs/tutorial.md) with synthetic FASTA and FASTQ files, worked calculations, and expected results.

## Development transparency

The initial implementation, documentation and automated tests were produced with substantial ChatGPT assistance and published at the repository owner's direction. GitHub CI verifies the included tests; it does not establish independent scientific validation, adoption or clinical suitability. The owner and any future contributors should independently review and understand changes before relying on them.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Meaningful bug reports, reproducible test cases and scientifically justified improvements are welcome.

## References

- NCBI, [The Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — standard translation table 1.
- IUBMB, [Nomenclature for incompletely specified bases in nucleic acid sequences](https://iubmb.qmul.ac.uk/misc/naseq.html) — ambiguity codes.
- NCBI, [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) — background on FASTA sequence representation.

## License

This independent toolkit is released under the MIT license. See [LICENSE](LICENSE). It contains no Helix or SmartLib source code and uses only synthetic examples.
