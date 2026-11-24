# BioToolkit — DNA sequence utilities

A small, dependency-free Python toolkit for students and researchers learning sequence analysis. BioToolkit emphasizes explicit biological assumptions, readable algorithms and reproducible command-line outputs.

Python 3.10+ | MIT license

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
biotoolkit --version
```

If your terminal cannot find the command, use `python -m biotoolkit` in its place.
No third-party runtime packages are required. Python 3.10 or newer is needed.
The [GitHub Actions builds](https://github.com/123452944693/BioToolkit/actions/workflows/tests.yml) provide a wheel and source distribution under `biotoolkit-distributions` for signed-in GitHub users.

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
| Alphabet | ASCII DNA IUPAC symbols A C G T R Y S W K M B D H V N, case-insensitive; no RNA U or gaps |
| GC percentage | 100 × (G + C) / (A + C + G + T), excluding ambiguous positions; undefined as `None` / `NA` when the denominator is zero |
| Reverse complement | IUPAC complements, output in 5′ to 3′ orientation |
| Translation | Standard genetic code, no start codon inference; stops shown as `*`; ambiguous codons resolve only when all possibilities encode the same symbol, otherwise `X` |
| Reading frame | Zero-based offset 0, 1 or 2; reverse-strand translation requires reverse complement first |
| Partial codons | Error by default; optional trimming of incomplete trailing codons |
| FASTA | Empty records and missing identifiers rejected; duplicate identifiers retained |

The parser stores one full record at a time. Use trusted, small inputs, especially when decompressing files or scanning sequences with many nested ORFs. When opening TSV output in a spreadsheet, import the identifier column as text. BioToolkit is designed for education and small-scale analyses. Clinical work, production-scale processing, splicing and alternative genetic codes are outside its scope.

### New analysis APIs

```python
from biotoolkit import find_orfs, global_align
print(find_orfs("ATGGCCTAA"))
print(global_align("ACGT", "AGT"))
```

FASTQ accepts unwrapped four-line records and Phred+33 scores (ASCII 33–126). Mean quality is the arithmetic mean of Phred scores. Q20/Q30 percentages count bases meeting or exceeding those thresholds. FASTA and FASTQ CLI commands support `.gz` files.

Candidate ORFs begin at exact ATG and end at a codon that can only encode a standard-code stop, such as TAA or TAR. A codon that could encode either a residue or a stop appears as `X` within the candidate protein. Results include nested starts. Coordinates are 0-based, half-open and relative to the forward reference, with the stop codon included. The returned protein excludes the terminal stop. Many nested starts can generate a large amount of output.

Alignment uses Needleman–Wunsch with linear gaps, literal IUPAC comparison and a one-million-cell limit.

## Tests

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
```

Tests cover all 64 standard codons, ambiguity handling, invalid inputs, FASTA and FASTQ parsing, command-line behavior, gzip input, ORF coordinates and global alignment. GitHub Actions runs the tests on Linux with Python 3.10–3.14 and Windows with Python 3.14.

## Tutorial

Try the [hands-on sequence analysis tutorial](docs/tutorial.md) with synthetic FASTA and FASTQ files, worked calculations, and expected results.

See the [benchmark methods and results](docs/benchmarks.md) for comparisons with Biopython and SeqKit.

## Development transparency

The initial implementation, documentation and automated tests were produced with substantial ChatGPT assistance and published at the repository owner's direction. GitHub CI verifies the included tests; it does not establish independent scientific validation, adoption or clinical suitability. The owner and any future contributors should independently review and understand changes before relying on them.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Meaningful bug reports, reproducible test cases and scientifically justified improvements are welcome.

## References

- NCBI, [The Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) — standard translation table 1.
- IUBMB, [Nomenclature for incompletely specified bases in nucleic acid sequences](https://iubmb.qmul.ac.uk/misc/naseq.html) — ambiguity codes.
- NCBI, [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) — background on FASTA sequence representation.

## License

This toolkit is released under the MIT license. See [LICENSE](LICENSE). All examples are synthetic.
