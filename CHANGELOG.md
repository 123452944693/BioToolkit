# Changelog

## 0.2.1

- Stop candidate ORFs at ambiguous codons that always translate to a stop.
- Reject non-ASCII characters before normalizing DNA, include input context in parser errors, and handle corrupt gzip input through the CLI error path.
- Add command descriptions and a `--version` flag.
- Build and check installable distributions in CI, and add Windows and Python 3.14 test coverage.
- Record reproducible synthetic-data comparisons with Biopython and SeqKit.

## 0.2.0

- Strict four-line FASTQ parsing, Phred+33 summary statistics and gzip file input.
- Complete six-frame ATG-to-stop candidate ORF discovery with explicit coordinates.
- Needleman–Wunsch global alignment API for small DNA sequences.
- Additional synthetic examples and regression tests.

## 0.1.0

- Initial FASTA, GC content, reverse complement and translation utilities.
