# Comparing BioToolkit with Biopython and SeqKit

Run from the repository root in a virtual environment with [SeqKit](https://bioinf.shenwei.me/seqkit/download/) installed.

```sh
python -m pip install . biopython==1.88
python benchmarks/compare.py
```

Use `--seqkit /path/to/seqkit` if the executable is outside your PATH. The script generates its own synthetic fixtures with a fixed random seed and exits with an error on a comparison failure. No sequence downloads are needed.

## Correctness checks

| Method | Reference and result |
| --- | --- |
| Standard translation | All 3,375 three-base combinations over the 15-symbol DNA alphabet agreed with [Biopython table 1 translation](https://biopython.org/docs/latest/api/Bio.Seq.html). Biopython's B, J and Z ambiguity symbols were mapped to BioToolkit's X. |
| FASTA, GC and reverse complement | 201 records matched Biopython. GC used only A/C/G/T in the denominator. The all-ambiguity fixture correctly returned an undefined GC percentage. |
| Candidate ORFs | 208 checks matched a separate reference implementation using Biopython for strand orientation and translation. The fixtures covered both strands, nested starts, exact stops, guaranteed ambiguous stops and possible stops. |
| FASTQ | All 20,000 read IDs and lengths matched Biopython and [SeqKit](https://bioinf.shenwei.me/seqkit/usage/). Mean Phred and Q20/Q30 matched calculations from Phred scores decoded by Biopython. SeqKit's average-quality definition uses mean error probabilities. |
| Global alignment | All 30 optimal scores matched [Biopython's PairwiseAligner](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html) with match 2, mismatch −1 and linear gap −2. Both outputs reconstructed the original sequences and rescored to the reported result. Different tracebacks can tie for the same optimal score. |

## Runtime on these inputs

BioToolkit 0.2.1, Biopython 1.88, SeqKit 2.14.0, CPython 3.12.14, Linux x86-64 on an AMD EPYC 7763. Times are the median of five runs after one warmup.

| Task | BioToolkit | Biopython | SeqKit |
| --- | ---: | ---: | ---: |
| Extract IDs and lengths from 20,000 gzip-compressed 75–150 base reads | 423.46 ms | 64.14 ms | 73.27 ms |
| Score and reconstruct 30 global alignments of 150 × 140 bases | 207.28 ms | 12.07 ms | — |

The FASTQ comparison includes each parser's usual validation and the SeqKit subprocess startup, gzip input and TSV decoding. SeqKit used one thread. The alignment comparison includes traceback in both implementations. These are local timings on small synthetic datasets. Use the reference tools for throughput-oriented work.
