# Hands-on sequence analysis tutorial

This tutorial uses synthetic sequences only. No sequencing instrument or biological sample is required. The goal is to learn the conventions behind each output, not to make clinical or gene-annotation decisions.

## 1. Install and check the environment

From the repository root, with Python 3.10 or later:

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
```

## 2. Calculate sequence statistics

```sh
biotoolkit stats examples/sample.fasta
```

For the first record, ATGGCNTAA, length is 9; counts are A=3, C=1, G=2, T=2 and ambiguous=1. GC percentage is 100 × 3/8 = 37.5%. Ambiguous bases are excluded from the denominator. Compare this value with the output. For the second record, calculate the count of ambiguous symbols by hand.

## 3. Translate a coding sequence

```sh
biotoolkit translate examples/sample.fasta --partial trim
```

The first record produces MA*: ATG encodes methionine, GCN encodes alanine for all four possible third bases, and TAA is a stop. BioToolkit does not infer a start or apply alternative translation tables.

## 4. Inspect read quality

```sh
biotoolkit fastq-stats examples/sample.fastq
```

The quality string !5?I corresponds to Phred scores 0, 20, 30 and 40. The arithmetic mean is 22.5, Q20 is 75%, and Q30 is 50%. Note that an arithmetic mean of Phred scores is not the same as an average error probability.

## 5. Explore candidate ORFs and alignment

```python
from biotoolkit import find_orfs, global_align
print(find_orfs("ATGGCCTAA"))
print(global_align("ACGT", "AGT"))
```

The forward ORF has coordinates [0, 9), frame 0, and protein MA. The alignment has optimal score 4 under the default scoring scheme (three matches and one gap). Candidate ORFs do not establish that a gene is expressed. Alignment ties can produce more than one optimal solution.

## 6. Give useful feedback

If you are testing the project, please note your Python version, command, expected behavior, actual output and a minimal synthetic example. Report reproducible issues using the repository's issue tracker. Do not upload private or human sequence data.

Useful feedback includes confusing documentation, scientific conventions that need clarification, incorrect output, or missing teaching use cases. Feedback and contributions are entirely voluntary.
