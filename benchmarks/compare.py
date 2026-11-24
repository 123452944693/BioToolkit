"""Validate BioToolkit and time equivalent tasks on fixed synthetic inputs."""

import argparse
import gzip
import platform
import random
import subprocess
import tempfile
from itertools import product
from math import isclose
from pathlib import Path
from statistics import median
from timeit import repeat

import Bio
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils import gc_fraction

from biotoolkit import (
    __version__,
    find_orfs,
    gc_content,
    global_align,
    quality_stats,
    read_fasta,
    read_fastq,
    reverse_complement,
    translate,
)


def reference_orfs(sequence, min_aa):
    """Check every ATG independently using Biopython for strand and codon translations."""
    result = []
    for strand, oriented in (("+", sequence), ("-", str(Seq(sequence).reverse_complement()))):
        for start in range(len(oriented) - 2):
            if oriented[start:start + 3] != "ATG":
                continue
            complete_end = len(oriented) - (len(oriented) - start) % 3
            translated = str(Seq(oriented[start:complete_end]).translate(table=1))
            stop = translated.find("*")
            if stop < min_aa:
                continue
            end = start + 3 * (stop + 1)
            left, right = (start, end) if strand == "+" else (len(sequence) - end, len(sequence) - start)
            result.append((left, right, strand, start % 3, translated[:stop].translate(str.maketrans("BJZ", "XXX"))))
    return sorted(result)


def validate_sequences(fasta_path, dna):
    alphabet = "ACGTRYSWKMBDHVN"
    ambiguity = str.maketrans("BJZ", "XXX")
    codons = ["".join(bases) for bases in product(alphabet, repeat=3)]
    for codon in codons:
        assert translate(codon) == str(Seq(codon).translate(table=1)).translate(ambiguity), codon
    with fasta_path.open() as handle:
        expected = [(name.split()[0], sequence) for name, sequence in SimpleFastaParser(handle)]
    with fasta_path.open() as handle:
        actual = [(record.id, record.sequence) for record in read_fasta(handle)]
    assert actual == expected
    for _, sequence in expected:
        canonical = "".join(base for base in sequence if base in "ACGT")
        measured_gc = gc_content(sequence)
        if canonical:
            assert isclose(measured_gc, 100 * gc_fraction(canonical))
        else:
            assert measured_gc is None
        assert reverse_complement(sequence) == str(Seq(sequence).reverse_complement())
    fixtures = ["ATGTARTAA", "ATGTANTAA", "TTAYTACAT", "CATGATGGGCTGAC"]
    fixtures.extend(dna(120, alphabet) for _ in range(100))
    for sequence in fixtures:
        for minimum in (1, 5):
            hits = [(h.start, h.end, h.strand, h.frame, h.protein) for h in find_orfs(sequence, minimum)]
            assert sorted(hits) == reference_orfs(sequence, minimum), (sequence, minimum)
    print(f"Validated {len(codons):,} codons, {len(expected):,} FASTA records and {2 * len(fixtures)} ORF checks.")


def validate_quality(path):
    with gzip.open(path, "rt") as handle:
        expected = []
        for record in SeqIO.parse(handle, "fastq"):
            scores = record.letter_annotations["phred_quality"]
            length = len(scores)
            expected.append((record.id, {
                "length": length,
                "mean_phred": sum(scores) / length,
                "q20_percent": 100 * sum(score >= 20 for score in scores) / length,
                "q30_percent": 100 * sum(score >= 30 for score in scores) / length,
            }))
    with gzip.open(path, "rt") as handle:
        actual = [(record.id, quality_stats(record)) for record in read_fastq(handle)]
    assert actual == expected
    print(f"Validated quality statistics for {len(expected):,} reads against decoded Biopython Phred scores.")


def biotoolkit_read_lengths(path):
    with gzip.open(path, "rt") as handle:
        return [(record.id, len(record.sequence)) for record in read_fastq(handle)]


def biopython_read_lengths(path):
    with gzip.open(path, "rt") as handle:
        return [(title.split()[0], len(sequence)) for title, sequence, _ in FastqGeneralIterator(handle)]


def seqkit_read_lengths(executable, path):
    result = subprocess.run(
        [executable, "fx2tab", "-j", "1", "-n", "-i", "-l", str(path)],
        check=True, stdout=subprocess.PIPE, text=True,
    )
    rows = (line.split("\t") for line in result.stdout.splitlines())
    return [(identifier, int(length)) for identifier, length in rows]


def biopython_global(aligner, first, second):
    aligned = aligner.align(first, second)[0]
    return aligned[0], aligned[1], aligned.score


def validate_alignment(result, first, second):
    aligned_first, aligned_second, score = result
    assert aligned_first.replace("-", "") == first and aligned_second.replace("-", "") == second
    assert len(aligned_first) == len(aligned_second)
    assert all((a, b) != ("-", "-") for a, b in zip(aligned_first, aligned_second))
    assert sum(2 if a == b else -2 if "-" in (a, b) else -1
               for a, b in zip(aligned_first, aligned_second)) == score


def benchmark(label, operation):
    operation()
    elapsed = median(repeat(operation, number=1, repeat=5))
    print(f"{label:<38} {elapsed * 1_000:>8.2f} ms")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seqkit", default="seqkit", help="SeqKit executable or full path")
    args = parser.parse_args()
    seqkit_version = subprocess.check_output([args.seqkit, "version"], text=True).strip()
    print(f"BioToolkit {__version__}, Biopython {Bio.__version__}, {seqkit_version}, "
          f"Python {platform.python_version()}, {platform.system()} {platform.machine()}")
    rng = random.Random(2817)

    def dna(length, alphabet="ACGT"):
        return "".join(rng.choices(alphabet, k=length))

    with tempfile.TemporaryDirectory() as temp:
        fasta = Path(temp) / "synthetic.fasta"
        fastq = Path(temp) / "synthetic.fastq.gz"
        with fasta.open("w") as handle:
            for i in range(200):
                seq = dna(173, "ACGTACGTNSWK")
                handle.write(f">seq{i} synthetic\n{seq[:80]}\n{seq[80:]}\n")
            handle.write(">ambiguous_only\nNNSSWW\n")
        with gzip.open(fastq, "wt") as handle:
            for i in range(20_000):
                length = rng.randint(75, 150)
                quality = "".join(chr(q + 33) for q in rng.choices((0, 19, 20, 29, 30, 40), k=length))
                handle.write(f"@read{i} synthetic\n{dna(length)}\n+\n{quality}\n")

        validate_sequences(fasta, dna)
        validate_quality(fastq)
        assert biotoolkit_read_lengths(fastq) == biopython_read_lengths(fastq) == seqkit_read_lengths(args.seqkit, fastq)
        print("Validated all 20,000 FASTQ identifiers and lengths against both reference tools.")

        pairs = [(dna(150), dna(140)) for _ in range(30)]
        aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, gap_score=-2)
        for first, second in pairs:
            result = global_align(first, second)
            reference = biopython_global(aligner, first, second)
            validate_alignment((result.first, result.second, result.score), first, second)
            validate_alignment(reference, first, second)
            assert result.score == reference[2]
        print("Validated 30 global alignment scores and reconstructed alignments.")
        print("\nMedian of five runs (one warmup per tool)")
        benchmark("FASTQ IDs/lengths, BioToolkit", lambda: biotoolkit_read_lengths(fastq))
        benchmark("FASTQ IDs/lengths, Biopython", lambda: biopython_read_lengths(fastq))
        benchmark("FASTQ IDs/lengths, SeqKit", lambda: seqkit_read_lengths(args.seqkit, fastq))
        benchmark("30 global alignments, BioToolkit", lambda: [global_align(a, b) for a, b in pairs])
        benchmark("30 global alignments, Biopython", lambda: [biopython_global(aligner, a, b) for a, b in pairs])


if __name__ == "__main__":
    main()
