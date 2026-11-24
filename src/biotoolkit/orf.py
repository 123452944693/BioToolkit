"""Find complete ATG-to-stop candidate ORFs in six reading frames."""
from dataclasses import dataclass
from .sequence import _translate_codon, normalize_dna, reverse_complement, translate

@dataclass(frozen=True)
class ORF:
    start: int
    end: int
    strand: str
    frame: int
    protein: str

def find_orfs(sequence: str, min_aa: int = 1) -> list[ORF]:
    """Return complete candidate ORFs, including nested starts.

    Coordinates are 0-based, half-open on the original forward sequence,
    including the stop codon. Protein excludes the terminal stop. Frame is
    the zero-based offset on the scanned strand. Starts must be exact ATG.
    A codon stops an ORF if every possible reading is a standard-code stop.
    """
    if not isinstance(min_aa, int) or min_aa < 1:
        raise ValueError("min_aa must be a positive integer")
    seq = normalize_dna(sequence)
    n = len(seq)
    result = []
    for strand, oriented in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            starts = []
            for i in range(frame, n - 2, 3):
                codon = oriented[i:i+3]
                if codon == "ATG":
                    starts.append(i)
                elif starts and _translate_codon(codon) == "*":
                    for start in starts:
                        if (i - start)//3 >= min_aa:
                            end = i+3
                            a, b = (start, end) if strand == "+" else (n-end, n-start)
                            result.append(ORF(a, b, strand, frame, translate(oriented[start:i])))
                    starts = []
    return sorted(result, key=lambda x: (x.start, x.end, x.strand, x.frame))
