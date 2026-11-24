"""DNA utilities. Ambiguity codes follow the IUPAC nucleotide alphabet."""
from collections import Counter
from itertools import product

IUPAC = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"AG", "Y":"CT", "S":"CG", "W":"AT", "K":"GT", "M":"AC", "B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG", "N":"ACGT"}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
# NCBI translation table 1, codon order T,C,A,G in each position.
_CODONS = dict(zip(("".join(c) for c in product("TCAG", repeat=3)),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def normalize_dna(sequence: str) -> str:
    """Uppercase DNA and remove whitespace; reject non-IUPAC symbols."""
    seq = "".join(sequence.split()).upper()
    invalid = sorted(set(seq) - IUPAC.keys())
    if invalid:
        raise ValueError("Invalid DNA symbols: " + ", ".join(invalid))
    return seq


def gc_content(sequence: str) -> float | None:
    """GC percentage among unambiguous A/C/G/T bases; None if none exist."""
    seq = normalize_dna(sequence)
    counts = Counter(seq)
    denominator = sum(counts[b] for b in "ACGT")
    return 100 * (counts["G"] + counts["C"]) / denominator if denominator else None


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement, preserving IUPAC ambiguity codes."""
    return normalize_dna(sequence).translate(_COMPLEMENT)[::-1]


def translate(sequence: str, frame: int = 0, *, partial: str = "error") -> str:
    """Translate DNA using standard code (table 1), no start-codon inference.

    frame: 0, 1 or 2 from the supplied 5' end. Stop codons are '*'.
    Ambiguous codons resolve to a residue only if all possibilities agree,
    otherwise 'X'. Incomplete trailing codons raise unless partial='trim'.
    """
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1 or 2")
    if partial not in ("error", "trim"):
        raise ValueError("partial must be 'error' or 'trim'")
    seq = normalize_dna(sequence)[frame:]
    remainder = len(seq) % 3
    if remainder and partial == "error":
        raise ValueError("Incomplete trailing codon; use partial='trim' to discard it")
    result = []
    for i in range(0, len(seq) - remainder, 3):
        residues = {_CODONS["".join(c)] for c in product(*(IUPAC[b] for b in seq[i:i+3]))}
        result.append(next(iter(residues)) if len(residues) == 1 else "X")
    return "".join(result)


def sequence_stats(sequence: str) -> dict:
    """Length, base counts, ambiguity count, and GC percentage."""
    seq = normalize_dna(sequence)
    counts = Counter(seq)
    bases = {b: counts[b] for b in "ACGT"}
    return {"length": len(seq), **bases, "ambiguous": len(seq) - sum(bases.values()), "gc_percent": gc_content(seq)}
