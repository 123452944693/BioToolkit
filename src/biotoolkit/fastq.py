"""Strict four-line DNA FASTQ parsing with Phred+33 quality scores."""
from dataclasses import dataclass
from typing import Iterable, Iterator
from .sequence import normalize_dna

@dataclass(frozen=True)
class FastqRecord:
    id: str
    sequence: str
    quality: str

    @property
    def scores(self) -> tuple[int, ...]:
        return tuple(ord(char) - 33 for char in self.quality)

def read_fastq(lines: Iterable[str]) -> Iterator[FastqRecord]:
    """Parse unwrapped four-line records; reject truncated or invalid input."""
    stream = iter(lines)
    for number, header in enumerate(stream, 1):
        header = header.rstrip("\r\n")
        if not header.startswith("@") or not header[1:].strip():
            raise ValueError(f"FASTQ record {number}: invalid header")
        identifier = header[1:].split()[0]
        try:
            seq, plus, quality = (next(stream).rstrip("\r\n") for _ in range(3))
        except (StopIteration, RuntimeError) as exc:
            raise ValueError(f"FASTQ record {number}: truncated record") from exc
        if not plus.startswith("+") or (plus[1:].strip() and plus[1:].split()[0] != identifier):
            raise ValueError(f"FASTQ record {number}: invalid plus line")
        if not seq or any(c.isspace() for c in seq):
            raise ValueError(f"FASTQ record {number}: empty or whitespace-containing sequence")
        seq = normalize_dna(seq)
        if len(seq) != len(quality) or any(not 33 <= ord(c) <= 126 for c in quality):
            raise ValueError(f"FASTQ record {number}: invalid quality length or symbols")
        yield FastqRecord(identifier, seq, quality)

def quality_stats(record: FastqRecord) -> dict:
    """Report length, mean Phred, and percentages at or above Q20 and Q30."""
    scores = record.scores
    n = len(scores)
    if not n or n != len(record.sequence) or any(not 0 <= q <= 93 for q in scores):
        raise ValueError("Invalid FASTQ record")
    return {"length": n, "mean_phred": sum(scores)/n,
            "q20_percent": 100*sum(q >= 20 for q in scores)/n,
            "q30_percent": 100*sum(q >= 30 for q in scores)/n}
