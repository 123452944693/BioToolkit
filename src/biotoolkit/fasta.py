"""Streaming parser for DNA FASTA files."""
from dataclasses import dataclass
from typing import Iterable, Iterator
from .sequence import normalize_dna

@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def read_fasta(lines: Iterable[str]) -> Iterator[FastaRecord]:
    """Yield records, reject missing headers, empty IDs and empty sequences.

    Duplicate IDs are retained; this parser does not impose uniqueness.
    """
    header = None
    chunks: list[str] = []

    def record() -> FastaRecord:
        assert header is not None
        if not chunks:
            raise ValueError(f"Empty sequence for {header!r}")
        parts = header.split(maxsplit=1)
        return FastaRecord(parts[0], parts[1] if len(parts) > 1 else "", "".join(chunks))

    for n, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield record()
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Missing FASTA identifier at line {n}")
            chunks = []
        else:
            if header is None:
                raise ValueError(f"Sequence before FASTA header at line {n}")
            try:
                chunks.append(normalize_dna(line))
            except ValueError as exc:
                raise ValueError(f"FASTA line {n} for {header!r}: {exc}") from exc
    if header is None:
        raise ValueError("No FASTA records found")
    yield record()
