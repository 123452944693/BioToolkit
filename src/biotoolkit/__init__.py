"""Transparent DNA sequence utilities for teaching and reproducible analysis."""
from .sequence import gc_content, reverse_complement, translate, sequence_stats
from .fasta import FastaRecord, read_fasta

__all__ = ["gc_content", "reverse_complement", "translate", "sequence_stats", "FastaRecord", "read_fasta"]
__version__ = "0.1.0"
