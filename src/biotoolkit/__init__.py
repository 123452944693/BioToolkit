"""Transparent DNA sequence utilities for teaching and reproducible analysis."""
from .sequence import gc_content, reverse_complement, translate, sequence_stats
from .fasta import FastaRecord, read_fasta

__all__ = ["gc_content", "reverse_complement", "translate", "sequence_stats", "FastaRecord", "read_fasta"]
from .fastq import FastqRecord, read_fastq, quality_stats
from .orf import ORF, find_orfs
from .alignment import Alignment, global_align

__all__ += ["FastqRecord", "read_fastq", "quality_stats", "ORF", "find_orfs", "Alignment", "global_align"]
__version__ = "0.2.1"
