"""Command-line interface. Output goes to stdout, errors to stderr."""
import argparse
import csv
import gzip
import sys
import zlib
from contextlib import nullcontext
from . import __version__
from .fasta import read_fasta
from .fastq import read_fastq, quality_stats
from .orf import find_orfs
from .sequence import reverse_complement, sequence_stats, translate


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="biotoolkit", description="Analyze DNA FASTA and FASTQ files")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    commands = {
        "stats": "Summarize base counts and GC percentage in FASTA records",
        "revcomp": "Reverse complement FASTA sequences",
        "translate": "Translate FASTA sequences using the standard genetic code",
        "fastq-stats": "Summarize per-read Phred+33 quality in FASTQ records",
        "orfs": "Find candidate open reading frames in FASTA sequences",
    }
    for command, description in commands.items():
        p = sub.add_parser(command, help=description, description=description)
        p.add_argument("input", help="Input path (.gz supported), or '-' for standard input")
        if command == "orfs":
            p.add_argument("--min-aa", type=int, default=1, help="Minimum protein length (default: 1)")
        if command == "translate":
            p.add_argument("--frame", type=int, choices=(0, 1, 2), default=0, help="Zero-based reading frame")
            p.add_argument("--partial", choices=("error", "trim"), default="error", help="Handle a trailing partial codon")
    args = parser.parse_args(argv)
    try:
        if args.input == "-":
            context = nullcontext(sys.stdin)
        elif args.input.endswith(".gz"):
            context = gzip.open(args.input, "rt", encoding="utf-8")
        else:
            context = open(args.input, encoding="utf-8")
        with context as stream:
            if args.command == "fastq-stats":
                out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
                out.writerow(["id", "length", "mean_phred", "q20_percent", "q30_percent"])
                for record in read_fastq(stream):
                    q = quality_stats(record)
                    out.writerow([record.id, q["length"], *[f"{q[key]:.4f}" for key in ("mean_phred", "q20_percent", "q30_percent")]])
                return 0
            records = read_fasta(stream)
            if args.command == "orfs":
                if args.min_aa < 1:
                    raise ValueError("min-aa must be positive")
                out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
                out.writerow(["id", "start", "end", "strand", "frame", "protein"])
                for record in records:
                    for hit in find_orfs(record.sequence, args.min_aa):
                        out.writerow([record.id, hit.start, hit.end, hit.strand, hit.frame, hit.protein])
                return 0
            if args.command == "stats":
                out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
                out.writerow(["id", "length", "A", "C", "G", "T", "ambiguous", "gc_percent"])
                for record in records:
                    stats = sequence_stats(record.sequence)
                    gc = stats.pop("gc_percent")
                    out.writerow([record.id, *stats.values(), "NA" if gc is None else f"{gc:.4f}"])
            else:
                for record in records:
                    seq = reverse_complement(record.sequence) if args.command == "revcomp" else translate(record.sequence, args.frame, partial=args.partial)
                    print(f">{record.id} {record.description}".rstrip())
                    print("\n".join(seq[i:i+80] for i in range(0, len(seq), 80)))
        return 0
    except (ValueError, OSError, EOFError, zlib.error) as exc:
        print(f"biotoolkit: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
