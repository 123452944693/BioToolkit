import unittest
import gzip
import tempfile
import zlib
from pathlib import Path
from io import StringIO
from contextlib import redirect_stderr, redirect_stdout
from unittest.mock import patch
from biotoolkit import read_fastq, quality_stats, find_orfs, global_align
from biotoolkit.cli import main

class FastqTests(unittest.TestCase):
    def test_quality(self):
        rec = next(read_fastq(StringIO("@r desc\nACGT\n+r\n!5?I\n")))
        self.assertEqual(rec.scores, (0,20,30,40))
        self.assertEqual(quality_stats(rec), {"length":4,"mean_phred":22.5,"q20_percent":75,"q30_percent":50})
    def test_invalid(self):
        for text in ("@r\nA\n+\n", "@r\nA\n+\n!!\n", "@r\nA\n+wrong\n!\n", "@r\nA\n+\n \n", "@r\nA C\n+\n!!!\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                list(read_fastq(StringIO(text)))
    def test_invalid_sequence_context(self):
        data = "@first\nA\n+\n!\n@second\nſ\n+\n!\n"
        with self.assertRaisesRegex(ValueError, "FASTQ record 2 .*Invalid DNA"):
            list(read_fastq(StringIO(data)))
    def test_gzip_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"example.fastq.gz"
            with gzip.open(path,"wt",encoding="utf-8") as out:
                out.write("@r\nAT\n+\nII\n")
            output = StringIO()
            with redirect_stdout(output):
                self.assertEqual(main(["fastq-stats", str(path)]),0)
            self.assertIn("r\t2\t40.0000\t100.0000\t100.0000",output.getvalue())
    def test_invalid_gzip_reports_cli_error(self):
        with patch("biotoolkit.cli.gzip.open", side_effect=zlib.error("invalid compressed data")):
            error = StringIO()
            with redirect_stderr(error):
                self.assertEqual(main(["fastq-stats", "reads.fastq.gz"]), 1)
            self.assertEqual(error.getvalue(), "biotoolkit: invalid compressed data\n")

class ORFTests(unittest.TestCase):
    def test_forward(self):
        hits = [x for x in find_orfs("CATGGCCTAAC") if x.strand == "+"]
        self.assertEqual([(x.start,x.end,x.frame,x.protein) for x in hits],[(1,10,1,"MA")])
    def test_reverse(self):
        hits = [x for x in find_orfs("GTTAGGCCATG") if x.strand == "-"]
        self.assertEqual([(x.start,x.end,x.frame,x.protein) for x in hits],[(1,10,1,"MA")])
    def test_nested_and_minimum(self):
        hits = [x for x in find_orfs("ATGATGTAA") if x.strand == "+"]
        self.assertEqual([x.protein for x in hits],["MM","M"])
        self.assertEqual(len([x for x in find_orfs("ATGATGTAA",2) if x.strand == "+"]),1)
    def test_unambiguous_stop_from_ambiguous_bases(self):
        for seq, strand in (("ATGTARTAA", "+"), ("TTAYTACAT", "-")):
            with self.subTest(strand=strand):
                hits = [hit for hit in find_orfs(seq) if hit.strand == strand]
                start, end = (0, 6) if strand == "+" else (3, 9)
                self.assertEqual([(h.start, h.end, h.protein) for h in hits], [(start, end, "M")])
    def test_possible_stop_is_translated_as_unknown(self):
        hits = [hit for hit in find_orfs("ATGTANTAA") if hit.strand == "+"]
        self.assertEqual([(h.start, h.end, h.protein) for h in hits], [(0, 9, "MX")])
    def test_incomplete(self):
        self.assertEqual(find_orfs("ATGCCC"),[])
        with self.assertRaises(ValueError): find_orfs("ATG",0)

class AlignmentTests(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(global_align("ACGT","ACGT").score,8)
    def test_indel(self):
        result = global_align("ACGT","AGT")
        self.assertEqual(result.score,4)
        self.assertEqual(result.first.replace("-",""),"ACGT")
        self.assertEqual(result.second.replace("-",""),"AGT")
        self.assertEqual(len(result.first),len(result.second))
    def test_empty(self):
        self.assertEqual(global_align("","AC").score,-4)
    def test_limits(self):
        with self.assertRaises(ValueError): global_align("A"*1000,"A"*1000)
        with self.assertRaises(ValueError): global_align("A","A",gap=0)
