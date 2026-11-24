import unittest
import gzip
import tempfile
from pathlib import Path
from io import StringIO
from contextlib import redirect_stdout
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
    def test_gzip_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"example.fastq.gz"
            with gzip.open(path,"wt",encoding="utf-8") as out:
                out.write("@r\nAT\n+\nII\n")
            output = StringIO()
            with redirect_stdout(output):
                self.assertEqual(main(["fastq-stats", str(path)]),0)
            self.assertIn("r\t2\t40.0000\t100.0000\t100.0000",output.getvalue())

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
