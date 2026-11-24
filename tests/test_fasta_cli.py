import unittest
from io import StringIO
from contextlib import redirect_stdout, redirect_stderr
from unittest.mock import patch
from biotoolkit import read_fasta
from biotoolkit.cli import main

class FastaTests(unittest.TestCase):
    def test_multiline(self):
        records = list(read_fasta(StringIO(">a description\nac\ngt\n\n>b\nNN\n")))
        self.assertEqual([(r.id,r.description,r.sequence) for r in records], [("a","description","ACGT"),("b","","NN")])
    def test_invalid(self):
        for data in ("", "ACGT\n", ">\nAC\n", ">a\n>b\nAT\n", ">a\nACX\n"):
            with self.subTest(data=data), self.assertRaises(ValueError): list(read_fasta(StringIO(data)))
    def test_duplicate_ids(self):
        self.assertEqual(len(list(read_fasta(StringIO(">a\nA\n>a\nT\n")))), 2)

class CliTests(unittest.TestCase):
    def run_cli(self, args, data):
        out, err = StringIO(), StringIO()
        with patch("sys.stdin", StringIO(data)), redirect_stdout(out), redirect_stderr(err):
            code = main(args)
        return code, out.getvalue(), err.getvalue()
    def test_stats(self):
        code, out, err = self.run_cli(["stats","-"], ">sample\nACGTNN\n>unknown\nNN\n")
        self.assertEqual(code, 0)
        self.assertIn("sample\t6\t1\t1\t1\t1\t2\t50.0000",out)
        self.assertIn("unknown\t2\t0\t0\t0\t0\t2\tNA",out)
    def test_revcomp(self):
        self.assertEqual(self.run_cli(["revcomp","-"], ">a\nATGC\n")[:2], (0,">a\nGCAT\n"))
    def test_translate(self):
        self.assertEqual(self.run_cli(["translate","-"], ">a\nATGGCNTAA\n")[:2], (0,">a\nMA*\n"))
    def test_bad_input(self):
        code, out, err = self.run_cli(["stats","-"], ">a\nACZ\n")
        self.assertEqual(code,1)
        self.assertIn("Invalid DNA",err)
