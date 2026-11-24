import unittest
from itertools import product
from biotoolkit import gc_content, reverse_complement, translate, sequence_stats
from biotoolkit.sequence import normalize_dna

class SequenceTests(unittest.TestCase):
    def test_normalize(self):
        self.assertEqual(normalize_dna("a c\nGt n"), "ACGTN")
    def test_invalid(self):
        for seq in ("ACGU", "AC-G", "123", "?"):
            with self.subTest(seq=seq), self.assertRaises(ValueError):
                normalize_dna(seq)
    def test_gc(self):
        self.assertEqual(gc_content("GGACTN"), 60.0)
        self.assertIsNone(gc_content("NNRY"))
        self.assertIsNone(gc_content(""))
    def test_stats(self):
        self.assertEqual(sequence_stats("ACGTNR"), {"length":6,"A":1,"C":1,"G":1,"T":1,"ambiguous":2,"gc_percent":50.0})
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGCRY"), "RYGCAT")
        self.assertEqual(reverse_complement("RYSWKMBDHVN"), "NBDHVKMWSRY")
        seq = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(reverse_complement(seq)), seq)
    def test_translation(self):
        self.assertEqual(translate("ATGGCTTGATAA"), "MA**")
        self.assertEqual(translate("AATGGCT", frame=1), "MA")
        self.assertEqual(translate("AAATGGCT", frame=2), "MA")
    def test_ambiguity(self):
        self.assertEqual(translate("GCNTGYTARNNN"), "AC*X")
    def test_partial(self):
        with self.assertRaises(ValueError): translate("ATGA")
        self.assertEqual(translate("ATGA", partial="trim"), "M")
        with self.assertRaises(ValueError): translate("ATG", 3)
        with self.assertRaises(ValueError): translate("ATG", partial="wrong")
    def test_all_codons(self):
        # Independent standard-code table grouped by residue.
        groups = {"F":"TTT TTC", "L":"TTA TTG CTT CTC CTA CTG", "I":"ATT ATC ATA", "M":"ATG", "V":"GTT GTC GTA GTG", "S":"TCT TCC TCA TCG AGT AGC", "P":"CCT CCC CCA CCG", "T":"ACT ACC ACA ACG", "A":"GCT GCC GCA GCG", "Y":"TAT TAC", "*":"TAA TAG TGA", "H":"CAT CAC", "Q":"CAA CAG", "N":"AAT AAC", "K":"AAA AAG", "D":"GAT GAC", "E":"GAA GAG", "C":"TGT TGC", "W":"TGG", "R":"CGT CGC CGA CGG AGA AGG", "G":"GGT GGC GGA GGG"}
        expected = {codon: aa for aa, codons in groups.items() for codon in codons.split()}
        self.assertEqual(len(expected), 64)
        for codon in product("ACGT", repeat=3):
            c = "".join(codon)
            self.assertEqual(translate(c), expected[c], c)
