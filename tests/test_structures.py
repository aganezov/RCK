import unittest

from core.structures import PositionCluster, SegmentCopyNumberRecord
from dassp.core.structures import Strand, Position, Segment, Adjacency


class StrandTestCase(unittest.TestCase):
    def test_strand_str(self):
        self.assertEqual(str(Strand.REVERSE), "-")
        self.assertEqual(str(Strand.FORWARD), "+")

    def test_from_pm_string(self):
        self.assertEqual(Strand.REVERSE, Strand.from_pm_string(string="-"))
        self.assertEqual(Strand.FORWARD, Strand.from_pm_string(string="+"))
        with self.assertRaises(ValueError):
            Strand.from_pm_string(string="?")


class PositionTestCase(unittest.TestCase):
    def setUp(self):
        self.position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.FORWARD)
        self.position2 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.position3 = Position(chromosome="chr1", coordinate=2, strand=Strand.FORWARD)
        self.position4 = Position(chromosome="chr2", coordinate=1, strand=Strand.FORWARD)

    def test_empty_extra_creation(self):
        self.assertDictEqual(Position(chromosome="chrom1", coordinate=1, strand=Strand.FORWARD).extra, {})

    def test_eq(self):
        non_position = "?"
        eq_position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.FORWARD)
        self.assertNotEqual(self.position1, self.position2)
        self.assertNotEqual(self.position1, self.position3)
        self.assertNotEqual(self.position1, self.position4)
        self.assertNotEqual(self.position1, non_position)
        self.assertEqual(self.position1, eq_position1)

    def test_lt(self):
        self.assertLess(self.position2, self.position1)
        self.assertGreater(self.position1, self.position2)
        self.assertLess(self.position1, self.position3)
        chr5_position = Position(chromosome="chr5", coordinate=5, strand=Strand.FORWARD)
        chr10_position = Position(chromosome="chr10", coordinate=1, strand=Strand.REVERSE)
        self.assertLess(self.position1, self.position4)
        self.assertLess(chr5_position, chr10_position)


class SegmentTestCase(unittest.TestCase):
    def setUp(self):
        self.position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.position2 = Position(chromosome="chr1", coordinate=2, strand=Strand.FORWARD)

    def test_creation(self):
        position3 = Position(chromosome="chr1", coordinate=1, strand=Strand.FORWARD)
        position4 = Position(chromosome="chr2", coordinate=2, strand=Strand.FORWARD)
        position5 = Position(chromosome="chr1", coordinate=0, strand=Strand.FORWARD)
        for pos in [position3, position4, position5]:
            with self.assertRaises(ValueError):
                Segment(start_position=self.position1, end_position=pos)
        Segment(start_position=self.position1, end_position=self.position2)

    def test_idx(self):
        s = Segment(start_position=self.position1, end_position=self.position2)
        self.assertIsNone(s._idx)
        self.assertEqual(s.idx, "chr1:1-2")
        s = Segment(start_position=self.position1, end_position=self.position2, idx="idx")
        self.assertEqual(s.idx, "idx")
        s.idx = "idx2"
        self.assertEqual(s.idx, "idx2")

    def test_str(self):
        s = Segment(start_position=self.position1, end_position=self.position2)
        self.assertEqual(str(s), "chr1:1-2")
        s.idx = "idx"
        self.assertEqual(str(s), "idx")

    def test_chromosome(self):
        s = Segment(start_position=self.position1, end_position=self.position2)
        self.assertEqual(s.chromosome, self.position1.chromosome)
        self.assertEqual(s.chromosome, self.position2.chromosome)


class AdjacencyTestCase(unittest.TestCase):
    def setUp(self):
        self.position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.position2 = Position(chromosome="chr1", coordinate=2, strand=Strand.FORWARD)

    def test_creation(self):
        a = Adjacency(position1=self.position2, position2=self.position1)
        self.assertEqual(a.position1, self.position1)
        self.assertEqual(a.position2, self.position2)
        a = Adjacency(position1=self.position1, position2=self.position2)
        self.assertEqual(a.position1, self.position1)
        self.assertEqual(a.position2, self.position2)

    def test_idx(self):
        s = Adjacency(position1=self.position1, position2=self.position2, idx="idx")
        self.assertEqual(s.idx, "idx")
        s.idx = None
        self.assertEqual(s.idx, "[" + str(self.position1) + "]-[" + str(self.position2) + "]")


class PositionClusterTestCase(unittest.TestCase):
    def setUp(self):
        self.position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.FORWARD)
        self.position2 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.position3 = Position(chromosome="chr1", coordinate=2, strand=Strand.FORWARD)
        self.position4 = Position(chromosome="chr2", coordinate=1, strand=Strand.FORWARD)

    def test_interning_sorting_on_creation(self):
        pc = PositionCluster(positions=[self.position3, self.position1, self.position2])
        self.assertListEqual(pc.positions, [self.position2, self.position1, self.position3])


class SegmentCopyNumberRecordTestCase(unittest.TestCase):
    def setUp(self):
        self.position1 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.position2 = Position(chromosome="chr1", coordinate=2, strand=Strand.FORWARD)

        self.s1 = Segment(start_position=self.position1, end_position=self.position2)

    def test_transform_to_non_allele_specific(self):
        maj_a_cn = 1
        maj_b_cn = 1
        scnr = SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=maj_a_cn, maj_b_cn=maj_b_cn)
        nas_cnr = SegmentCopyNumberRecord.to_non_allele_specific(scnr=scnr)
        self.assertEqual(nas_cnr.maj_a_segmentcn.segment, self.s1)
        self.assertFalse(nas_cnr.allele_specific)
        self.assertEqual(nas_cnr.maj_b_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.min_a_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.min_b_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.maj_a_segmentcn.cn, maj_a_cn + maj_b_cn)
        self.assertEqual(nas_cnr.maj_b_segmentcn.cn, 0)
        self.assertEqual(nas_cnr.min_a_segmentcn.cn, maj_a_cn + maj_b_cn)
        self.assertEqual(nas_cnr.min_b_segmentcn.cn, 0)
        self.assertFalse(nas_cnr.distinct_min_cn)

    def test_to_non_allele_specific_distinct_minority(self):
        maj_a_cn = 2
        maj_b_cn = 1
        min_a_cn = 1
        min_b_cn = 1
        scnr = SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=maj_a_cn, maj_b_cn=maj_b_cn, min_a_cn=min_a_cn, min_b_cn=min_b_cn)
        nas_cnr = SegmentCopyNumberRecord.to_non_allele_specific(scnr=scnr)
        self.assertFalse(nas_cnr.allele_specific)
        self.assertEqual(nas_cnr.maj_a_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.maj_b_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.min_a_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.min_b_segmentcn.segment, self.s1)
        self.assertEqual(nas_cnr.maj_a_segmentcn.cn, maj_a_cn + maj_b_cn)
        self.assertEqual(nas_cnr.maj_b_segmentcn.cn, 0)
        self.assertEqual(nas_cnr.min_a_segmentcn.cn, min_a_cn + min_b_cn)
        self.assertEqual(nas_cnr.min_b_segmentcn.cn, 0)
        self.assertTrue(nas_cnr.distinct_min_cn)

    def test_to_non_allele_specific_non_allele_specific(self):
        maj_a_cn = 2
        maj_b_cn = 1
        min_a_cn = 1
        min_b_cn = 1
        scnr = SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=maj_a_cn, maj_b_cn=maj_b_cn, min_a_cn=min_a_cn, min_b_cn=min_b_cn)
        nas_cnr = SegmentCopyNumberRecord.to_non_allele_specific(scnr=scnr)
        self.assertFalse(nas_cnr.allele_specific)
        nas_cnr2 = SegmentCopyNumberRecord.to_non_allele_specific(scnr=nas_cnr)
        self.assertEqual(nas_cnr, nas_cnr2)
        self.assertFalse(nas_cnr2.allele_specific)


if __name__ == '__main__':
    unittest.main()
