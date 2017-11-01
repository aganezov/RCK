import os
import unittest
from dassp.core.io import PCAWGVCFNovelAdjacencyReader, BattenbergSegmentReader
from dassp.core.structures import Position, Strand


class PCAWGVCFNovelAdjacencyReaderTestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = Position(chromosome="chr1", coordinate=1, strand=Strand.FORWARD)
        self.p2 = Position(chromosome="chr2", coordinate=2, strand=Strand.FORWARD)
        self.p3 = Position(chromosome="chr2", coordinate=3, strand=Strand.FORWARD)
        self.p4 = Position(chromosome="chr3", coordinate=2, strand=Strand.FORWARD)
        self.p1.extra["self_id"] = "SVMERGE1_1"
        self.p2.extra["self_id"] = "SVMERGE1_2"
        self.p3.extra["self_id"] = "SVMERGE2_1"
        self.p4.extra["self_id"] = "SVMERGE2_2"
        self.p1.extra["mate_id"] = "SVMERGE1_2"
        self.p2.extra["mate_id"] = "SVMERGE1_1"
        self.p3.extra["mate_id"] = "SVMERGE2_2"
        self.p4.extra["mate_id"] = "SVMERGE2_1"
        self.reader = PCAWGVCFNovelAdjacencyReader(file_path=os.path.join("..", "test_data", "0009b464-b376-4fbc-8a56-da538269a02f.pcawg_consensus_1.6.161116.somatic.sv.vcf"))

    def test_na_from_positions(self):
        na_positions = [self.p1, self.p2, self.p3, self.p4]
        adjacencies = PCAWGVCFNovelAdjacencyReader.get_adjacencies_from_positions(na_positions=na_positions)
        self.assertEqual(len(adjacencies), 2)
        self.assertSetEqual({a.idx for a in adjacencies}, {"1", "2"})

    def test_na_from_positions_value_errors(self):
        self.p1.extra["mate_id"] = "SVMERGE2_1"
        with self.assertRaises(ValueError):
            PCAWGVCFNovelAdjacencyReader.get_adjacencies_from_positions(na_positions=[self.p1, self.p2, self.p3, self.p4])
        self.p1.extra["mate_id"] = "SVMERGE1_2"
        self.p1.extra["self_id"] = "SVMERGE1_11"
        with self.assertRaises(ValueError):
            PCAWGVCFNovelAdjacencyReader.get_adjacencies_from_positions(na_positions=[self.p1, self.p2, self.p3, self.p4])

    def test_file_parse_pcawg_position_from_test_data(self):
        positions = self.reader.parse_pcawg_vcf_position_file()
        self.assertEqual(len(positions), 666)
        adjacencies = PCAWGVCFNovelAdjacencyReader.get_adjacencies_from_positions(na_positions=positions)
        self.assertEqual(len(adjacencies), 333)


class BattenbergSegmentReaderTestCase(unittest.TestCase):
    def setUp(self):
        self.battenberg_string_record_1 = "1	1	783071	17373317	0.503606020372035	1	0.260389636832818	3.93085389414586	2	2	1	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA"
        self.battenberg_string_record_2 = "2	1	17373751	27264269	0.575606746031517	0	0.0543201014005013	3.38013749740953	2	1	0.552546527501979	2	2	0.447453472498021	0.0070476350401058	0.006862365350927	0.538663907729184	0.565668929691937	2	0	0.276273263750989	2	2	0.723726736249011	0.00352381752005285	0.0037356565796258	0.26913462359894	0.283778618192779	2	1	0.776273263750989	2	3	0.223726736249011	0.00352381752005293	0.0037356565796258	0.76913462359894	0.783778618192779	1	1	0.606883438642538	2	1	0.393116561357462	0.00680071620748574	0.00720896748029918	0.592181279780749	0.620357765306304	0	1	0.303441719321269	2	1	0.696558280678731	0.00340035810374287	0.00330141329110712	0.296897586506694	0.309811306162619	1	1	0.803441719321269	3	1	0.196558280678731	0.00340035810374256	0.00330141329110714	0.796897586506693	0.809811306162619"

    def test_parse_bsr1(self):
        bpr1 = BattenbergSegmentReader.parse(battenberg_record=self.battenberg_string_record_1)
        self.assertEqual(bpr1.extra["battenberg_id"], "1")
        self.assertAlmostEqual(bpr1.extra["maj_frac"], float(1), places=4)
        self.assertAlmostEqual(bpr1.extra["min_frac"], float(0), places=4)
        self.assertEqual(bpr1.segment.start_position.coordinate, 783071)
        self.assertEqual(bpr1.segment.end_position.coordinate, 17373317)
        self.assertEqual(bpr1.maj_a_segmentcn.segment, bpr1.segment)
        self.assertEqual(bpr1.maj_b_segmentcn.segment, bpr1.segment)
        self.assertEqual(bpr1.min_a_segmentcn.segment, bpr1.segment)
        self.assertEqual(bpr1.min_b_segmentcn.segment, bpr1.segment)
        self.assertEqual(bpr1.maj_a_segmentcn.cn, 2)
        self.assertEqual(bpr1.maj_b_segmentcn.cn, 2)
        self.assertEqual(bpr1.min_a_segmentcn.cn, 2)
        self.assertEqual(bpr1.min_b_segmentcn.cn, 2)
        self.assertFalse(bpr1.distinct_min_cn)

    def test_parse_bsr2(self):
        bpr2 = BattenbergSegmentReader.parse(battenberg_record=self.battenberg_string_record_2)
        self.assertEqual(bpr2.extra["battenberg_id"], "2")
        self.assertAlmostEqual(bpr2.extra["maj_frac"], float(0.5525), places=4)
        self.assertAlmostEqual(bpr2.extra["min_frac"], float(0.4475), places=4)
        self.assertEqual(bpr2.segment.start_position.coordinate, 17373751)
        self.assertEqual(bpr2.segment.start_position.strand, Strand.REVERSE)
        self.assertEqual(bpr2.segment.start_position.chromosome, "1")
        self.assertEqual(bpr2.segment.end_position.coordinate, 27264269)
        self.assertEqual(bpr2.segment.end_position.strand, Strand.FORWARD)
        self.assertEqual(bpr2.segment.end_position.chromosome, "1")
        self.assertEqual(bpr2.maj_a_segmentcn.segment, bpr2.segment)
        self.assertEqual(bpr2.maj_b_segmentcn.segment, bpr2.segment)
        self.assertEqual(bpr2.min_a_segmentcn.segment, bpr2.segment)
        self.assertEqual(bpr2.min_b_segmentcn.segment, bpr2.segment)
        self.assertEqual(bpr2.maj_a_segmentcn.cn, 2)
        self.assertEqual(bpr2.maj_b_segmentcn.cn, 1)
        self.assertEqual(bpr2.min_a_segmentcn.cn, 2)
        self.assertEqual(bpr2.min_b_segmentcn.cn, 2)
        self.assertTrue(bpr2.distinct_min_cn)

    def test_file_parse_battenberg_test_data(self):
        parser = BattenbergSegmentReader(file_path=os.path.join("..", "test_data", "0009b464-b376-4fbc-8a56-da538269a02f_subclones.txt"))
        records = parser.parse_battenberg_file()
        self.assertEqual(len(records), 251)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "1"]), 11)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "2"]), 9)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "3"]), 21)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "4"]), 12)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "5"]), 9)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "6"]), 18)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "7"]), 13)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "8"]), 30)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "9"]), 6)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "10"]), 11)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "11"]), 13)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "12"]), 15)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "13"]), 7)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "14"]), 5)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "15"]), 1)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "16"]), 11)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "17"]), 7)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "18"]), 4)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "19"]), 19)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "20"]), 10)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "21"]), 8)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "22"]), 1)
        self.assertEqual(len([r for r in records if r.segment.chromosome == "X"]), 10)


if __name__ == '__main__':
    unittest.main()
