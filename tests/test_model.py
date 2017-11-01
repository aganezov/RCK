import unittest

from dassp.core.model import DASSPModel
from dassp.core.structures import Position, Segment, Adjacency, Strand, SegmentCopyNumberRecord, AdjacencyType


class ModelTestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = Position(chromosome="1", coordinate=1, strand=Strand.REVERSE)
        self.p2 = Position(chromosome="1", coordinate=2, strand=Strand.FORWARD)
        self.p3 = Position(chromosome="1", coordinate=3, strand=Strand.REVERSE)
        self.p4 = Position(chromosome="1", coordinate=4, strand=Strand.FORWARD)
        self.p5 = Position(chromosome="1", coordinate=5, strand=Strand.REVERSE)
        self.p6 = Position(chromosome="1", coordinate=6, strand=Strand.FORWARD)

        self.p7 = Position(chromosome="2", coordinate=1, strand=Strand.REVERSE)
        self.p8 = Position(chromosome="2", coordinate=2, strand=Strand.FORWARD)
        self.p9 = Position(chromosome="2", coordinate=3, strand=Strand.REVERSE)
        self.p10 = Position(chromosome="2", coordinate=4, strand=Strand.FORWARD)
        self.p11 = Position(chromosome="2", coordinate=5, strand=Strand.REVERSE)
        self.p12 = Position(chromosome="2", coordinate=6, strand=Strand.FORWARD)

        self.s1 = Segment(start_position=self.p1, end_position=self.p2)
        self.s2 = Segment(start_position=self.p3, end_position=self.p4)
        self.s3 = Segment(start_position=self.p5, end_position=self.p6)

        self.s4 = Segment(start_position=self.p7, end_position=self.p8)
        self.s5 = Segment(start_position=self.p9, end_position=self.p10)
        self.s6 = Segment(start_position=self.p11, end_position=self.p12)

        self.ra1 = Adjacency(position1=self.s1.end_position, position2=self.s2.start_position, adjacency_type=AdjacencyType.REFERENCE)
        self.ra2 = Adjacency(position1=self.s2.end_position, position2=self.s3.start_position, adjacency_type=AdjacencyType.REFERENCE)

    def test_model_creation_not_enough_data(self):
        adjacencies = [self.ra1, self.ra2]
        segments = [self.s1, self.s2, self.s3]
        segment_copy_number_records = [
            SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=1, maj_b_cn=1),
            SegmentCopyNumberRecord(segment=self.s2, maj_a_cn=2, maj_b_cn=2),
            SegmentCopyNumberRecord(segment=self.s3, maj_a_cn=3, maj_b_cn=2)
        ]
        with self.assertRaises(ValueError):
            # missing segments
            DASSPModel(segments=[], adjacencies=adjacencies, scnr=[], n=3, allele_specific=False)
        with self.assertRaises(ValueError):
            # missing segment copy number records
            DASSPModel(segments=segments, adjacencies=adjacencies, scnr=[], n=3, allele_specific=False)
        with self.assertRaises(ValueError):
            # bad number of clones
            DASSPModel(segments=segments, adjacencies=adjacencies, scnr=segment_copy_number_records,
                       n=-1, allele_specific=False)

    def test_model_creation_adjacencies_groups_inconsistent_with_adjacencies(self):
        adjacencies = [self.ra1]
        segments = [self.s1, self.s2, self.s3]
        segment_copy_number_records = [
            SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=1, maj_b_cn=1),
            SegmentCopyNumberRecord(segment=self.s2, maj_a_cn=2, maj_b_cn=2),
            SegmentCopyNumberRecord(segment=self.s3, maj_a_cn=2, maj_b_cn=2)
        ]
        with self.assertRaises(ValueError):
            DASSPModel(segments=segments, adjacencies=adjacencies, scnr=segment_copy_number_records, n=3,
                       adjacencies_groups=[[self.ra1, self.ra2]])

    def test_model_creation_non_allele_specific(self):
        na = Adjacency(position1=self.s2.start_position, position2=self.s2.end_position)
        adjacencies = [self.ra1, self.ra2, na]
        segments = [self.s1, self.s2, self.s3]
        segment_copy_number_records = [
            SegmentCopyNumberRecord(segment=self.s1, maj_a_cn=1, maj_b_cn=0, min_a_cn=2),
            SegmentCopyNumberRecord(segment=self.s2, maj_a_cn=3, maj_b_cn=0, min_a_cn=1),
            SegmentCopyNumberRecord(segment=self.s3, maj_a_cn=1, maj_b_cn=0, min_a_cn=2)
        ]
        model = DASSPModel(segments=segments, adjacencies=adjacencies, scnr=segment_copy_number_records, n=3)
        model.build_gurobi_model()
        # model.solve_gurobi_model()


if __name__ == '__main__':
    unittest.main()
