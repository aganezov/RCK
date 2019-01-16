import unittest

from rck.core.graph import IntervalAdjacencyGraph
from rck.core.structures import Position, Strand, Segment, Adjacency, AdjacencyType


class TestIntervalAdjacencyGraph(unittest.TestCase):
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

    def test_construction_no_adjacencies(self):
        segments = [self.s1, self.s2, self.s3]
        iag = IntervalAdjacencyGraph(segments=segments, adjacencies=[])
        nodes = list(iag.nodes(data=True))
        self.assertEqual(len(nodes), 6)
        self.assertEqual(len(list(iag.edges())), 3)
        segment_edges = list(iag.segment_edges(data=True))
        self.assertEqual(len(segment_edges), 3)
        self.assertEqual(len(list(iag.adjacency_edges(data=True))), 0)
        segments_edges_as_objects_on_edges = [e[2]["object"] for e in segment_edges]
        for s in segments:
            self.assertIn(s.start_position.idx, {n[0] for n in nodes})
            self.assertIn(s.end_position.idx, {n[0] for n in nodes})
            self.assertIn(s, segments_edges_as_objects_on_edges)

    def test_construction_only_ref_adjacencies(self):
        segments = [self.s1, self.s2, self.s3]
        ra1 = Adjacency(position1=self.s1.end_position, position2=self.s2.start_position, adjacency_type=AdjacencyType.REFERENCE)
        ra2 = Adjacency(position1=self.s2.end_position, position2=self.s3.start_position, adjacency_type=AdjacencyType.REFERENCE)
        adjacencies = [ra1, ra2]
        iag = IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)
        nodes = list(iag.nodes(data=True))
        self.assertEqual(len(nodes), 6)
        edges = list(iag.edges(data=True))
        segment_edges = list(iag.segment_edges(data=True))
        adjacency_edges = list(iag.adjacency_edges(data=True))
        self.assertEqual(len(edges), 5)
        self.assertEqual(len(segment_edges), 3)
        self.assertEqual(len(adjacency_edges), 2)
        for s in segments:
            self.assertIn(s.start_position.idx, {n[0] for n in nodes})
            self.assertIn(s.end_position.idx, {n[0] for n in nodes})
            self.assertIn(s, {e[2]["object"] for e in segment_edges})
        for a in adjacencies:
            self.assertIn(a, {e[2]["object"] for e in adjacency_edges})

    def test_construction_ref_and_nov_adjacencies(self):
        segments = [self.s1, self.s2, self.s3]
        ra1 = Adjacency(position1=self.s1.end_position, position2=self.s2.start_position, adjacency_type=AdjacencyType.REFERENCE)
        ra2 = Adjacency(position1=self.s2.end_position, position2=self.s3.start_position, adjacency_type=AdjacencyType.REFERENCE)
        na1 = Adjacency(position1=self.s1.end_position, position2=self.s3.start_position, adjacency_type=AdjacencyType.NOVEL)
        adjacencies = [ra1, ra2, na1]
        iag = IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)
        nodes = list(iag.nodes(data=True))
        edges = list(iag.edges(data=True))
        segment_edges = list(iag.segment_edges(data=True))
        adjacency_edges = list(iag.adjacency_edges(data=True))
        r_adjacency_edges = list(iag.ref_adjacency_edges(data=True))
        n_adjacency_edges = list(iag.nov_adjacency_edges(data=True))
        self.assertSetEqual({(e[0], e[1]) for e in adjacency_edges},
                            {(e[0], e[1]) for e in r_adjacency_edges}.union({(e[0], e[1]) for e in n_adjacency_edges}))
        self.assertSetEqual({e[2]["object"] for e in segment_edges}, set(segments))
        self.assertEqual(len(nodes), 6)
        self.assertEqual(len(edges), 6)
        self.assertEqual(len(segment_edges), 3)
        self.assertEqual(len(adjacency_edges), 3)
        self.assertEqual(len(r_adjacency_edges), 2)
        self.assertEqual(len(n_adjacency_edges), 1)

    def test_construction_invalid_consistency_check_ref_from_different_chromosomes(self):
        invalid_ra = Adjacency(position1=self.s2.start_position, position2=self.s4.end_position, adjacency_type=AdjacencyType.REFERENCE)
        segments = [self.s1, self.s2, self.s3, self.s4]
        adjacencies = [invalid_ra]
        with self.assertRaises(ValueError):
            IntervalAdjacencyGraph.check_consistency(segments=segments, adjacencies=adjacencies)
        with self.assertRaises(ValueError):
            IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)

    def test_construction_invalid_consistency_check_adjacency_with_position_not_from_segments(self):
        adjacency = Adjacency(position1=self.s4.end_position, position2=self.s1, adjacency_type=AdjacencyType.NOVEL)
        segments = [self.s1, self.s2, self.s3]
        adjacencies = [adjacency]
        with self.assertRaises(ValueError):
            IntervalAdjacencyGraph.check_consistency(segments=segments, adjacencies=adjacencies)
        with self.assertRaises(ValueError):
            IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)
