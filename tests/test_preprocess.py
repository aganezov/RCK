import unittest

from core.preprocess import cluster_positions
from dassp.core.structures import Position, Strand, PositionCluster


class TestPositionsPreprocessor(unittest.TestCase):
    def setUp(self):
        self.p0 = Position(chromosome="chr1", coordinate=1, strand=Strand.REVERSE)
        self.p1 = Position(chromosome="chr1", coordinate=40, strand=Strand.FORWARD)
        self.p2 = Position(chromosome="chr1", coordinate=100, strand=Strand.FORWARD)
        self.p3 = Position(chromosome="chr1", coordinate=200, strand=Strand.FORWARD)
        self.p4 = Position(chromosome="chr2", coordinate=1, strand=Strand.FORWARD)
        self.p5 = Position(chromosome="chr2", coordinate=25, strand=Strand.FORWARD)
        self.p6 = Position(chromosome="chr3", coordinate=1, strand=Strand.FORWARD)
        self.p7 = Position(chromosome="chr4", coordinate=1, strand=Strand.FORWARD)
        self.p8 = Position(chromosome="chr4", coordinate=10, strand=Strand.FORWARD)
        self.p9 = Position(chromosome="chr4", coordinate=1000, strand=Strand.FORWARD)
        self.all_positions = [self.p0, self.p1, self.p2, self.p3, self.p4, self.p5, self.p6, self.p7, self.p8, self.p9]

    def test_positions_clustering_individual_cluster(self):
        clusters = cluster_positions(positions=[self.p3, self.p1, self.p2, self.p0], sigma=50)
        self.assertEqual(len(clusters), 3)
        for c in clusters:
            self.assertIsInstance(c, PositionCluster)
        self.assertListEqual(clusters[0].positions, [self.p0, self.p1])
        self.assertListEqual(clusters[1].positions, [self.p2])
        self.assertListEqual(clusters[2].positions, [self.p3])
        clusters = cluster_positions(positions=[self.p3, self.p1, self.p0, self.p2], sigma=100)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].positions, [self.p0, self.p1, self.p2, self.p3])
        clusters = cluster_positions(positions=[self.p3, self.p1, self.p0, self.p2], sigma=10)
        self.assertEqual(len(clusters), 4)
        self.assertListEqual(clusters[0].positions, [self.p0])
        self.assertListEqual(clusters[1].positions, [self.p1])
        self.assertListEqual(clusters[2].positions, [self.p2])
        self.assertListEqual(clusters[3].positions, [self.p3])

    def test_positions_clustering_multiple_chromosomes(self):
        clusters = cluster_positions(positions=self.all_positions, sigma=25)
        self.assertEqual(len(clusters), 8)
        self.assertListEqual(clusters[0].positions, [self.p0])
        self.assertListEqual(clusters[1].positions, [self.p1])
        self.assertListEqual(clusters[2].positions, [self.p2])
        self.assertListEqual(clusters[3].positions, [self.p3])
        self.assertListEqual(clusters[4].positions, [self.p4, self.p5])
        self.assertListEqual(clusters[5].positions, [self.p6])
        self.assertListEqual(clusters[6].positions, [self.p7, self.p8])
        self.assertListEqual(clusters[7].positions, [self.p9])


if __name__ == '__main__':
    unittest.main()
