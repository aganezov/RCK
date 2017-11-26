import unittest

from dassp.core.structures import Haplotype, Segment, Position, Strand, PositionType, AdjacencyType
from dassp.simulation.parts import ChromosomeGenerator, reverse_segment, delete_segments, duplicate_segments, translocation_segments, reverse_segments
from simulation.manager import SimulationManager, generate_mutated_genome, get_adjacencies_from_genome, get_scn_profile_from_genome


class ChromosomeGeneratorTestCase(unittest.TestCase):
    def test_generate_single_mp_chromosome_pair(self):
        size = 10
        chr_cnt = 1
        chromosomes = ChromosomeGenerator.generate_chromosomes(chromosome_size=10, chromosomes_cnt=1, ab=True)
        self.assertEqual(len(chromosomes), chr_cnt * 2)
        self.assertEqual(len(chromosomes[0]), size)
        self.assertEqual(len(chromosomes[1]), size)
        for chr_a_s, chr_b_s in zip(chromosomes[0], chromosomes[1]):
            self.assertIsInstance(chr_a_s, Segment)
            self.assertIsInstance(chr_b_s, Segment)
            self.assertEqual(chr_a_s.idx, chr_b_s.idx)
            self.assertEqual(chr_a_s.extra["haplotype"], Haplotype.A)
            self.assertEqual(chr_b_s.extra["haplotype"], Haplotype.B)

    def test_generate_multiple_chromosome_pairs(self):
        size = 100
        for chr_cnt in [1, 5, 10, 22]:
            chromosomes = ChromosomeGenerator.generate_chromosomes(chromosome_size=size, chromosomes_cnt=chr_cnt, ab=True)
            self.assertEqual(len(chromosomes), 2 * chr_cnt)
            for chr_a, chr_b in zip(chromosomes[::2], chromosomes[1::2]):
                self.assertEqual(len(chr_a), size)
                self.assertEqual(len(chr_b), size)
                for chr_a_s, chr_b_s in zip(chr_a, chr_b):
                    self.assertIsInstance(chr_a_s, Segment)
                    self.assertIsInstance(chr_b_s, Segment)
                    self.assertEqual(chr_a_s.idx, chr_b_s.idx)
                    self.assertIsNot(chr_a_s, chr_b_s)
                    self.assertEqual(chr_a_s.extra["haplotype"], Haplotype.A)
                    self.assertEqual(chr_b_s.extra["haplotype"], Haplotype.B)


class OperationTestCase(unittest.TestCase):
    def setUp(self):
        self.chr_size = 100
        self.chromosome = ChromosomeGenerator.generate_chromosome(chromosome_size=self.chr_size, chr_name="1")

        self.spans = [1, 5, 10, 15, 20, 25, 30, 50, 70, 80]

    def test_reverse_segment(self):
        p1 = Position(chromosome="1", coordinate=1, strand=Strand.REVERSE, ptype=PositionType.ARTIFICIAL)
        p2 = Position(chromosome="1", coordinate=2, strand=Strand.FORWARD, ptype=PositionType.ARTIFICIAL)
        s = Segment(start_position=p1, end_position=p2)
        reversed_segment = reverse_segment(segment=s)
        self.assertIs(s, reversed_segment)
        self.assertIs(reversed_segment.start_position, p2)
        self.assertIs(reversed_segment.end_position, p1)


class DeletionTestCase(OperationTestCase):
    def test_deletion_size_one(self):
        for i in range(self.chr_size - 1):
            del_start_index = i
            del_end_index = i + 1
            new_chromosome = delete_segments(chromosome=self.chromosome, start_segment_index=del_start_index, end_segment_index=del_end_index)
            self.assertEqual(len(new_chromosome), len(self.chromosome) - 1)
            self.assertListEqual([s.idx for s in new_chromosome], [s.idx for s in self.chromosome[:del_start_index] + self.chromosome[del_end_index:]])

    def test_deletion_various_sizes(self):
        for del_size in self.spans:
            del_start_index = 2
            del_end_index = del_start_index + del_size
            new_chromosome = delete_segments(chromosome=self.chromosome, start_segment_index=del_start_index, end_segment_index=del_end_index)
            self.assertEqual(len(new_chromosome), len(self.chromosome) - del_size)

    def test_arm_deletion_beginning(self):
        new_chromosome = delete_segments(chromosome=self.chromosome, start_segment_index=0, end_segment_index=20)
        self.assertEqual(len(new_chromosome), self.chr_size - 20)
        self.assertListEqual([s.idx for s in new_chromosome], [s.idx for s in self.chromosome[20:]])

    def test_arm_deletion_ending(self):
        new_chromosome = delete_segments(chromosome=self.chromosome, start_segment_index=self.chr_size - 20, end_segment_index=self.chr_size)
        self.assertEqual(len(new_chromosome), self.chr_size - 20)
        self.assertListEqual([s.idx for s in new_chromosome], [s.idx for s in self.chromosome[:-20]])


class DuplicationTestCase(OperationTestCase):
    def test_duplication_single_segment(self):
        for i in range(len(self.chromosome)):
            new_chromosome = duplicate_segments(chromosome=self.chromosome, start_segment_index=i, end_segment_index=i + 1)
            self.assertEqual(len(new_chromosome), len(self.chromosome) + 1)
            self.assertEqual(new_chromosome[i].idx, new_chromosome[i + 1].idx)
            self.assertListEqual([s.idx for s in self.chromosome[:i]], [s.idx for s in new_chromosome[:i]])
            self.assertListEqual([s.idx for s in self.chromosome[i:]], [s.idx for s in new_chromosome[i + 1:]])

    def test_duplication_span_of_segments(self):
        for dupl_span in self.spans:
            dupl_start_index = 2
            dupl_end_index = dupl_start_index + dupl_span
            new_chromosome = duplicate_segments(chromosome=self.chromosome, start_segment_index=dupl_start_index, end_segment_index=dupl_end_index)
            self.assertEqual(len(new_chromosome), len(self.chromosome) + dupl_span)
            self.assertListEqual([s.idx for s in self.chromosome[:dupl_start_index]],
                                 [s.idx for s in new_chromosome[:dupl_start_index]])
            self.assertListEqual([s.idx for s in new_chromosome[dupl_start_index:dupl_start_index + dupl_span]],
                                 [s.idx for s in new_chromosome[dupl_end_index:dupl_end_index + dupl_span]])
            self.assertListEqual([s.idx for s in self.chromosome[dupl_end_index:]],
                                 [s.idx for s in new_chromosome[dupl_end_index + dupl_span:]])

    def test_duplication_arm(self):
        new_chromosome = duplicate_segments(chromosome=self.chromosome, start_segment_index=0, end_segment_index=20)
        self.assertEqual(len(new_chromosome), len(self.chromosome) + 20)
        self.assertListEqual([s.idx for s in new_chromosome[:20]], [s.idx for s in new_chromosome[20:40]])
        self.assertListEqual([s.idx for s in self.chromosome[20:]], [s.idx for s in new_chromosome[40:]])


class TranslocationTestCase(OperationTestCase):
    def setUp(self):
        super(TranslocationTestCase, self).setUp()
        self.chromosome2 = ChromosomeGenerator.generate_chromosome(chromosome_size=self.chr_size, chr_name="2")

    def test_translocation_cc(self):
        for chr1_segment_index in self.spans:
            for chr2_segment_index in self.spans:
                new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                            chromosome1_segment_index=chr1_segment_index, chromosome2_segment_index=chr2_segment_index,
                                                            cc=True)
                self.assertEqual(len(new_chr1), chr1_segment_index + len(self.chromosome2) - chr2_segment_index)
                self.assertEqual(len(new_chr2), chr2_segment_index + len(self.chromosome) - chr1_segment_index)
                self.assertListEqual([s.idx for s in self.chromosome[:chr1_segment_index]],
                                     [s.idx for s in new_chr1[:chr1_segment_index]])
                self.assertListEqual([s.idx for s in self.chromosome2[chr2_segment_index:]],
                                     [s.idx for s in new_chr1[chr1_segment_index:]])
                self.assertListEqual([s.idx for s in self.chromosome2[:chr2_segment_index]],
                                     [s.idx for s in new_chr2[:chr2_segment_index]])
                self.assertListEqual([s.idx for s in self.chromosome[chr1_segment_index:]],
                                     [s.idx for s in new_chr2[chr2_segment_index:]])

    def test_translocation_non_cc(self):
        for chr1_segment_index in self.spans:
            for chr2_segment_index in self.spans:
                new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                            chromosome1_segment_index=chr1_segment_index, chromosome2_segment_index=chr2_segment_index,
                                                            cc=False)
                self.assertEqual(len(new_chr1), chr1_segment_index + chr2_segment_index)
                self.assertEqual(len(new_chr2), len(self.chromosome2) - chr2_segment_index + len(self.chromosome) - chr1_segment_index)
                self.assertListEqual([s.idx for s in self.chromosome[:chr1_segment_index]],
                                     [s.idx for s in new_chr1[:chr1_segment_index]])
                self.assertListEqual([reverse_segment(s).idx for s in self.chromosome2[:chr2_segment_index][::-1]],
                                     [s.idx for s in new_chr1[chr1_segment_index:]])
                self.assertListEqual([reverse_segment(s).idx for s in self.chromosome[chr1_segment_index:][::-1]],
                                     [s.idx for s in new_chr2[:len(self.chromosome) - chr1_segment_index]])
                self.assertListEqual([s.idx for s in self.chromosome2[chr2_segment_index:]],
                                     [s.idx for s in new_chr2[len(self.chromosome) - chr1_segment_index:]])

    def test_translocation_one_arm_empty(self):
        for chr2_segment_index in self.spans:
            new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                        chromosome1_segment_index=0, chromosome2_segment_index=chr2_segment_index, cc=True)
            self.assertEqual(len(new_chr1), len(self.chromosome2) - chr2_segment_index)
            self.assertEqual(len(new_chr2), len(self.chromosome) + chr2_segment_index)
            self.assertListEqual([s.idx for s in new_chr1], [s.idx for s in self.chromosome2[chr2_segment_index:]])
            self.assertListEqual([s.idx for s in new_chr2], [s.idx for s in self.chromosome2[:chr2_segment_index] + self.chromosome])

            new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                        chromosome1_segment_index=0, chromosome2_segment_index=chr2_segment_index, cc=False)
            self.assertEqual(len(new_chr1), chr2_segment_index)
            self.assertEqual(len(new_chr2), len(self.chromosome) + len(self.chromosome2) - chr2_segment_index)
            self.assertListEqual([s.idx for s in new_chr1], [reverse_segment(s).idx for s in self.chromosome2[:chr2_segment_index][::-1]])
            self.assertListEqual([s.idx for s in new_chr2[:len(self.chromosome)]], [reverse_segment(s).idx for s in self.chromosome][::-1])
            self.assertListEqual([s.idx for s in new_chr2[len(self.chromosome):]], [s.idx for s in self.chromosome2[chr2_segment_index:]])

        for chr1_segment_index in self.spans:
            new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                        chromosome1_segment_index=chr1_segment_index, chromosome2_segment_index=0, cc=True)
            self.assertEqual(len(new_chr1), chr1_segment_index + len(self.chromosome2))
            self.assertEqual(len(new_chr2), len(self.chromosome) - chr1_segment_index)
            self.assertListEqual([s.idx for s in new_chr1], [s.idx for s in self.chromosome[:chr1_segment_index] + self.chromosome2])
            self.assertListEqual([s.idx for s in new_chr2], [s.idx for s in self.chromosome[chr1_segment_index:]])
            new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                        chromosome1_segment_index=chr1_segment_index, chromosome2_segment_index=0, cc=False)
            self.assertEqual(len(new_chr1), chr1_segment_index)
            self.assertEqual(len(new_chr2), len(self.chromosome) + len(self.chromosome2) - chr1_segment_index)
            self.assertListEqual([s.idx for s in self.chromosome[:chr1_segment_index]],
                                 [s.idx for s in new_chr1])
            self.assertListEqual([reverse_segment(s).idx for s in self.chromosome[chr1_segment_index:][::-1]],
                                 [s.idx for s in new_chr2[:len(self.chromosome) - chr1_segment_index]])
            self.assertListEqual([s.idx for s in self.chromosome2], [s.idx for s in new_chr2[len(self.chromosome) - chr1_segment_index:]])

    def test_translocation_both_arm_empty(self):
        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=0, chromosome2_segment_index=0, cc=True)
        self.assertEqual(len(new_chr1), len(self.chromosome2))
        self.assertEqual(len(new_chr2), len(self.chromosome))
        self.assertListEqual([s.idx for s in self.chromosome], [s.idx for s in new_chr2])
        self.assertListEqual([s.idx for s in self.chromosome2], [s.idx for s in new_chr1])

        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=0, chromosome2_segment_index=0, cc=False)
        self.assertEqual(len(new_chr1), 0)
        self.assertEqual(len(new_chr2), len(self.chromosome) + len(self.chromosome2))
        self.assertListEqual([reverse_segment(s).idx for s in self.chromosome[::-1]], [s.idx for s in new_chr2[:len(self.chromosome)]])
        self.assertListEqual([s.idx for s in self.chromosome2], [s.idx for s in new_chr2[len(self.chromosome):]])

        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=len(self.chromosome), chromosome2_segment_index=0, cc=True)
        self.assertEqual(len(new_chr1), len(self.chromosome) + len(self.chromosome2))
        self.assertEqual(len(new_chr2), 0)
        self.assertListEqual([s.idx for s in self.chromosome + self.chromosome2], [s.idx for s in new_chr1])

        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=len(self.chromosome), chromosome2_segment_index=0, cc=False)
        self.assertEqual(len(new_chr1), len(self.chromosome))
        self.assertEqual(len(new_chr2), len(self.chromosome2))
        self.assertListEqual([s.idx for s in self.chromosome], [s.idx for s in new_chr1])
        self.assertListEqual([s.idx for s in self.chromosome2], [s.idx for s in new_chr2])
        #
        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=0, chromosome2_segment_index=len(self.chromosome2), cc=True)
        self.assertEqual(len(new_chr2), len(self.chromosome) + len(self.chromosome2))
        self.assertEqual(len(new_chr1), 0)
        self.assertListEqual([s.idx for s in self.chromosome2 + self.chromosome], [s.idx for s in new_chr2])

        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=0, chromosome2_segment_index=len(self.chromosome2), cc=False)
        self.assertEqual(len(new_chr1), len(self.chromosome2))
        self.assertEqual(len(new_chr2), len(self.chromosome))
        self.assertListEqual([reverse_segment(s).idx for s in self.chromosome2[::-1]], [s.idx for s in new_chr1])
        self.assertListEqual([reverse_segment(s).idx for s in self.chromosome[::-1]], [s.idx for s in new_chr2])
        ##
        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=len(self.chromosome), chromosome2_segment_index=len(self.chromosome2), cc=True)
        self.assertEqual(len(new_chr1), len(self.chromosome))
        self.assertEqual(len(new_chr1), len(self.chromosome2))
        self.assertListEqual([s.idx for s in self.chromosome], [s.idx for s in new_chr1])
        self.assertListEqual([s.idx for s in self.chromosome2], [s.idx for s in new_chr2])

        new_chr1, new_chr2 = translocation_segments(chromosome1=self.chromosome, chromosome2=self.chromosome2,
                                                    chromosome1_segment_index=len(self.chromosome), chromosome2_segment_index=len(self.chromosome2), cc=False)
        self.assertEqual(len(new_chr1), len(self.chromosome2) + len(self.chromosome))
        self.assertEqual(len(new_chr2), 0)
        self.assertListEqual([s.idx for s in self.chromosome], [s.idx for s in new_chr1[:len(self.chromosome)]])
        self.assertListEqual([reverse_segment(s).idx for s in self.chromosome2[::-1]], [s.idx for s in new_chr1[len(self.chromosome):]])


class ReversalTestCase(OperationTestCase):
    def test_reversal_single_segment(self):
        for i in range(len(self.chromosome) - 1):
            new_chr = reverse_segments(chromosome=self.chromosome, start_segment_index=i, end_segment_index=i+1)
            self.assertEqual(len(self.chromosome), len(new_chr))
            self.assertEqual(reverse_segment(self.chromosome[i]).idx, new_chr[i].idx)
            self.assertListEqual([s.idx for s in self.chromosome[:i]], [s.idx for s in new_chr[:i]])
            self.assertListEqual([s.idx for s in self.chromosome[i+1:]], [s.idx for s in new_chr[i+1:]])

    def test_reversal_span_of_segments(self):
        for span in self.spans:
            rev_start_index = 2
            rev_end_index = rev_start_index + span
            new_chr = reverse_segments(chromosome=self.chromosome, start_segment_index=rev_start_index, end_segment_index=rev_end_index)
            self.assertEqual(len(self.chromosome), len(new_chr))
            self.assertListEqual([s.idx for s in self.chromosome[:rev_start_index]], [s.idx for s in new_chr[:rev_start_index]])
            self.assertListEqual([s.idx for s in self.chromosome[rev_end_index:]], [s.idx for s in new_chr[rev_end_index:]])
            self.assertListEqual([reverse_segment(s).idx for s in self.chromosome[rev_start_index:rev_end_index][::-1]],
                                 [s.idx for s in new_chr[rev_start_index:rev_end_index]])

    def test_full_chromosomal_reversal(self):
        new_chr = reverse_segments(chromosome=self.chromosome, start_segment_index=0, end_segment_index=len(self.chromosome))
        self.assertListEqual([s.idx for s in new_chr], [reverse_segment(s).idx for s in self.chromosome[::-1]])


class SimulationManagerTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_generate_mutated_genome_3_ab_chromosomes_no_translocation_in_the_mutation_config(self):
        self.genome = ChromosomeGenerator.generate_chromosomes(chromosomes_cnt=3, ab=True, chromosome_size=100)
        history = generate_mutated_genome(starting_genome=self.genome, mutation_cnt=9)
        self.assertEqual(len(history["genomes"]), 10)
        self.assertEqual(len(history["mutations"]), 9)
        for genome in history["genomes"]:
            for chromosome in genome:
                self.assertTrue(len(chromosome) > 0)

    def test_get_adjacencies_from_genome_reference(self):
        genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=100, chromosomes_cnt=3, ab=True)
        adjacencies = get_adjacencies_from_genome(genome=genome, is_reference=True)
        self.assertEqual(len(adjacencies), 297)
        for adjacency_idx in adjacencies:
            self.assertEqual(len(adjacencies[adjacency_idx]), 2)
            self.assertIn((Haplotype.A, Haplotype.A), adjacencies[adjacency_idx])
            self.assertIn((Haplotype.B, Haplotype.B), adjacencies[adjacency_idx])
            self.assertEqual(len(adjacencies[adjacency_idx][(Haplotype.A, Haplotype.A)]), 1)
            self.assertEqual(len(adjacencies[adjacency_idx][(Haplotype.B, Haplotype.B)]), 1)
            self.assertEqual(adjacencies[adjacency_idx][(Haplotype.A, Haplotype.A)][0].adjacency_type,
                             AdjacencyType.REFERENCE)
            self.assertEqual(adjacencies[adjacency_idx][(Haplotype.B, Haplotype.B)][0].adjacency_type,
                             AdjacencyType.REFERENCE)

    def test_get_adjacencies_from_genome_novel(self):
        genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=100, chromosomes_cnt=3, ab=True)
        adjacencies = get_adjacencies_from_genome(genome=genome, is_reference=False)
        self.assertEqual(len(adjacencies), 297)
        for adjacency_idx in adjacencies:
            self.assertEqual(len(adjacencies[adjacency_idx]), 2)
            self.assertIn((Haplotype.A, Haplotype.A), adjacencies[adjacency_idx])
            self.assertIn((Haplotype.B, Haplotype.B), adjacencies[adjacency_idx])
            self.assertEqual(len(adjacencies[adjacency_idx][(Haplotype.A, Haplotype.A)]), 1)
            self.assertEqual(len(adjacencies[adjacency_idx][(Haplotype.B, Haplotype.B)]), 1)
            self.assertEqual(adjacencies[adjacency_idx][(Haplotype.A, Haplotype.A)][0].adjacency_type,
                             AdjacencyType.NOVEL)
            self.assertEqual(adjacencies[adjacency_idx][(Haplotype.B, Haplotype.B)][0].adjacency_type,
                             AdjacencyType.NOVEL)

    def test_get_allele_specific_scn_profiles_from_genome(self):
        genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=100, chromosomes_cnt=3, ab=True)
        scn_profile_by_idx = get_scn_profile_from_genome(genome=genome)
        self.assertEqual(len(scn_profile_by_idx), 300)
        for scnr in scn_profile_by_idx.values():
            self.assertEqual(len(scnr[Haplotype.A]), 1)
            self.assertEqual(len(scnr[Haplotype.B]), 1)
            self.assertNotIn(Haplotype.UNKNOWN, scnr)


