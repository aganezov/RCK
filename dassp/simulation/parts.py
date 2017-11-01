from copy import deepcopy

from dassp.core.structures import Haplotype, Position, Strand, PositionType, Segment

CHROMOSOME_SIZE = 100
CHROMOSOMES_CNT = 3


class ChromosomeGenerator(object):
    @classmethod
    def generate_chromosomes(cls,
                             chromosome_size=CHROMOSOME_SIZE,
                             chromosomes_cnt=CHROMOSOMES_CNT,
                             ab=True):
        result = []
        for chr_cnt in range(chromosomes_cnt):
            a_chromosome = cls.generate_chromosome(chromosome_size=chromosome_size, chr_name=str(chr_cnt))
            result.append(a_chromosome)
            if ab:
                b_chromosome = cls.get_b_chromosome_from_a(a_chromosome=a_chromosome)
                result.append(b_chromosome)
        return result

    @classmethod
    def generate_chromosome(cls, chromosome_size, chr_name, haplotype=Haplotype.A):
        result = []
        for coordinate in range(chromosome_size):
            star_position = Position(chromosome=chr_name, coordinate=coordinate * 2, strand=Strand.REVERSE,
                                     ptype=PositionType.ARTIFICIAL)
            end_position = Position(chromosome=chr_name, coordinate=coordinate * 2 + 1, strand=Strand.FORWARD,
                                    ptype=PositionType.ARTIFICIAL)
            segment = Segment(start_position=star_position, end_position=end_position, extra={"haplotype": haplotype})
            result.append(segment)
        return result

    @classmethod
    def get_b_chromosome_from_a(cls, a_chromosome):
        result = deepcopy(a_chromosome)
        for s in result:
            s.extra["haplotype"] = Haplotype.B
        return result


def reverse_segment(segment):
    segment.start_position, segment.end_position = segment.end_position, segment.start_position
    return segment


def delete_segments(chromosome, start_segment_index, end_segment_index):
    if start_segment_index >= len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "deletion is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in deletion is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    return deepcopy(chromosome[:start_segment_index]) + deepcopy(chromosome[end_segment_index:])


def reverse_segments(chromosome, start_segment_index, end_segment_index):
    if len(chromosome) == 0:
        return []
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in reversal is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    if start_segment_index >= len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "reversal is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    span_to_reverse = chromosome[start_segment_index:end_segment_index]
    reversed_span = []
    for s in span_to_reverse[::-1]:
        reversed_span.append(reverse_segment(deepcopy(s)))
    return deepcopy(chromosome[:start_segment_index]) + reversed_span + deepcopy(chromosome[end_segment_index:])


def duplicate_segments(chromosome, start_segment_index, end_segment_index):
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in duplication is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    if start_segment_index >= len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "duplication is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    duplicated_span = deepcopy(chromosome[start_segment_index:end_segment_index])
    return deepcopy(chromosome[:end_segment_index]) + duplicated_span + deepcopy(chromosome[end_segment_index:])


def translocation_segments(chromosome1,
                           chromosome2,
                           chromosome1_segment_index,
                           chromosome2_segment_index,
                           cc=True):
    if chromosome1_segment_index > len(chromosome1):
        raise ValueError("Translocation index {ti} is greater than the length of the first chromosome {chr_l}"
                         "".format(ti=chromosome1_segment_index, chr_l=len(chromosome1)))
    if chromosome2_segment_index > len(chromosome2):
        raise ValueError("Translocation index {ti} is greater than the length of the second chromosome {chr_l}"
                         "".format(ti=chromosome2_segment_index, chr_l=len(chromosome2)))
    chr1_beg, chr1_end = deepcopy(chromosome1[:chromosome1_segment_index]), deepcopy(chromosome1[chromosome1_segment_index:])
    chr2_beg, chr2_end = deepcopy(chromosome2[:chromosome2_segment_index]), deepcopy(chromosome2[chromosome2_segment_index:])
    if cc:
        new_chr1 = chr1_beg + chr2_end
        new_chr2 = chr2_beg + chr1_end
    else:
        new_chr1 = chr1_beg + reverse_segments(chromosome=chr2_beg, start_segment_index=0, end_segment_index=len(chr2_beg))
        new_chr2 = reverse_segments(chromosome=chr1_end, start_segment_index=0, end_segment_index=len(chr1_end)) + chr2_end
    return new_chr1, new_chr2
