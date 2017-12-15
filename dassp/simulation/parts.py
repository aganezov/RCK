from collections import defaultdict, Counter
from copy import deepcopy
from enum import Enum

from dassp.core.structures import Haplotype, Position, Strand, PositionType, Segment

CHROMOSOME_SIZE = 100
CHROMOSOMES_CNT = 3

HAPLOTYPE = "haplotype"


class ChromosomeGenerator(object):
    @classmethod
    def generate_chromosomes(cls,
                             chromosome_size=CHROMOSOME_SIZE,
                             chromosomes_cnt=CHROMOSOMES_CNT,
                             ab=True):
        result = []
        for chr_cnt in range(1, chromosomes_cnt + 1):
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
            segment = Segment(start_position=star_position, end_position=end_position, extra={HAPLOTYPE: haplotype})
            result.append(segment)
        return result

    @classmethod
    def get_b_chromosome_from_a(cls, a_chromosome):
        result = deepcopy(a_chromosome)
        for s in result:
            s.extra[HAPLOTYPE] = Haplotype.B
        return result


def reverse_segment(segment, copy=True):
    if not copy:
        segment.start_position, segment.end_position = segment.end_position, segment.start_position
        return segment
    result = deepcopy(segment)
    result.start_position, result.end_position = result.end_position, result.start_position
    return result


def delete_segments(chromosome, start_segment_index, end_segment_index, copy=True):
    if start_segment_index >= len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "deletion is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in deletion is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    if copy:
        return deepcopy(chromosome[:start_segment_index]) + deepcopy(chromosome[end_segment_index:])
    return chromosome[:start_segment_index] + chromosome[:end_segment_index]


def reverse_segments(chromosome, start_segment_index, end_segment_index, copy=True):
    if len(chromosome) == 0:
        return []
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in reversal is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    if start_segment_index > len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "reversal is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    span_to_reverse = chromosome[start_segment_index:end_segment_index]
    reversed_span = []
    for s in span_to_reverse[::-1]:
        reversed_span.append(reverse_segment(s, copy=True))
    if copy:
        deepcopy(chromosome[:start_segment_index]) + reversed_span + deepcopy(chromosome[end_segment_index:])
    return chromosome[:start_segment_index] + reversed_span + chromosome[:end_segment_index:]


def duplicate_segments(chromosome,
                       start_segment_index,  # beginning of a span to be duplicated
                       end_segment_index,  # end of a span to be duplicated
                       insert_shift=0,  # a shift from the beginning (if <0) of the end (if >=0) to the location where the inserted span to be placed
                       reverse_copy=False,  # whether to insert an original / reversed copy of a span of segments
                       copy=True):
    if start_segment_index > end_segment_index >= 0:
        raise ValueError("Index {li} of the last segment in duplication is smaller than the start index {si}"
                         "".format(li=end_segment_index, si=start_segment_index))
    if start_segment_index >= len(chromosome):
        raise ValueError("An index {si} for the segment corresponding to the beginning of the "
                         "duplication is larger or equal to than the total number of segments {cnt} in the chromosomes"
                         "".format(si=start_segment_index, cnt=len(chromosome)))
    if insert_shift < 0:
        insert_segment_index = start_segment_index + insert_shift
    else:
        insert_segment_index = end_segment_index + insert_shift
    if insert_segment_index < 0:
        raise ValueError("A shift of ({i_s}) for the inserted sequence of segments attempts to places the copy outside the chromosome ({s_i}{i_s}={i_si})"
                         "".format(i_s=insert_shift, s_i=start_segment_index, i_si=insert_segment_index))
    elif insert_segment_index > len(chromosome):
        raise ValueError("A shift of ({i_s}) for the inserted sequence of segments attempts to places the copy outside the chromosome ({e_i}+{i_s}={i_si})"
                         "".format(i_s=insert_shift, e_i=end_segment_index, i_si=insert_segment_index))
    duplicated_span = deepcopy(chromosome[start_segment_index:end_segment_index])
    if reverse_copy:
        duplicated_span = reverse_segments(chromosome=duplicated_span, start_segment_index=0, end_segment_index=len(duplicated_span), copy=False)
    if copy:
        deepcopy(chromosome[:insert_segment_index]) + duplicated_span + deepcopy(chromosome[insert_segment_index:])
    return chromosome[:insert_segment_index] + duplicated_span + chromosome[insert_segment_index:]


def translocate_segments(chromosome1,
                         chromosome2,
                         chromosome1_segment_index,
                         chromosome2_segment_index,
                         cc=True,
                         copy=True):
    if chromosome1_segment_index > len(chromosome1):
        raise ValueError("Translocation index {ti} is greater than the length of the first chromosome {chr_l}"
                         "".format(ti=chromosome1_segment_index, chr_l=len(chromosome1)))
    if chromosome2_segment_index > len(chromosome2):
        raise ValueError("Translocation index {ti} is greater than the length of the second chromosome {chr_l}"
                         "".format(ti=chromosome2_segment_index, chr_l=len(chromosome2)))

    chr1_beg, chr1_end = chromosome1[:chromosome1_segment_index], chromosome1[chromosome1_segment_index:]
    chr2_beg, chr2_end = chromosome2[:chromosome2_segment_index], chromosome2[chromosome2_segment_index:]
    if copy:
        chr1_beg, chr1_end, chr2_beg, chr2_end = deepcopy(chr1_beg), deepcopy(chr1_end), deepcopy(chr2_beg), deepcopy(chr2_end)
    if cc:
        new_chr1 = chr1_beg + chr2_end
        new_chr2 = chr2_beg + chr1_end
    else:
        new_chr1 = chr1_beg + reverse_segments(chromosome=chr2_beg, start_segment_index=0, end_segment_index=len(chr2_beg), copy=False)
        new_chr2 = reverse_segments(chromosome=chr1_end, start_segment_index=0, end_segment_index=len(chr1_end), copy=False) + chr2_end
    return new_chr1, new_chr2


def split_and_reassemble_chromosomes(chromosomes,
                                     per_chrs_segment_indexes,
                                     target_structures,
                                     ensure_cut_indexes_valid=True,
                                     ensure_target_structure_valid=True,
                                     ensure_no_deletions=False,
                                     ensure_no_duplications=False,
                                     ensure_telomere_preservation=False,
                                     ensre_chromosomes_cnt_preservation=False,
                                     copy=True):
    """

    :param chromosomes: a list of chromosomes (i.e., lists of segments)
    :param per_chrs_segment_indexes: a list of lists, where each internal list contains a indexes of where respective chromosomes hav to be "chopped"
    :param target_structures: a list of new chromosome configurations, each configuration is represented as a list of entries, where each entry is a triplet determining
        a "span" of segments obtained from the original chromosomal chopping as well as whether to use its original version or the reversed one
    :param ensure_cut_indexes_valid:
    :param ensure_target_structure_valid:
    :param ensure_no_deletions:
    :param ensure_no_duplications:
    :param ensure_telomere_preservation:
    :param ensre_chromosomes_cnt_preservation:
    :param copy:
    :return:
    """
    if ensure_cut_indexes_valid:
        if not chromothripsis_cut_indexes_are_valid(chromosomes=chromosomes,
                                                    per_chrs_segment_indexes=per_chrs_segment_indexes,
                                                    target_structure=target_structures):
            raise ValueError()
    if ensure_target_structure_valid:
        if not chromothripsis_target_structure_is_valid(chromosomes=chromosomes,
                                                        per_chrs_segment_indexes=per_chrs_segment_indexes,
                                                        target_structure=target_structures):
            raise ValueError()
    if ensure_no_deletions:
        if not chromothripsis_no_deletions(chromosomes=chromosomes,
                                           per_chrs_segment_indexes=per_chrs_segment_indexes,
                                           target_structure=target_structures):
            raise ValueError()
    if ensure_no_duplications:
        if not chromothripsis_no_duplications(chromosomes=chromosomes,
                                              per_chrs_segment_indexes=per_chrs_segment_indexes,
                                              target_structure=target_structures):
            raise ValueError()
    if ensure_telomere_preservation:
        if not chromothripsis_preserves_telomeres(chromosomes=chromosomes,
                                                  per_chrs_segments_indexes=per_chrs_segment_indexes,
                                                  target_structure=target_structures):
            raise ValueError()
    if ensre_chromosomes_cnt_preservation:
        if not chromothripsis_preserves_number_of_chromosomes(chromosomes=chromosomes,
                                                              per_chrs_segment_indexes=per_chrs_segment_indexes,
                                                              target_structure=target_structures):
            raise ValueError()
    parts = defaultdict(dict)                       # shallow copies of parts of the original chromosomes
    ###
    # chop chromosomes into parts
    ###
    for cnt, chromosome in enumerate(chromosomes):
        current_part_number = 0
        current_l_segment_index = 0
        for current_r_segment_index in per_chrs_segment_indexes[0]:
            parts[cnt][current_part_number] = chromosome[current_l_segment_index:current_r_segment_index]
            current_l_segment_index = current_r_segment_index
            current_part_number += 1
        parts[cnt][current_part_number] = chromosome[current_l_segment_index:len(chromosome)]
    ###
    # reassemble chromosomes from the parts as specified by the target_structure info
    ###
    result = []
    for structure in target_structures:
        result_chromosome = []
        for origin_chromosome_i, span_index, reverse_use in structure:
            span = parts[origin_chromosome_i][span_index]
            if reverse_use:
                span = reverse_segments(chromosome=span, start_segment_index=0, end_segment_index=len(span), copy=copy)
            if copy:
                result_chromosome.extend(deepcopy(span))
            else:
                result_chromosome.extend(span)
        result.append(result_chromosome)
    return result


def chromothripsis_cut_indexes_are_valid(chromosomes, per_chrs_segment_indexes, target_structure):
    """
    number of index groups has to be equal to the number of chromosomes
    all indexes within a group for each chromosome have to be distinct
    for each chromosome indexes have to be compliant with chromosomal length
    """
    if len(per_chrs_segment_indexes) > len(chromosomes):
        return False
    for chromosome, per_chr_indexes in zip(chromosomes, per_chrs_segment_indexes):
        if len(per_chr_indexes) == 0:
            return False
        if len(set(per_chr_indexes)) != len(per_chr_indexes):
            return False
        for index in per_chr_indexes:
            if index < 0 or index > len(chromosome):
                return False
    return True


def chromothripsis_target_structure_is_valid(chromosomes, per_chrs_segment_indexes, target_structure):
    """
    all entries across each chromosome group in the target structure have to:
        * have a chromosome index that references one of the input chromosomes
        * reference a span index within a referenced chromosome such that it complies with how the referenced chromosome was "chopped"
        * flag of whether to use the original span or its reversed version has to be a boolean
    """
    possible_parts = {}
    for chr_cnt, chromosome in enumerate(chromosomes):
        possible_parts[chr_cnt] = len(per_chrs_segment_indexes[chr_cnt]) + 1
    for chr_structure in target_structure:
        if len(chr_structure) == 0:
            return False
        for entry in chr_structure:
            chromosome_index, span_index, reverse_use = entry
            if chromosome_index not in possible_parts:
                return False
            if span_index < 0 or span_index >= possible_parts[chromosome_index]:
                return False
            if not isinstance(reverse_use, bool):
                return False
    return True


def chromothripsis_no_deletions(chromosomes, per_chrs_segment_indexes, target_structure):
    to_be_used = set()
    for cnt, chromosome in enumerate(chromosomes):
        for i in range(per_chrs_segment_indexes[cnt] + 1):
            to_be_used.add((cnt, i))
    used = set()
    for structure in target_structure:
        for entry in structure:
            chromosome_index, span_index, _ = entry
            used.add((chromosome_index, span_index))
    return len(to_be_used - used) == 0


def chromothripsis_no_duplications(chromosomes, per_chrs_segment_indexes, target_structure):
    used = []
    for structure in target_structure:
        for entry in structure:
            chromosome_index, span_index, _ = entry
            used.append((chromosome_index, span_index))
    return len(set(used)) == len(used)


def chromothripsis_is_cnn(chromosomes, per_chrs_segment_indexes, target_structure):
    no_deletions = chromothripsis_no_deletions(chromosomes=chromosomes, per_chrs_segment_indexes=per_chrs_segment_indexes, target_structure=target_structure)
    no_duplications = chromothripsis_no_duplications(chromosomes=chromosomes, per_chrs_segment_indexes=per_chrs_segment_indexes, target_structure=target_structure)
    return no_deletions and no_duplications


def chromothripsis_preserves_telomeres(chromosomes, per_chrs_segments_indexes, target_structure):
    original_left_telomeres = set()
    original_right_telomeres = set()
    for chr_cnt, chromosome in enumerate(chromosomes):
        original_left_telomeres.add((chr_cnt, 0))
        original_right_telomeres.add((chr_cnt, len(per_chrs_segments_indexes[chr_cnt])))
    used_telomeres = set()
    for structure in target_structure:
        ###
        # checking novel left telomere
        ###
        chr_index, span_index, reverse_use = structure[0]
        if reverse_use:
            if (chr_index, span_index) not in original_right_telomeres:
                return False
            used_telomeres.add((chr_index, span_index))
        else:
            if (chr_index, span_index) not in original_left_telomeres:
                return False
            used_telomeres.add((chr_index, span_index))
        ###
        # checking novel right telomere
        ###
        chr_index, span_index, reverse_use = structure[-1]
        if reverse_use:
            if (chr_index, span_index) not in original_left_telomeres:
                return False
            used_telomeres.add((chr_index, span_index))
        else:
            if (chr_index, span_index) not in original_right_telomeres:
                return False
            used_telomeres.add((chr_index, span_index))
    return len(original_right_telomeres + original_left_telomeres - used_telomeres) == 0


def chromothripsis_preserves_number_of_chromosomes(chromosomes, per_chrs_segment_indexes, target_structure):
    return len(target_structure) == len(chromosomes)


class MutationType(Enum):
    DELETION = 0
    DUPLICATION = 1
    REVERSAL = 2
    TRANSLOCATION = 3
    CHROMOTHRIPSIS = 6

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class Mutation(object):
    def __init__(self, mutation_type, mutation_data):
        self.mutation_type = mutation_type
        self.mutation_data = mutation_data
