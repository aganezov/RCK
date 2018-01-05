# -*- coding: utf-8 -*-
from collections import Counter
from copy import deepcopy

import networkx as nx

from enum import Enum

HAPLOTYPE = "haplotype"
PHASING = "phasing"


class Strand(Enum):
    REVERSE = 0
    FORWARD = 1

    def __str__(self):
        if self.value == self.REVERSE.value:
            return "-"
        return "+"

    @classmethod
    def from_pm_string(cls, string):
        if string not in ["+", "-"]:
            raise ValueError("STRAND string has to be either \"+\" or \"-\". \"{v}\" was supplied".format(v=string))
        return cls.REVERSE if string == "-" else cls.FORWARD


class PositionType(Enum):
    REAL = 0
    ARTIFICIAL = 1


class Position(object):
    STR_SEPARATOR = ":"

    def __init__(self, chromosome, coordinate, strand, ptype=PositionType.REAL, extra=None, idx=None):
        self.chromosome = chromosome
        self.coordinate = coordinate
        self.strand = strand
        self.ptype = ptype
        self.extra = extra if extra is not None else {}
        self._idx = idx

    def __eq__(self, other):
        if not isinstance(other, Position):
            return False
        return str(self) == str(other)

    @property
    def idx(self):
        if self._idx is None:
            return str(self)
        return self._idx

    def __lt__(self, other):
        if not isinstance(other, Position):
            return False
        if self.chromosome != other.chromosome:
            if len(self.chromosome) != len(other.chromosome):
                return len(self.chromosome) < len(other.chromosome)
            return self.chromosome < other.chromosome
        if self.coordinate != other.coordinate:
            return self.coordinate < other.coordinate
        return self.strand.value < other.strand.value

    def __le__(self, other):
        return self == other or self < other

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        if self.is_haplotype_specific:
            haplotype_str = ", haplotype={haplotype}{sep}".format(haplotype=str(self.extra[HAPLOTYPE]), sep=self.STR_SEPARATOR)
        else:
            haplotype_str = ""
        return "Position(chromosome={chrom}, coordinate={coord}, strand={strand}{hap_str})".format(chrom=self.chromosome,
                                                                                                   coord=self.coordinate,
                                                                                                   strand=str(self.strand),
                                                                                                   hap_str=haplotype_str)

    def __str__(self):
        if self.is_haplotype_specific:
            haplotype_str = "{haplotype}{sep}".format(haplotype=str(self.extra[HAPLOTYPE]), sep=self.STR_SEPARATOR)
        else:
            haplotype_str = ""
        return "{chrom}{sep}{haplotype_str}{strand}{coord}".format(chrom=self.chromosome,
                                                                   strand=str(self.strand),
                                                                   coord=self.coordinate,
                                                                   sep=self.STR_SEPARATOR,
                                                                   haplotype_str=haplotype_str)

    @property
    def is_haplotype_specific(self):
        return HAPLOTYPE in self.extra


class PositionCluster(object):
    def __init__(self, positions):
        self.positions = positions
        self.sort_positions()

    def sort_positions(self):
        self.positions = sorted(self.positions)

    @property
    def leftmost_position(self):
        return self.positions[0]

    @property
    def rightmost_position(self):
        return self.positions[-1]

    @property
    def length(self):
        return self.rightmost_position.coordinate - self.leftmost_position.coordinate

    @property
    def fs_cnt(self):
        return len([p for p in self.positions if p.strand == Strand.FORWARD])

    @property
    def rs_cnt(self):
        return len([p for p in self.positions if p.strand == Strand.REVERSE])


class Haplotype(Enum):
    A = "A"
    B = "B"
    UNKNOWN = "?"

    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, string):
        if string not in ["A", "B", "?"]:
            raise ValueError("Haplotype string has to be either \"A\", \"B\", or \"?\". \"{v}\" was supplied".format(v=string))
        if string == "A":
            return cls.A
        if string == "B":
            return cls.B
        if string == "?":
            return cls.UNKNOWN


class Phasing(Enum):
    AA = (Haplotype.A, Haplotype.A)
    AB = (Haplotype.A, Haplotype.B)
    BA = (Haplotype.B, Haplotype.A)
    BB = (Haplotype.B, Haplotype.B)
    AU = (Haplotype.A, Haplotype.UNKNOWN)
    BU = (Haplotype.B, Haplotype.UNKNOWN)
    UA = (Haplotype.UNKNOWN, Haplotype.A)
    UB = (Haplotype.UNKNOWN, Haplotype.B)
    UU = (Haplotype.UNKNOWN, Haplotype.UNKNOWN)

    def __str__(self):
        return "{h1}{h2}".format(h1=str(self.value[0]), h2=str(self.value[1]))


HAPLOTYPE_PAIRS_TO_PHASING = {entry.value: entry for entry in Phasing}


def haplotype_pair_to_phasing(h1, h2):
    return HAPLOTYPE_PAIRS_TO_PHASING[(h1, h2)]


def phasing_to_haplotype_pair(phasing):
    return phasing.value


class Segment(object):
    STR_REPR_SEPARATOR = ":"
    STR_COORD_SEPARATOR = "-"

    def __init__(self, start_position, end_position, extra=None, idx=None):
        if start_position.strand != Strand.REVERSE or end_position.strand != Strand.FORWARD:
            raise ValueError("START position has to come from RC strand, "
                             "while END position has to come from the F strand.")
        if start_position.coordinate >= end_position.coordinate:
            raise ValueError("START position coordinate has to be less than the END position coordinate")
        if start_position.chromosome != end_position.chromosome:
            raise ValueError("Segment's START and END position has to come from the same chromosome")
        self.start_position = start_position
        self.end_position = end_position
        self.extra = extra if extra is not None else {}
        self._idx = idx

    @property
    def idx(self):
        if self._idx is None:
            haplotype_str = "" if HAPLOTYPE not in self.extra else (str(self.extra[HAPLOTYPE]) + self.STR_REPR_SEPARATOR)
            return "{chrom}{repr_separator}{haplotype}{start}{coord_separator}{end}".format(chrom=self.chromosome,
                                                                                            start=self.start_position.coordinate,
                                                                                            end=self.end_position.coordinate,
                                                                                            haplotype=haplotype_str,
                                                                                            repr_separator=self.STR_REPR_SEPARATOR,
                                                                                            coord_separator=self.STR_COORD_SEPARATOR)
        return self._idx

    def __hash__(self):
        return hash(self.idx)

    def __eq__(self, other):
        if not isinstance(other, Segment):
            return False
        return self.idx == other.idx

    @property
    def chromosome(self):
        return self.start_position.chromosome

    @property
    def is_reversed(self):
        return self.start_position >= self.end_position

    @idx.setter
    def idx(self, value):
        self._idx = value

    def __str__(self):
        return self.idx

    def __repr__(self):
        return "Segment(id={idx}, start_position={sp}, end_position={ep}, extra={ex})" \
               "".format(sp=repr(self.start_position), ep=repr(self.end_position),
                         ex=repr(self.extra), idx=self._idx)

    @classmethod
    def from_string(cls, string):
        data = string.split(cls.STR_REPR_SEPARATOR)
        if len(data) < 2:
            raise ValueError()
        chromosome = data[0]
        if len(data) == 2:  # withOUT haplotype info
            start_coord, end_coord = data[1].split(cls.STR_COORD_SEPARATOR)
            start_coord, end_coord = int(start_coord), int(end_coord)
            extra = None
        else:  # with haplotype info
            start_coord, end_coord = data[2].split(cls.STR_COORD_SEPARATOR)
            start_coord, end_coord = int(start_coord), int(end_coord)
            haplotype = data[1]
            haplotype = Haplotype.from_string(string=haplotype)
            extra = {HAPLOTYPE: haplotype}
        segment_is_reversed = start_coord > end_coord
        if segment_is_reversed:
            start_position = Position(chromosome=chromosome, coordinate=end_coord, strand=Strand.REVERSE)
            end_position = Position(chromosome=chromosome, coordinate=start_coord, strand=Strand.FORWARD)
        else:
            start_position = Position(chromosome=chromosome, coordinate=start_coord, strand=Strand.REVERSE)
            end_position = Position(chromosome=chromosome, coordinate=end_coord, strand=Strand.FORWARD)
        segment = cls(start_position=start_position, end_position=end_position, extra=extra)
        if segment_is_reversed:
            segment.start_position, segment.end_position = segment.end_position, segment.start_position
        return segment


class PhylogenyTree(object):
    def __init__(self, tree, tree_root):
        self.root = "R"
        self.tree = tree
        self.mutation_tree_root = tree_root
        self.tree.add_edge((self.root, self.mutation_tree_root))


class AdjacencyType(Enum):
    REFERENCE = 0
    NOVEL = 1


class Adjacency(object):
    def __init__(self, position1, position2, adjacency_type=None, idx=None, extra=None):
        if adjacency_type is None:
            adjacency_type = AdjacencyType.NOVEL
        if position1 <= position2:
            self.position1 = position1
            self.position2 = position2
        else:
            self.position1 = position2
            self.position2 = position1
        self._idx = idx
        self.extra = extra if extra is not None else {}
        self.adjacency_type = adjacency_type

    @property
    def idx(self):
        if self._idx is None:
            return str(self)
        return self._idx

    @idx.setter
    def idx(self, value):
        self._idx = value

    def __eq__(self, other):
        if not isinstance(other, Adjacency):
            return False
        return self.idx == other.idx

    def __hash__(self):
        return hash(self.idx)

    def __repr__(self):
        return "Adjacency(id={idx}, position1={p1}, position2={p2})" \
               "".format(p1=repr(self.position1), p2=repr(self.position2), idx=self._idx)

    def __str__(self):
        return "[{p1}]-[{p2}]".format(p1=str(self.position1), p2=str(self.position2))


class SegmentCopyNumber(object):
    def __init__(self, segment, cn):
        self.segment = segment
        self.cn = cn

    def __eq__(self, other):
        if not isinstance(other, SegmentCopyNumber):
            return False
        return self.segment == other.segment and self.cn == other.cn


class AdjacencyCopyNumber(object):
    def __init__(self, adjacency, cn):
        self.adjacency = adjacency
        self.cn = cn


SegmentCN = SegmentCopyNumber
AdjacencyCN = AdjacencyCopyNumber


class SegmentCopyNumberRecord(object):
    def __init__(self, segment, maj_a_cn, maj_b_cn, min_a_cn=None, min_b_cn=None, allele_specific=True, extra=None):
        self.segment = segment
        self.allele_specific = allele_specific
        self.maj_a_segmentcn = SegmentCopyNumber(segment=segment, cn=maj_a_cn)
        self.maj_b_segmentcn = SegmentCopyNumber(segment=segment, cn=maj_b_cn)
        if min_a_cn is None:
            self.min_a_segmentcn = self.maj_a_segmentcn
        else:
            self.min_a_segmentcn = SegmentCopyNumber(segment=segment, cn=min_a_cn)
        if min_b_cn is None:
            self.min_b_segmentcn = self.maj_b_segmentcn
        else:
            self.min_b_segmentcn = SegmentCopyNumber(segment=segment, cn=min_b_cn)
        if not self.allele_specific and self.maj_b_segmentcn.cn != 0:
            raise ValueError("Segment copy number record is specified as 'non-allele-specific', but value B allele ({b_cn}) is non zero.".format(b_cn=self.maj_b_segmentcn.cn))
        self.extra = extra if extra is not None else {}

    @property
    def distinct_min_cn(self):
        return self.maj_a_segmentcn != self.min_a_segmentcn or self.maj_b_segmentcn != self.min_b_segmentcn

    @classmethod
    def to_non_allele_specific(cls, scnr):
        return cls(segment=scnr.segment,
                   maj_a_cn=scnr.maj_a_segmentcn.cn + scnr.maj_b_segmentcn.cn,
                   maj_b_cn=0,
                   min_a_cn=scnr.min_a_segmentcn.cn + scnr.min_b_segmentcn.cn,
                   min_b_cn=0,
                   allele_specific=False,
                   extra={"source_scnr": scnr})

    def __eq__(self, other):
        if not isinstance(other, SegmentCopyNumberRecord):
            return False
        return self.maj_a_segmentcn == other.maj_a_segmentcn and \
               self.maj_b_segmentcn == other.maj_b_segmentcn and \
               self.min_a_segmentcn == other.min_a_segmentcn and \
               self.min_b_segmentcn == other.min_b_segmentcn and \
               self.allele_specific == other.allele_specific


SegmentCNRecord = SegmentCopyNumberRecord
SCNR = SegmentCNRecord


class AdjacencyGroup(object):
    def __init__(self, adjacencies, idx, fp=0.0):
        self.adjacencies = adjacencies
        self.fp = fp
        self.idx = idx


def get_segments_from_genome(genome, copy=True, make_all_non_reversed=True):
    result = []
    for chromosome in genome:
        for s in chromosome:
            v = s
            if copy:
                v = deepcopy(s)
            if v.is_reversed and make_all_non_reversed:
                reverse_segment(segment=v, copy=False)
            result.append(v)
    return result


def get_telomeres_from_genome(genome, copy=True):
    result = []
    for chromosome in genome:
        ss, es = chromosome[0], chromosome[-1]
        lt = ss.start_position
        if copy:
            lt = deepcopy(ss.start_position)
        rt = es.end_position
        if copy:
            rt = deepcopy(es.end_position)
        if HAPLOTYPE in ss.extra:
            lt.extra[HAPLOTYPE] = ss.extra[HAPLOTYPE]
        if HAPLOTYPE in es.extra:
            rt.extra[HAPLOTYPE] = es.extra[HAPLOTYPE]
        result.append(lt)
        result.append(rt)
    return result


def strip_phasing_from_adjacencies(adjacencies, inplace=True, strip_positions_haplotypes=True):
    result = []
    for a in adjacencies:
        v = a
        if not inplace:
            v = deepcopy(a)
        if HAPLOTYPE in v.position1.extra and strip_positions_haplotypes:
            del v.position1.extra[HAPLOTYPE]
        if HAPLOTYPE in v.position2.extra and strip_positions_haplotypes:
            del v.position2.extra[HAPLOTYPE]
        if PHASING in v.extra:
            del v.extra[PHASING]
        result.append(v)
    return result


def strip_haplotype_from_segments(segments, inplace=True, strip_positions_haplotypes=True):
    result = []
    for s in segments:
        v = s
        if not inplace:
            v = deepcopy(s)
        if HAPLOTYPE in v.extra:
            del v.extra[HAPLOTYPE]
        if strip_positions_haplotypes and HAPLOTYPE in v.start_position.extra:
            del v.start_position.extra[HAPLOTYPE]
        if strip_positions_haplotypes and HAPLOTYPE in v.end_position.extra:
            del v.end_position.extra[HAPLOTYPE]
        result.append(v)
    return result


def strip_haplotype_from_positions(positions, inplace=True):
    result = []
    for p in positions:
        v = p
        if not inplace:
            v = deepcopy(p)
        if HAPLOTYPE in v.extra:
            del v.extra[HAPLOTYPE]
        result.append(v)
    return result


def get_adjacencies_from_genome(genome, copy=True, inherit_haplotypes=True, default_adjacency_type=AdjacencyType.NOVEL):
    result = []
    for chromosome in genome:
        for s1, s2 in zip(chromosome[:-1], chromosome[1:]):
            p1 = s1.end_position
            if copy:
                p1 = deepcopy(s1.end_position)
            p2 = s2.start_position
            if copy:
                p2 = deepcopy(s2.start_position)
            if inherit_haplotypes:
                if HAPLOTYPE not in p1.extra and HAPLOTYPE in s1.extra:
                    p1.extra[HAPLOTYPE] = s1.extra[HAPLOTYPE]
                if HAPLOTYPE not in p2.extra and HAPLOTYPE in s2.extra:
                    p2.extra[HAPLOTYPE] = s2.extra[HAPLOTYPE]
            adjacency = Adjacency(position1=p1, position2=p2, adjacency_type=default_adjacency_type)
            result.append(adjacency)
    return result


def get_unique_adjacencies(adjacencies, copy=True):
    unique = set()
    for a in adjacencies:
        v = a
        if copy:
            v = deepcopy(a)
        unique.add(v)
    return list(unique)


def assign_ref_adjacency_status(adjacencies, ref_adjacencies=None, inplace=True):
    result = []
    for a in adjacencies:
        v = a
        if not inplace:
            v = deepcopy(a)
        if ref_adjacencies is None or v in ref_adjacencies:
            v.adjacency_type = AdjacencyType.REFERENCE
        result.append(v)
    return result


def assign_adjacency_status(adjacencies, ref_adjacencies=None, inplace=True):
    result = []
    for a in adjacencies:
        v = a
        if not inplace:
            v = deepcopy(a)
        if ref_adjacencies is not None:
            if v in ref_adjacencies:
                v.adjacency_type = AdjacencyType.REFERENCE
            else:
                v.adjacency_type = AdjacencyType.NOVEL
        result.append(v)
    return result


def assign_nov_adjacency_status(adjacencies, ref_adjacencies=None, inplace=True):
    result = []
    for a in adjacencies:
        v = a
        if not inplace:
            v = deepcopy(a)
        if ref_adjacencies is None or v not in ref_adjacencies:
            v.adjacency_type = AdjacencyType.NOVEL
        result.append(v)
    return result


def reverse_segment(segment, copy=True):
    result = segment
    if copy:
        result = deepcopy(segment)
    result.start_position, result.end_position = result.end_position, result.start_position
    return result


def propagate_haplotype_segment_to_positions(segment, inplace=True):
    result = segment
    if not inplace:
        result = deepcopy(segment)
    if HAPLOTYPE in result.extra:
        result.start_position.extra[HAPLOTYPE] = result.extra[HAPLOTYPE]
        result.end_position.extra[HAPLOTYPE] = result.extra[HAPLOTYPE]
    return result


def propagate_haplotype_chrom_to_positions(chromosome, inplace=True):
    result = []
    for s in chromosome:
        v = propagate_haplotype_segment_to_positions(segment=s, inplace=inplace)
        result.append(v)
    return result


def propagate_haplotype_genome_to_positions(genome, inplace=True):
    result = []
    for chromosome in genome:
        v = propagate_haplotype_chrom_to_positions(chromosome=chromosome, inplace=inplace)
        result.append(v)
    return result


def get_segments_copy_number_profile(segments):
    return Counter(segments)


def get_adjacencies_copy_number_profile(adjacencies):
    return Counter(adjacencies)
