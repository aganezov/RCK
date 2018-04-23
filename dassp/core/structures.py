# -*- coding: utf-8 -*-
import numpy as np
import math
from collections import Counter, namedtuple, defaultdict
from copy import deepcopy
import networkx as nx
import itertools

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

    def __lt__(self, other):
        if not isinstance(other, Haplotype):
            return False
        return self.value < other.value


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
        return hash(self) == hash(other)

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
        if self.is_haplotype_specific and other.is_haplotype_specific and self.extra[HAPLOTYPE] != other.extra[HAPLOTYPE]:
            return self.extra[HAPLOTYPE] < other.extra[HAPLOTYPE]
        if self.coordinate != other.coordinate:
            return self.coordinate < other.coordinate
        return self.strand.value < other.strand.value

    def less_than__non_hap(self, other):
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
        if self.is_haplotype_specific:
            return hash(self.stable_id_hap)
        return hash(self.stable_id_non_hap)

    def __repr__(self):
        if self.is_haplotype_specific:
            haplotype_str = ", haplotype={haplotype}".format(haplotype=str(self.extra[HAPLOTYPE]), sep=self.STR_SEPARATOR)
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

    def get_haplotype(self, default=Haplotype.A):
        if self.is_haplotype_specific:
            return self.haplotype
        return default

    @property
    def haplotype(self):
        if not self.is_haplotype_specific:
            raise ValueError()
        return self.extra[HAPLOTYPE]

    @property
    def stable_id_non_hap(self):
        return "{chrom}{sep}{strand}{coord}".format(chrom=self.chromosome,
                                                    strand=str(self.strand),
                                                    coord=self.coordinate,
                                                    sep=self.STR_SEPARATOR)

    @property
    def stable_id_hap(self):
        if not self.is_haplotype_specific:
            raise ValueError()
        return self.id_hap_string_from_elements(chrom=self.chromosome,
                                                strand=str(self.strand),
                                                coord=self.coordinate,
                                                sep=self.STR_SEPARATOR,
                                                hap=str(self.haplotype))

    @classmethod
    def id_hap_string_from_elements(cls, chrom, hap, strand, coord, sep=None):
        if sep is None:
            sep = cls.STR_SEPARATOR
        return "{chrom}{sep}{haplotype}{sep}{strand}{coord}".format(chrom=str(chrom),
                                                                    strand=str(strand),
                                                                    coord=str(coord),
                                                                    sep=sep,
                                                                    haplotype=str(hap))

    def get_non_hap_copy(self):
        result = deepcopy(self)
        if result.is_haplotype_specific:
            del result.extra[HAPLOTYPE]
        return result


class PositionCluster(object):
    def __init__(self, positions):
        self.positions = positions
        self.sort_positions()

    def sort_positions(self, coordinate_only=True):
        if coordinate_only:
            key = lambda p: p.coordinate
        else:
            key = lambda p: p
        self.positions = sorted(self.positions, key=key)

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

    @property
    def cluster_mean_coordinate(self):
        return int(np.mean([p.coordinate for p in self.positions]))

    @property
    def cluster_pos_coordinate(self):
        return []


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

    @staticmethod
    def flip(phasing):
        h1, h2 = phasing.value
        return haplotype_pair_to_phasing(h1=h2, h2=h1)


HAPLOTYPE_PAIRS_TO_PHASING = {entry.value: entry for entry in Phasing}


def haplotype_pair_to_phasing(h1, h2):
    return HAPLOTYPE_PAIRS_TO_PHASING[(h1, h2)]


def phasing_to_haplotype_pair(phasing):
    return phasing.value


def flipped_phasing(phasing):
    return Phasing.flip(phasing=phasing)


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
        if start_position.is_haplotype_specific and end_position.is_haplotype_specific and start_position.haplotype != end_position.haplotype:
            raise ValueError("Haplotypes for segment do not match")
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
        if self.is_haplotype_specific:
            return hash(self.stable_id_hap)
        return hash(self.stable_id_non_hap)

    def __eq__(self, other):
        if not isinstance(other, Segment):
            return False
        return hash(self) == hash(other)

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

    @property
    def is_haplotype_specific(self):
        if HAPLOTYPE in self.extra:
            return True
        if self.start_position.is_haplotype_specific and \
                self.end_position.is_haplotype_specific and \
                self.start_position.haplotype == self.end_position.haplotype:
            return True
        return False

    @property
    def haplotype(self):
        if not self.is_haplotype_specific:
            raise ValueError()
        if HAPLOTYPE in self.extra:
            return self.extra[HAPLOTYPE]
        return self.start_position.haplotype

    def get_haplotype(self, default=Haplotype.A):
        if self.is_haplotype_specific:
            return self.haplotype
        return default

    @property
    def stable_id_non_hap(self):
        sp, ep = self.start_position, self.end_position
        if ep.less_than__non_hap(sp):
            sp, ep = ep, sp
        return "{chrom}{repr_separator}{start}{coord_separator}{end}".format(chrom=self.chromosome,
                                                                             start=sp.coordinate,
                                                                             end=ep.coordinate,
                                                                             repr_separator=self.STR_REPR_SEPARATOR,
                                                                             coord_separator=self.STR_COORD_SEPARATOR)

    @property
    def stable_id_hap(self):
        if not self.is_haplotype_specific:
            raise ValueError()
        sp, ep = tuple(sorted([self.start_position, self.end_position]))
        hap = self.haplotype
        return "{chrom}{repr_separator}{haplotype}{repr_separator}{start}{coord_separator}{end}".format(chrom=self.chromosome,
                                                                                                        start=sp.coordinate,
                                                                                                        end=ep.coordinate,
                                                                                                        haplotype=str(hap),
                                                                                                        repr_separator=self.STR_REPR_SEPARATOR,
                                                                                                        coord_separator=self.STR_COORD_SEPARATOR)

    def get_non_hap_copy(self):
        is_reversed = self.is_reversed
        result = self.__class__(start_position=self.start_position.get_non_hap_copy(),
                                end_position=self.end_position.get_non_hap_copy())
        if is_reversed:
            reverse_segment(segment=result, copy=False)
        return result

    @property
    def length(self):
        return int(self.end_position.coordinate - self.start_position.coordinate)

    @property
    def length_1000(self):
        return int(math.ceil(int(self.end_position.coordinate - self.start_position.coordinate) / 1000.0))


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
    def __init__(self, position1, position2, adjacency_type=AdjacencyType.NOVEL, idx=None, extra=None):

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

    @property
    def is_sorted_non_phased(self):
        return not self.position2.less_than__non_hap(other=self.position1)

    @property
    def is_sorted_phased(self):
        if not self.is_phased:
            raise ValueError()
        if self.position1.chromosome != self.position2.chromosome:
            if len(self.position1.chromosome) != len(self.position2.chromosome):
                return len(self.position1.chromosome) < len(self.position2.chromosome)
            return self.position1.chromosome < self.position2.chromosome
        phasing = self.phasing
        if phasing.value[0] != phasing.value[1]:
            return phasing.value[0] < phasing.value[1]
        if self.position1.coordinate != self.position2.coordinate:
            return self.position1.coordinate < self.position2.coordinate
        return self.position1.strand.value < self.position2.strand.value

    @property
    def is_phased(self):
        return PHASING in self.extra or (self.position1.is_haplotype_specific and self.position2.is_haplotype_specific)

    @property
    def phasing(self):
        if not self.is_phased:
            raise ValueError()
        if PHASING in self.extra:
            return self.extra[PHASING]
        return haplotype_pair_to_phasing(h1=self.position1.haplotype, h2=self.position2.haplotype)

    @property
    def stable_phasing(self):
        phasing = self.phasing
        if self.position2.less_than__non_hap(self.position1):
            return flipped_phasing(phasing=phasing)
        return phasing

    def get_phasing(self, sort=True, default=Phasing.AA):
        if not self.is_phased:
            return default
        if sort:
            return self.stable_phasing
        return self.phasing

    @property
    def stable_id_non_phased(self):
        p1_id, p2_id = self.position1.stable_id_non_hap, self.position2.stable_id_non_hap
        if self.position2.less_than__non_hap(self.position1):
            p1_id, p2_id = p2_id, p1_id
        return "[{p1}]-[{p2}]".format(p1=p1_id, p2=p2_id)

    @property
    def id_non_phased(self):
        p1_id, p2_id = self.position1.stable_id_non_hap, self.position2.stable_id_non_hap
        return "[{p1}]-[{p2}]".format(p1=p1_id, p2=p2_id)

    @property
    def stable_id_phased(self):
        if not self.is_phased:
            raise ValueError()
        phasing = self.phasing
        p1, p2 = self.position1, self.position2
        if p1 > p2:
            phasing = flipped_phasing(phasing=phasing)
            p1, p2 = p2, p1
        p1_id = Position.id_hap_string_from_elements(chrom=p1.chromosome, hap=phasing.value[0], strand=p1.strand, coord=p1.coordinate)
        p2_id = Position.id_hap_string_from_elements(chrom=p2.chromosome, hap=phasing.value[1], strand=p2.strand, coord=p2.coordinate)
        return "[{p1}]-[{p2}]".format(p1=p1_id, p2=p2_id)

    @property
    def id_phased(self):
        if not self.is_phased:
            raise ValueError()
        phasing = self.phasing
        p1, p2 = self.position1, self.position2
        p1_id = Position.id_hap_string_from_elements(chrom=p1.chrom, hap=phasing.value[0], strand=p1.strand, coord=p1.coord)
        p2_id = Position.id_hap_string_from_elements(chrom=p2.chrom, hap=phasing.value[1], strand=p2.strand, coord=p2.coord)
        return "[{p1}]-[{p2}]".format(p1=p1_id, p2=p2_id)

    def get_non_phased_copy(self):
        return self.__class__(position1=self.position1.get_non_hap_copy(),
                              position2=self.position2.get_non_hap_copy(),
                              adjacency_type=self.adjacency_type)

    @property
    def is_self_loop_hapl(self):
        return self.position1.stable_id_non_hap == self.position2.stable_id_non_hap


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


class SegmentCopyNumberProfile(object):

    def __init__(self):
        self.records = defaultdict(dict)

    def set_cn_record(self, sid, hap, cn, check_cn_value=False):
        if check_cn_value and cn < 0:
            raise ValueError()
        self.records[sid][hap] = cn

    def set_cn_record_for_segment(self, segment, cn, haplotype=Haplotype.A, check_cn_value=False):
        haplotype = segment.get_haplotype(default=haplotype)
        sid = segment.stable_id_non_hap
        self.set_cn_record(sid=sid, hap=haplotype, cn=cn, check_cn_value=check_cn_value)

    def get_hap_aware_cn_by_seg_and_hap(self, segment, haplotype, default=0):
        return self.get_cn(sid=segment.stable_id_non_hap, haplotype=haplotype, default=default)

    def get_hap_aware_cn_by_hap_aware_seg(self, segment, default=0):
        if not segment.is_haplotype_specific:
            raise ValueError()
        return self.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=segment.haplotype, default=default)

    def get_non_hap_aware_cn_by_seg(self, segment, default=0):
        return self.get_combined_cn(sid=segment.stable_id_non_hap, default=default)

    def get_cn(self, sid, haplotype, default=0):
        if not self.has_record(sid=sid, haplotype=haplotype):
            return default
        return self.records[sid][haplotype]

    def get_combined_cn(self, sid, default=0):
        if not self.has_any_hap_record(sid=sid):
            return default
        result = 0
        for haplotype in Haplotype:
            if haplotype in self.records[sid]:
                result += self.records[sid][haplotype]
        return result

    def has_any_hap_record(self, sid):
        if sid not in self.records:
            return False
        for haplotype in Haplotype:
            if haplotype in self.records[sid]:
                return True
        return False

    def has_record(self, sid, haplotype):
        return sid in self.records and haplotype in self.records[sid]

    def has_record_for_non_hap_aware_seg(self, segment):
        return self.has_any_hap_record(sid=segment.stable_id_non_hap)

    def has_record_for_hap_aware_seg_by_hap_aware_seg(self, segment):
        return self.has_record(sid=segment.stable_id_non_hap, haplotype=segment.haplotype)

    def has_record_for_hap_aware_seg_by_seg_and_hap(self, segment, haplotype):
        return self.has_record(sid=segment.stable_id_non_hap, haplotype=haplotype)

    @classmethod
    def from_genome(cls, genome, default_haplotype=Haplotype.A):
        result = cls()
        for segment in genome.iter_segments():
            haplotype = segment.get_haplotype(default=default_haplotype)
            current_value = result.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=haplotype, default=0)
            result.set_cn_record_for_segment(segment=segment, cn=current_value + 1, haplotype=haplotype)
        return result

    def sid_keys(self):
        return self.records.keys()

    def sid_hap_pairs(self):
        for sid in self.records:
            for hap in self.records[sid]:
                yield sid, hap


class AdjacencyCopyNumberProfile(object):
    def __init__(self):
        self.records = defaultdict(dict)

    def _set_cn_record(self, aid, phasing, cn, check_cn_value=False):
        if check_cn_value and cn < 0:
            raise ValueError()
        self.records[aid][phasing] = cn

    def set_cnt_record_for_adjacency(self, adjacency, cn, phasing=Phasing.AA, check_cn_value=False):
        phasing = adjacency.get_phasing(sort=True, default=phasing)
        aid = adjacency.stable_id_non_phased
        self._set_cn_record(aid=aid, phasing=phasing, cn=cn, check_cn_value=check_cn_value)

    def get_phase_aware_cn_by_adj_and_phasing(self, adjacency, phasing, default=0):
        aid = adjacency.stable_id_non_phased
        phasing = phasing if adjacency.is_sorted_non_phased else flipped_phasing(phasing=phasing)
        return self.get_cn(aid=aid, phasing=phasing, default=default)

    def get_phase_aware_cn_by_phased_adj(self, adjacency, default=0):
        if not adjacency.is_phased:
            raise ValueError()
        aid = adjacency.stable_id_non_phased
        phasing = adjacency.stable_phasing
        return self.get_cn(aid=aid, phasing=phasing, default=default)

    def get_cn(self, aid, phasing, default=0):
        if not self.has_record(aid=aid, phasing=phasing):
            return default
        return self.records[aid][phasing]

    def get_combined_cn(self, aid, default=0):
        if not self.has_any_phase_record(aid=aid):
            return default
        result = 0
        for phasing in Phasing:
            if phasing in self.records[aid]:
                result += self.records[aid][phasing]
        return result

    def has_record(self, aid, phasing):
        return aid in self.records and phasing in self.records[aid]

    def has_any_phase_record(self, aid):
        if aid not in self.records:
            return False
        for phasing in Phasing:
            if phasing in self.records[aid]:
                return True
        return False

    def has_any_record_for_non_phased_adj(self, adjacency):
        return self.has_any_phase_record(aid=adjacency.stable_id_non_phased)

    def record_for_phased_adj_by_phased_adj(self, adjacency):
        if not adjacency.is_phased:
            raise ValueError()
        return self.has_record(aid=adjacency.stable_id_non_phased, phasing=adjacency.stable_phasing)

    def has_record_for_phased_adj_by_adj_and_phasing(self, adjacency, phasing):
        aid = adjacency.stable_id_non_phased
        phasing = phasing if adjacency.is_sorted_non_phased else flipped_phasing(phasing)
        return self.has_record(aid=aid, phasing=phasing)

    @classmethod
    def from_genome(cls, genome, default_phasing=Phasing.AA):
        result = cls()
        for adjacency in genome.iter_adjacencies():
            phasing = adjacency.get_phasing(sort=False, default=default_phasing)
            current_value = result.get_phase_aware_cn_by_adj_and_phasing(adjacency=adjacency, phasing=phasing, default=0)
            result.set_cnt_record_for_adjacency(adjacency=adjacency, cn=current_value + 1, phasing=phasing)
        return result

    def aid_keys(self):
        return self.records.keys()

    def aid_phase_pairs(self):
        for aid in self.records:
            for phasing in self.records[aid]:
                yield aid, phasing

    def haploid_adj_group_present(self, adj_group):
        for adj in adj_group.adjacencies:
            if self.get_combined_cn(aid=adj.stable_id_non_phased, default=0) == 0:
                return False
        return True


class StructureProfile(object):
    def __init__(self, scn_profile=None, acn_profile=None):
        self.scn_profile = scn_profile
        self.acn_profile = acn_profile

    @classmethod
    def from_genome(cls, genome, default_haplotype=Haplotype.A, default_phasing=Phasing.AA):
        result = cls()
        result.scn_profile = SegmentCopyNumberProfile.from_genome(genome=genome, default_haplotype=default_haplotype)
        result.acn_profile = AdjacencyCopyNumberProfile.from_genome(genome=genome, default_phasing=default_phasing)
        return result


class Chromosome(list):
    """just a standard list of segments with a possibility for extended functionality later on"""

    def iter_adjacencies(self, adjacency_type=AdjacencyType.NOVEL, inherit_phasing_from_positions=True):
        lsg, rsg = itertools.tee(self)
        next(rsg, None)
        for ls, rs in zip(lsg, rsg):
            p1, p2 = ls.end_position, rs.start_position
            result = Adjacency(position1=p1, position2=p2, adjacency_type=adjacency_type)
            if inherit_phasing_from_positions and p1.is_haplotype_specific and p2.is_haplotype_specific:
                if p1 > p2:
                    phasing = haplotype_pair_to_phasing(h1=p2.haplotype, h2=p1.haplotype)
                else:
                    phasing = haplotype_pair_to_phasing(h1=p1.haplotype, h2=p2.haplotype)
                result.extra[PHASING] = phasing
            yield result

    def __add__(self, other):
        return self.__class__(list.__add__(self, other))

    def __mul__(self, other):
        return self.__class__(list.__mul__(self, other))

    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        try:
            return self.__class__(result)
        except TypeError:
            return result


class Genome(list):
    """Just a standard list of chromosomes (i.e., lists) with possibility for extended functionality"""

    def iter_segments(self):
        for chromosome in self:
            for segment in chromosome:
                yield segment

    def iter_adjacencies(self, adjacency_type=AdjacencyType.NOVEL,
                         inherit_phasing_from_positions=True):
        for chromosome in self:
            for adjacency in chromosome.iter_adjacencies(adjacency_type=adjacency_type,
                                                         inherit_phasing_from_positions=inherit_phasing_from_positions):
                yield adjacency

    def __add__(self, other):
        return self.__class__(list.__add__(self, other))

    def __mul__(self, other):
        return self.__class__(list.__mul__(self, other))

    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        if isinstance(result, Chromosome):
            return result
        try:
            return self.__class__(result)
        except TypeError:
            return result

    def iter_telomeres(self):
        for chromosome in self:
            ss, es = chromosome[0], chromosome[-1]
            lt = ss.start_position
            rt = es.end_position
            yield lt
            yield rt


class AdjacencyGroup(object):
    def __init__(self, adjacencies, fp=0.0):
        self.adjacencies = adjacencies
        self.fp = fp

    @property
    def stable_id_non_phased(self):
        return "<>".join([adj.stable_id_non_phased for adj in self.adjacencies])


def get_segments_list_from_genome(genome, copy=True, make_all_non_reversed=True):
    result = []
    for s in genome.iter_segments():
        if copy:
            v = deepcopy(s)
        if v.is_reversed and make_all_non_reversed:
            reverse_segment(segment=v, copy=False)
        result.append(v)
    return result


def get_telomeres_from_genome(genome, copy=True, inherit_haplotypes=True):
    result = []
    for chromosome in genome:
        ss, es = chromosome[0], chromosome[-1]
        lt = ss.start_position
        if copy:
            lt = deepcopy(ss.start_position)
        rt = es.end_position
        if copy:
            rt = deepcopy(es.end_position)
        if HAPLOTYPE in ss.extra and inherit_haplotypes:
            lt.extra[HAPLOTYPE] = ss.extra[HAPLOTYPE]
        if HAPLOTYPE in es.extra and inherit_haplotypes:
            rt.extra[HAPLOTYPE] = es.extra[HAPLOTYPE]
        result.append(lt)
        result.append(rt)
    return result


def strip_phasing_from_adjacencies(adjacencies, inplace=True, strip_positions_haplotypes=True, sort=True):
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
        if sort and v.position1 > v.position2:
            v.position1, v.position2 = v.position2, v.position1
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
        result.start_position.extra[HAPLOTYPE] = result.haplotype
        result.end_position.extra[HAPLOTYPE] = result.haplotype
    return result


def propagate_haplotype_chrom_to_positions(chromosome, inplace=True):
    result = Chromosome()
    for s in chromosome:
        v = propagate_haplotype_segment_to_positions(segment=s, inplace=inplace)
        result.append(v)
    return result


def propagate_haplotype_genome_to_positions(genome, inplace=True):
    result = Genome()
    for chromosome in genome:
        v = propagate_haplotype_chrom_to_positions(chromosome=chromosome, inplace=inplace)
        result.append(v)
    return result


def propagate_phasing_adjacency_to_positions(adjacency, inplace=True):
    result = adjacency
    if not inplace:
        result = deepcopy(adjacency)
    if PHASING in result.extra:
        h1, h2 = phasing_to_haplotype_pair(phasing=result.extra[PHASING])
        result.position1.extra[HAPLOTYPE] = h1
        result.position2.extra[HAPLOTYPE] = h2
    return result


def get_segments_copy_number_profile(segments):
    return Counter(segments)


def get_adjacencies_copy_number_profile(adjacencies):
    return Counter(adjacencies)


def get_unique_forward_haploid_segments(genome):
    result = []
    seen = set()
    for segment in genome.iter_segments():
        if segment.stable_id_non_hap in seen:
            continue
        result.append(segment.get_non_hap_copy())
        seen.add(segment.stable_id_non_hap)
    return result


def get_unique_haploid_ref_and_novel_sorted_adjacencies(ref_genome, mut_genomes):
    ref_adjacencies = {}
    for adj in ref_genome.iter_adjacencies(adjacency_type=AdjacencyType.REFERENCE):
        aid = adj.stable_id_non_phased
        if aid not in ref_adjacencies:
            ref_adjacencies[aid] = adj.get_non_phased_copy()
    nov_adjacencies = {}
    for mut_genome in mut_genomes:
        for adj in mut_genome.iter_adjacencies(adjacency_type=AdjacencyType.NOVEL):
            aid = adj.stable_id_non_phased
            if aid not in ref_adjacencies and aid not in nov_adjacencies:
                nov_adjacencies[aid] = adj.get_non_phased_copy()
    return list(ref_adjacencies.values()), list(nov_adjacencies.values())


def get_novel_adjacencies_from_ref_and_mut_genomes(ref_genome, mut_genome):
    ref_adjacencies = {}
    for adj in ref_genome.iter_adjacencies(adjacency_type=AdjacencyType.REFERENCE):
        aid = adj.stable_id_phased
        if aid not in ref_adjacencies:
            ref_adjacencies[aid] = adj
    nov_adjacencies = {}
    for adj in mut_genome.iter_adjacencies(adjacency_type=AdjacencyType.NOVEL):
        aid = adj.stable_id_phased
        if aid not in ref_adjacencies and aid not in nov_adjacencies:
            nov_adjacencies[aid] = adj
    return list(ref_adjacencies.values()), list(nov_adjacencies.values())


def get_shuffled_scnp(segment_copy_number_profile, segments):
    result = SegmentCopyNumberProfile()
    for segment in segments:
        s_id = segment.stable_id_non_hap
        a_cn = segment_copy_number_profile.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A, default=0)
        b_cn = segment_copy_number_profile.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B, default=0)
        flip = np.random.choice([True, False])
        if flip:
            a_cn, b_cn = b_cn, a_cn
        result.set_cn_record(sid=s_id, hap=Haplotype.A, cn=a_cn)
        result.set_cn_record(sid=s_id, hap=Haplotype.B, cn=b_cn)
    return result


def get_shuffled_clone_specific_scnp(clone_specific_scnp, segments):
    result = {}
    clone_ids = sorted(clone_specific_scnp.keys())
    for clone_id in clone_ids:
        result[clone_id] = SegmentCopyNumberProfile()
    for segment in segments:
        s_id = segment.stable_id_non_hap
        a_cns = {clone_id: clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A, default=0)
                 for clone_id in clone_ids}
        b_cns = {clone_id: clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B, default=0)
                 for clone_id in clone_ids}
        flip = np.random.choice([True, False])
        if flip:
            a_cns, b_cns = b_cns, a_cns
        for clone_id in clone_ids:
            result[clone_id].set_cn_record(sid=s_id, hap=Haplotype.A, cn=a_cns[clone_id])
            result[clone_id].set_cn_record(sid=s_id, hap=Haplotype.B, cn=b_cns[clone_id])
    return result


def add_noise_to_clone_specific_scnp(tensor, segments, min_change=-1, max_change=1):
    clone_ids = sorted(tensor.keys())
    # TODO: finish


def segment_copy_number_tensors_are_compatible(tensor1, tensor2, segments):
    if set(tensor1.keys()) != set(tensor2.keys()):
        return False
    clone_ids = sorted(tensor1.keys())
    for segment in segments:
        t1_a_cns = {clone_id: tensor1[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
        t1_b_cns = {clone_id: tensor1[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
        t2_a_cns = {clone_id: tensor2[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
        t2_b_cns = {clone_id: tensor2[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
        for clone_id in clone_ids:
            compatible = (t1_a_cns[clone_id], t1_b_cns[clone_id]) == (t2_a_cns[clone_id], t2_b_cns[clone_id]) or \
                         (t1_a_cns[clone_id], t1_b_cns[clone_id]) == (t2_b_cns[clone_id], t2_a_cns[clone_id])
            if not compatible:
                return False
    return True


def get_nov_adjacency_group(ref_genome, mut_genome, group_size=5):
    result_dict = {}
    ref_adjacencies = {}
    for adj in ref_genome.iter_adjacencies(adjacency_type=AdjacencyType.REFERENCE):
        aid = adj.stable_id_non_phased
        if aid not in ref_adjacencies:
            ref_adjacencies[aid] = adj.get_non_phased_copy()
    nov_adjacencies = {}
    for adj in mut_genome.iter_adjacencies(adjacency_type=AdjacencyType.NOVEL):
        aid = adj.stable_id_non_phased
        if aid not in ref_adjacencies and aid not in nov_adjacencies:
            nov_adjacencies[aid] = adj.get_non_phased_copy()
    non_phased_novel_adjacencies = list(nov_adjacencies.values())
    group_size = min(len(non_phased_novel_adjacencies), group_size)
    selection = np.random.choice(a=non_phased_novel_adjacencies, size=group_size, replace=False)
    for adj in selection:
        result_dict[adj.stable_id_non_phased] = adj
    result = AdjacencyGroup(adjacencies=list(result_dict.values()))
    return result


def get_unique_haploid_telomeres(genome):
    result = []
    seen = set()
    for tel in genome.iter_telomeres():
        if tel.stable_id_non_hap in seen:
            continue
        result.append(tel.get_non_hap_copy())
    return result


def fragment_id_non_hap_stable_from_segments(segments_span):
    return "<{fsid}>-<{lsid}>".format(fsid=segments_span[0].stable_id_non_hap, lsid=segments_span[-1].stable_id_non_hap)


def check_and_fill_segments_to_fragments(segments, segments_to_fragments):
    if segments_to_fragments is None:
        segments_to_fragments = {}
    for segment in segments:
        sid = segment.stable_id_non_hap
        if sid not in segments_to_fragments:
            segments_to_fragments[sid] = fragment_id_non_hap_stable_from_segments(segments_span=[segment])
    return segments_to_fragments


def segments_to_fragments_reversed(segments_to_fragments):
    result = defaultdict(list)
    for key, value in segments_to_fragments.items():
        result[value].append(key)
    return result


def cn_distance_clone_specific_scnp(tensor1, tensor2, segments, segments_to_fragments=None):
    if set(tensor1.keys()) != set(tensor2.keys()):
        raise Exception()
    clone_ids = sorted(tensor1.keys())
    sids_to_segments = {segment.stable_id_non_hap: segment for segment in segments}
    segments_to_fragments = check_and_fill_segments_to_fragments(segments=segments, segments_to_fragments=segments_to_fragments)
    fragments_to_segments = segments_to_fragments_reversed(segments_to_fragments=segments_to_fragments)
    value = 0
    for fid, sids in fragments_to_segments.items():
        fid_value_ab = 0
        fid_value_ba = 0
        # computing \gamma = 1, \sigma = 2 (i.e., AB for tensor 2)
        for cid in clone_ids:
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                cn_dif_B = abs(tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                fid_value_ab += ((cn_dif_A + cn_dif_B) * segment_length)
        # computing \gamma = 2, \sigma = 1 (i.e, BA for tensor 2)
        for cid in clone_ids:
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                cn_dif_B = abs(tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                fid_value_ba += ((cn_dif_A + cn_dif_B) * segment_length)
        value += min(fid_value_ab, fid_value_ba)
    return value


def get_max_cn_value_from_scnp(segment_copy_number_profile, segments):
    result = -1
    scnp = segment_copy_number_profile
    for segment in segments:
        sid = segment.stable_id_non_hap
        for h in [Haplotype.A, Haplotype.B]:
            cn = scnp.get_cn(sid=sid, haplotype=h, default=0)
            if cn > result:
                result = cn
    return result


def get_max_cn_value_from_clone_specific_scnp(tensor, segments):
    result = -1
    for scnp in tensor.values():
        clone_specific_c_max = get_max_cn_value_from_scnp(segment_copy_number_profile=scnp, segments=segments)
        if clone_specific_c_max > result:
            result = clone_specific_c_max
    return result


REAL_COORDINATE = "real_coordinate"


def split_cluster_to_adhere_for_is(cluster):
    result = []
    current_pos_cnt = 0
    current_neg_cnt = 0
    first_position = cluster.leftmost_position
    current_cluster = PositionCluster(positions=[first_position])
    if first_position.strand == Strand.FORWARD:
        current_pos_cnt = 1
    else:
        current_neg_cnt = 1
    for p in cluster.positions[1:]:
        if p.strand == Strand.FORWARD:
            if current_pos_cnt == 1:
                result.append(current_cluster)
                current_cluster = PositionCluster(positions=[p])
                current_pos_cnt = 1
                current_neg_cnt = 0
            else:
                current_pos_cnt += 1
                current_cluster.positions.append(p)
        else:
            if current_neg_cnt == 1:
                result.append(current_cluster)
                current_cluster = PositionCluster([p])
                current_pos_cnt = 0
                current_neg_cnt = 1
            else:
                current_neg_cnt += 1
                current_cluster.positions.append(p)
    result.append(current_cluster)
    return result


def update_positions_in_cluster(cluster):
    mean_value = int(np.mean([p.coordinate for p in cluster.positions]))
    new_pos_coordinate = mean_value - 1
    new_neg_coordinate = mean_value
    for p in cluster.positions:
        if p.strand == Strand.FORWARD:
            p.extra[REAL_COORDINATE] = p.coordinate
            p.coordinate = new_pos_coordinate
        else:
            p.extra[REAL_COORDINATE] = p.coordinate
            p.coordinate = new_neg_coordinate


def merge_adjacencies_positions_with_window(adjacencies, window_size=300, copy=False, return_clusters_by_chr=True):
    clusters_by_chromosomes = defaultdict(list)
    positions_by_chromosomes = defaultdict(list)
    if copy:
        adjacencies = deepcopy(adjacencies)
    for adjacency in adjacencies:
        p1_chr = adjacency.position1.chromosome
        p2_chr = adjacency.position2.chromosome
        positions_by_chromosomes[p1_chr].append(adjacency.position1)
        positions_by_chromosomes[p2_chr].append(adjacency.position2)
    for chr_name in sorted(positions_by_chromosomes.keys()):
        positions_by_chromosomes[chr_name] = sorted(positions_by_chromosomes[chr_name], key=lambda a: a.coordinate)
    for chr_name in sorted(positions_by_chromosomes.keys()):
        positions = positions_by_chromosomes[chr_name]
        if len(positions) == 0:
            continue
        current_cluster = PositionCluster(positions=[positions[0]])
        current_cluster_added = False
        for position in positions[1:]:
            distance = position.coordinate - current_cluster.rightmost_position.coordinate
            if distance <= window_size:
                current_cluster.positions.append(position)
            else:
                positive_strand_positions = [p for p in current_cluster.positions if p.strand == Strand.FORWARD]
                negative_strand_positions = [p for p in current_cluster.positions if p.strand == Strand.REVERSE]
                if len(positive_strand_positions) > 1 or len(negative_strand_positions) > 1:
                    clusters = split_cluster_to_adhere_for_is(cluster=current_cluster)
                else:
                    clusters = [current_cluster]
                for cluster in clusters:
                    update_positions_in_cluster(cluster=cluster)
                    clusters_by_chromosomes[chr_name].append(cluster)
                    current_cluster_added = True
                current_cluster = PositionCluster([position])
                current_cluster_added = False
        if not current_cluster_added:
            positive_strand_positions = [p for p in current_cluster.positions if p.strand == Strand.FORWARD]
            negative_strand_positions = [p for p in current_cluster.positions if p.strand == Strand.REVERSE]
            if len(positive_strand_positions) > 1 or len(negative_strand_positions) > 1:
                clusters = split_cluster_to_adhere_for_is(cluster=current_cluster)
            else:
                clusters = [current_cluster]
            for cluster in clusters:
                update_positions_in_cluster(cluster=cluster)
                clusters_by_chromosomes[chr_name].append(cluster)
    if return_clusters_by_chr:
        return adjacencies, clusters_by_chromosomes
    return adjacencies


def sorted_segments_donot_overlap(segments):
    for ls, rs in zip(segments[:-1], segments[1:]):
        if ls.end_position.coordinate >= rs.start_position.coordinate:
            return False
    return True


def na_cluster_in_fragment(fragment, cluster):
    return fragment.start_position.coordinate <= cluster.leftmost_position.coordinate and fragment.end_position.coordinate >= cluster.rightmost_position.coordinate


def na_position_clusters_lie_within_fragments(fragments, clusters):
    for cluster in clusters:
        for fragment in fragments:
            if fragment.start_position.coordinate <= cluster.cluster_mean_coordinate <= fragment.end_position.coordinate:
                break
        else:
            return False
    return True


def partition_fragments_into_segments_by_na_clusters(hapl_fragments, na_clusters_by_chr):
    fragments_by_chr = defaultdict(list)
    segments_by_chr = defaultdict(list)
    fragments_to_segments = defaultdict(list)
    for fragment in hapl_fragments:
        chromosome = fragment.chromosome
        fragments_by_chr[chromosome].append(fragment)
    for chr_name in list(fragments_by_chr.keys()):
        fragments_by_chr[chr_name] = sorted(fragments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        if not sorted_segments_donot_overlap(segments=fragments_by_chr[chr_name]):
            raise Exception()
    for chr_name in list(fragments_by_chr.keys()):
        fragments = fragments_by_chr[chr_name]
        if chr_name not in na_clusters_by_chr:
            continue
        na_clusters = na_clusters_by_chr[chr_name]
        suitable = na_position_clusters_lie_within_fragments(fragments=fragments, clusters=na_clusters)
        for fragment in fragments:
            clusters_in_fragment = [c for c in na_clusters if na_cluster_in_fragment(fragment=fragment, cluster=c)]
            clusters = sorted(clusters_in_fragment, key=lambda c: (c.leftmost_position.coordinate, c.rightmost_position.coordinate))
            segments = []
            current_lc = fragment.start_position.coordinate
            for cluster in clusters:
                assert len(cluster.positions) <= 2
                if cluster.leftmost_position.strand == Strand.FORWARD:
                    rc = cluster.leftmost_position.coordinate
                else:
                    rc = cluster.leftmost_position.coordinate - 1
                start_position = Position(chromosome=chr_name, coordinate=current_lc, strand=Strand.REVERSE)
                end_position = Position(chromosome=chr_name, coordinate=rc, strand=Strand.FORWARD)
                segment = Segment(start_position=start_position, end_position=end_position)
                segments.append(segment)
                if cluster.rightmost_position.strand == Strand.FORWARD:
                    current_lc = cluster.rightmost_position.coordinate + 1
                else:
                    current_lc = cluster.rightmost_position.coordinate
            start_position = Position(chromosome=chr_name, coordinate=current_lc, strand=Strand.REVERSE)
            end_position = Position(chromosome=chr_name, coordinate=fragment.end_position.coordinate, strand=Strand.FORWARD)
            segment = Segment(start_position=start_position, end_position=end_position)
            segments.append(segment)
            for segment in segments:
                segments_by_chr[chr_name].append(segment)
                fragments_to_segments[fragment.stable_id_non_hap].append(segment)
    return segments_by_chr, fragments_to_segments


def refined_segments_and_tensor(segments, tensor, merge_fragments, max_merge_gap, fill_gaps, max_fill_gap):
    clone_ids = sorted(tensor.keys())
    new_fragments_by_chr = defaultdict(list)
    fragments_by_chr = defaultdict(list)
    new_tensor = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    for segment in segments:
        chromosome = segment.chromosome
        fragments_by_chr[chromosome].append(segment)
    for chr_name in sorted(fragments_by_chr.keys()):
        fragments_by_chr[chr_name] = sorted(fragments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        if not sorted_segments_donot_overlap(segments=fragments_by_chr[chr_name]):
            raise Exception()
    for chr_name in sorted(fragments_by_chr.keys()):
        chr_fragments = fragments_by_chr[chr_name]
        if len(chr_fragments) == 0:
            continue
        current_new_fragment = deepcopy(chr_fragments[0])
        current_new_f_id = current_new_fragment.stable_id_non_hap
        current_f_cnas = {clone_id: tensor[clone_id].get_cn(sid=current_new_f_id, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
        current_f_cnbs = {clone_id: tensor[clone_id].get_cn(sid=current_new_f_id, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
        for fragment in chr_fragments[1:]:
            distance = fragment.start_position.coordinate - current_new_fragment.end_position.coordinate
            f_cnas = {clone_id: tensor[clone_id].get_cn(sid=fragment.stable_id_non_hap, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
            f_cnbs = {clone_id: tensor[clone_id].get_cn(sid=fragment.stable_id_non_hap, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
            to_merge = merge_fragments and (distance <= max_merge_gap) and \
                       cns_match(cns1=current_f_cnas, cns2=f_cnas, clone_ids=clone_ids) and \
                       cns_match(cns1=current_f_cnbs, cns2=f_cnbs, clone_ids=clone_ids)
            if to_merge:
                current_new_fragment.end_position = fragment.end_position
            else:
                if fill_gaps and (2 <= distance <= max_fill_gap):
                    mid_coordinate = current_new_fragment.end_position.coordinate + int(distance / 2)
                    current_new_fragment.end_position.coordinate = mid_coordinate
                    current_new_f_id = current_new_fragment.stable_id_non_hap
                    fragment.start_position.coordinate = mid_coordinate + 1
                new_fragments_by_chr[chr_name].append(current_new_fragment)
                for clone_id in clone_ids:
                    new_tensor[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.A, cn=current_f_cnas[clone_id])
                    new_tensor[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.B, cn=current_f_cnbs[clone_id])
                current_new_fragment = deepcopy(fragment)
                current_new_f_id = current_new_fragment.stable_id_non_hap
                current_f_cnas = f_cnas
                current_f_cnbs = f_cnbs
        new_fragments_by_chr[chr_name].append(current_new_fragment)
        for clone_id in clone_ids:
            new_tensor[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.A, cn=current_f_cnas[clone_id])
            new_tensor[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.B, cn=current_f_cnbs[clone_id])
    new_fragments = []
    for chr_name in sorted(new_fragments_by_chr.keys()):
        new_fragments.extend(new_fragments_by_chr[chr_name])
    return new_fragments, new_tensor


def refined_nas(novel_adjacencies, window_size, inplace=False):
    if not inplace:
        novel_adjacencies = deepcopy(novel_adjacencies)
    positions_by_chr = defaultdict(list)
    for na in novel_adjacencies:
        p1_chr = na.position1.chromosome
        p2_chr = na.position2.chromosome
        positions_by_chr[p1_chr].append(na.position1)
        positions_by_chr[p2_chr].append(na.position2)
    for chr_name in sorted(positions_by_chr.keys()):
        positions_by_chr[chr_name] = sorted(positions_by_chr[chr_name], key=lambda p: p.coordinate)
    for chr_name in sorted(positions_by_chr.keys()):
        positions = positions_by_chr[chr_name]
        if len(positions) < 2:
            continue
        suitable_for_merging = True
        for lp, rp in zip(positions[:-1], positions[1:]):
            if not suitable_for_merging:
                suitable_for_merging = True
                continue
            if positions_are_reciprocal(p1=lp, p2=rp, window_size=window_size):
                merge_reciprocal_positions(p1=lp, p2=rp)
                suitable_for_merging = False
    return novel_adjacencies


def removed_short_nas(novel_adjacencies, min_size, inplace=False):
    if not inplace:
        novel_adjacencies = deepcopy(novel_adjacencies)
    result = []
    for na in novel_adjacencies:
        if na.position1.chromosome != na.position2.chromosome:
            result.append(na)
        elif abs(na.position1.coordinate - na.position2.coordinate) >= min_size:
            result.append(na)
    return result


def violates_is(novel_adjacencies):
    positions_by_chr_to_nas = defaultdict(lambda: defaultdict(set))
    for na in novel_adjacencies:
        p1_chr = na.position1.chromosome
        p2_chr = na.position2.chromosome
        positions_by_chr_to_nas[p1_chr][(na.position1.coordinate, na.position1.strand)].add(na.stable_id_non_phased)
        positions_by_chr_to_nas[p2_chr][(na.position2.coordinate, na.position2.strand)].add(na.stable_id_non_phased)
    for chr_name in sorted(positions_by_chr_to_nas.keys()):
        for coord, strand in positions_by_chr_to_nas[chr_name].keys():
            nas_ids = positions_by_chr_to_nas[chr_name][(coord, strand)]
            if len(nas_ids) > 1:
                return True
    return False


def positions_are_reciprocal(p1, p2, window_size):
    if p1.chromosome != p2.chromosome:
        return False
    return abs(p1.coordinate - p2.coordinate) <= window_size and p1.strand != p2.strand


def merge_reciprocal_positions(p1, p2):
    lp, rp = (p1, p2) if p1.coordinate < p2.coordinate else (p2, p1)
    fsp, rsp = (p1, p2) if p1.strand == Strand.FORWARD else (p2, p1)
    distance = rp.coordinate - lp.coordinate
    mid_coordinate = lp.coordinate + int(distance / 2)
    fsp.coordinate = mid_coordinate
    rsp.coordinate = mid_coordinate + 1


def cns_match(cns1, cns2, clone_ids):
    for clone_id in clone_ids:
        if cns1[clone_id] != cns2[clone_id]:
            return False
    return True


def refine_segments_and_nas_telomeres(segments, scnt, novel_adjacencies, telomeres=None,
                                      move_fragments_boundaries=False, fragments_boundaries_move_max_distance=1000,
                                      inplace=False):
    fragments = deepcopy(segments)
    if telomeres is None:
        telomeres = []
    refined_segments = []
    clone_ids = sorted(scnt.keys())
    refined_scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    fragments_by_chr = defaultdict(list)
    for fragment in fragments:
        fragments_by_chr[fragment.start_position.chromosome].append(fragment)
    positions_by_chr = defaultdict(list)
    for na in novel_adjacencies:
        p1 = na.position1
        p2 = na.position2
        positions_by_chr[p1.chromosome].append(p1)
        positions_by_chr[p2.chromosome].append(p2)
    for tel_position in telomeres:
        positions_by_chr[tel_position.chromsoome].append(tel_position)
    if not positions_within_segments(segments_by_chr=fragments_by_chr, positions_by_chr=positions_by_chr):
        raise Exception()
    for chr_name in sorted(fragments_by_chr.keys()):
        fragments_by_chr[chr_name] = sorted(fragments_by_chr[chr_name], key=lambda s: (s.start_position, s.end_position))
        if not sorted_segments_donot_overlap(segments=fragments_by_chr[chr_name]):
            raise Exception()
    for chr_name in sorted(fragments_by_chr.keys()):
        processed_positions_ids = set()
        chr_fragments = iter(fragments_by_chr[chr_name])
        chr_nas_positions = iter(sorted(positions_by_chr[chr_name], key=lambda p: p.coordinate))
        current_fragment = next(chr_fragments, None)
        current_position = next(chr_nas_positions, None)
        current_segment = deepcopy(current_fragment)
        refined_segments.append(current_segment)
        while current_fragment is not None and current_position is not None:
            if current_position.stable_id_non_hap in processed_positions_ids:
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate < current_fragment.start_position.coordinate:
                raise Exception()
            elif current_position.coordinate == current_segment.start_position.coordinate and current_position.strand == Strand.REVERSE:
                processed_positions_ids.add(current_position.stable_id_non_hap)
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate <= current_fragment.end_position.coordinate:
                if current_position.strand == Strand.FORWARD:
                    left_partition_coordinate = current_position.coordinate
                else:
                    left_partition_coordinate = current_position.coordinate - 1
                right_partition_coordinate = left_partition_coordinate + 1
                new_end_position = Position(chromosome=chr_name, coordinate=left_partition_coordinate, strand=Strand.FORWARD)
                new_start_position = Position(chromosome=chr_name, coordinate=right_partition_coordinate, strand=Strand.REVERSE)
                current_segment.end_position = new_end_position
                set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
                current_segment = deepcopy(current_fragment)
                current_segment.start_position = new_start_position
                refined_segments.append(current_segment)
                processed_positions_ids.add(current_position.stable_id_non_hap)
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate > current_fragment.end_position.coordinate:
                set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
                current_fragment = next(chr_fragments, None)
                current_segment = deepcopy(current_fragment)
                refined_segments.append(current_segment)
        set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
        current_fragment = next(chr_fragments, None)
        current_segment = deepcopy(current_fragment)
        while current_fragment is not None:
            refined_segments.append(current_segment)
            set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
            current_fragment = next(chr_fragments, None)
            current_segment = deepcopy(current_fragment)
    return fragments, refined_segments, refined_scnt


def positions_within_segments(segments_by_chr, positions_by_chr):
    segments_chr_names = set(segments_by_chr.keys())
    positions_chr_names = set(positions_by_chr.keys())
    if len(positions_chr_names - segments_chr_names) > 0:
        return False
    chr_names = sorted(positions_chr_names & segments_chr_names)
    for chr_name in chr_names:
        segments = iter(sorted(segments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate)))
        positions = iter(sorted(positions_by_chr[chr_name], key=lambda p: p.coordinate))
        current_segment = next(segments, None)
        current_position = next(positions, None)
        while current_segment is not None and current_position is not None:
            if current_position.coordinate < current_segment.start_position.coordinate:
                return False
            elif current_segment.start_position.coordinate <= current_position.coordinate <= current_segment.end_position.coordinate:
                current_position = next(positions, None)
            else:
                current_segment = next(segments, None)
        if current_position is not None:
            return False
    return True


def set_cnr(parent_fragment, fcnt, child_segment, scnt):
    clone_ids = set(fcnt.keys())
    child_clone_ids = set(scnt.keys())
    if len(clone_ids ^ child_clone_ids) != 0:
        raise Exception()
    for clone_id in sorted(clone_ids):
        pcna = fcnt[clone_id].get_cn(sid=parent_fragment.stable_id_non_hap, haplotype=Haplotype.A, default=0)
        pcnb = fcnt[clone_id].get_cn(sid=parent_fragment.stable_id_non_hap, haplotype=Haplotype.B, default=0)
        scnt[clone_id].set_cn_record(sid=child_segment.stable_id_non_hap, hap=Haplotype.A, cn=pcna)
        scnt[clone_id].set_cn_record(sid=child_segment.stable_id_non_hap, hap=Haplotype.B, cn=pcnb)


class SCNBoundariesStrategies(Enum):
    FIXED = "fix"
    UNIFORM_SPREAD = "uniform-spread"
    LENGTH_SPREAD = "length-spread"
    UNIFORM_MIN_MAX = "uniform-min-max"

    @classmethod
    def from_str_value(cls, string_value):
        if string_value == "fix":
            return cls.FIXED
        elif string_value == "uniform-spread":
            return cls.UNIFORM_SPREADя
        elif string_value == "length-spread":
            return cls.LENGTH_SPREAD
        elif string_value == "uniform-min-max":
            return cls.UNIFORM_MIN_MAX
        else:
            raise Exception()


class LengthSpreadRelationships(Enum):
    DUMMY = "dummy"

    @classmethod
    def dummy_lower_boundary(cls, segment, scn):
        s_length = segment.length
        if s_length < 100000:
            return 0
        if s_length < 1000000:
            return max(scn - 3, 0)
        if s_length < 10000000:
            return max(scn - 2, 0)
        return max(scn - 1, 0)

    @classmethod
    def dummy_upper_boundary(cls, segment, scn):
        s_length = segment.length
        if s_length < 100000:
            return 20
        if s_length < 1000000:
            return max(scn + 3, 0)
        if s_length < 10000000:
            return max(scn + 2, 0)
        return max(scn + 1, 0)


def get_lower_scn_bounds(segment, scn, strategy,
                         uniform_spread_size=None,
                         length_spread_relationship=None,
                         uniform_min=None):
    if strategy == SCNBoundariesStrategies.FIXED:
        return max(scn, 0)
    if strategy == SCNBoundariesStrategies.UNIFORM_MIN_MAX:
        if uniform_min is None:
            raise Exception()
        return int(uniform_min)
    if strategy == SCNBoundariesStrategies.UNIFORM_SPREAD:
        if uniform_spread_size is None:
            raise Exception()
        return max(scn - uniform_spread_size, 0)
    if strategy == SCNBoundariesStrategies.LENGTH_SPREAD:
        if length_spread_relationship is None:
            raise Exception()
        return LengthSpreadRelationships.dummy_lower_boundary(segment=segment, scn=scn)
    raise Exception()


def get_upper_scn_bounds(segment, scn, strategy,
                         uniform_spread_size=None,
                         length_spread_relationship=None,
                         uniform_max=None):
    if strategy == SCNBoundariesStrategies.FIXED:
        return max(scn, 0)
    if strategy == SCNBoundariesStrategies.UNIFORM_MIN_MAX:
        if uniform_max is None:
            raise Exception()
        return int(uniform_max)
    if strategy == SCNBoundariesStrategies.UNIFORM_SPREAD:
        if uniform_spread_size is None:
            raise Exception()
        return max(scn + uniform_spread_size, 0)
    if strategy == SCNBoundariesStrategies.LENGTH_SPREAD:
        if length_spread_relationship is None:
            raise Exception()
        return LengthSpreadRelationships.dummy_upper_boundary(segment=segment, scn=scn)
    raise Exception()


def get_scn_boundaries(segments, scnt, strategy,
                       min_allow_zero_for_positive=-1,
                       max_allow_zero_for_positive=1000000000,
                       min_allow_positive_for_zero=-1,
                       max_allow_positive_for_zero=1000000000,
                       uniform_spread_size=None,
                       length_spread_relation=None,
                       uniform_min=None,
                       uniform_max=None,
                       is_female=True,
                       ):
    clone_ids = sorted(scnt.keys())
    result = {clone_id: {s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in segments} for clone_id in clone_ids}
    for clone_id in clone_ids:
        for s in segments:
            sid = s.stable_id_non_hap
            for h in (Haplotype.A, Haplotype.B):
                cn = scnt[clone_id].get_cn(sid=sid, haplotype=h, default=0)
                lower = get_lower_scn_bounds(segment=s, scn=cn, strategy=strategy,
                                             uniform_spread_size=uniform_spread_size,
                                             length_spread_relationship=length_spread_relation,
                                             uniform_min=uniform_min)
                s_length = s.length
                if lower == 0 and cn > 0:
                    if s_length < min_allow_zero_for_positive or s_length > max_allow_zero_for_positive:
                        lower = 1
                if lower > 0 and cn == 0:
                    if s_length < min_allow_positive_for_zero or s_length > max_allow_positive_for_zero:
                        lower = 0
                upper = get_upper_scn_bounds(segment=s, scn=cn, strategy=strategy,
                                             uniform_spread_size=uniform_spread_size,
                                             length_spread_relationship=length_spread_relation,
                                             uniform_max=uniform_max)
                if upper == 0 and cn > 0:
                    if s_length < min_allow_zero_for_positive or s_length > max_allow_zero_for_positive:
                        upper = 1
                if upper > 0 and cn == 0:
                    if s_length < min_allow_positive_for_zero or s_length > max_allow_positive_for_zero:
                        upper = 0
                if "x" in s.chromosome.lower() and not is_female and h == Haplotype.B:
                    lower, upper = 0, 0
                result[clone_id][sid][h] = (lower, upper)
    return result


def boundaries_overlap(fragments, segments):
    fragments_it = iter(fragments)
    segments_it = iter(segments)
    current_fragment = next(fragments_it, None)
    current_segment = next(segments_it, None)
    while current_fragment is not None and current_segment is not None:
        if current_segment.start_position.coordinate < current_fragment.start_position.coordinate < current_segment.end_position.coordinate:
            return True
        if current_segment.start_position.coordinate < current_fragment.end_position.coordinate < current_segment.end_position.coordinate:
            return True
        if current_segment.start_position.coordinate > current_fragment.end_position.coordinate:
            current_fragment = next(fragments_it, None)
        else:
            current_segment = next(segments_it, None)
    return False


def get_segments_for_fragments_ids_dict(segments, fragments, allow_non_covered=True):
    result = dict()
    segments_by_chr = defaultdict(list)
    fragments_by_chr = defaultdict(list)
    for s in segments:
        segments_by_chr[s.chromosome].append(s)
    for f in fragments:
        fragments_by_chr[f.chromosome].append(f)
    for chr_name in sorted(segments_by_chr.keys()):
        if chr_name not in fragments_by_chr and not allow_non_covered:
            raise Exception()
        chr_segments = sorted(segments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        if not sorted_segments_donot_overlap(segments=chr_segments):
            raise Exception()
        chr_fragments = sorted(fragments_by_chr[chr_name], key=lambda f: (f.start_position.coordinate, f.end_position.coordinate))
        if not sorted_segments_donot_overlap(segments=chr_fragments):
            raise Exception()
        if boundaries_overlap(fragments=chr_fragments, segments=chr_segments):
            raise Exception()
        fragments_it = iter(chr_fragments)
        segments_it = iter(chr_segments)
        current_fragment = next(fragments_it, None)
        current_span = []
        current_segment = next(segments_it, None)
        while current_fragment is not None and current_segment is not None:
            while current_segment is not None and current_segment.end_position.coordinate < current_fragment.start_position.coordinate:
                fid = fragment_id_non_hap_stable_from_segments(segments_span=[current_segment])
                result[current_segment.stable_id_non_hap] = fid
                current_segment = next(segments_it, None)
            while current_segment is not None and current_fragment.start_position.coordinate <= current_segment.start_position.coordinate and current_segment.end_position.coordinate <= current_fragment.end_position.coordinate:
                current_span.append(current_segment)
                current_segment = next(segments_it, None)
            if len(current_span) > 0:
                fid = fragment_id_non_hap_stable_from_segments(segments_span=current_span)
                for s in current_span:
                    result[s.stable_id_non_hap] = fid
            current_fragment = next(fragments_it, None)
            current_span = []
        while current_segment is not None:
            fid = fragment_id_non_hap_stable_from_segments(segments_span=[current_segment])
            result[current_segment.stable_id_non_hap] = fid
            current_segment = next(segments_it, None)
    return result


def get_ref_adjacencies_from_segments(segments, assign_external_ids=True):
    from dassp.core.io import EXTERNAL_NA_ID
    current_external_id_cnt = 1
    result = []
    segments_by_chr = defaultdict(list)
    for s in segments:
        segments_by_chr[s.chromosome].append(s)
    for chr_name in sorted(segments_by_chr.keys()):
        segments_by_chr[chr_name] = sorted(segments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        segments = segments_by_chr[chr_name]
        for ls, rs in zip(segments[:-1], segments[1:]):
            adjacency = Adjacency(position1=ls.end_position, position2=rs.start_position, adjacency_type=AdjacencyType.REFERENCE)
            naid = "r{naid}".format(naid=str(current_external_id_cnt))
            adjacency.extra[EXTERNAL_NA_ID] = naid
            current_external_id_cnt += 1
            result.append(adjacency)
    return result


def get_ref_telomeres_from_segments(segments):
    result = []
    segments_by_chr = defaultdict(list)
    for s in segments:
        segments_by_chr[s.chromosome].append(s)
    for chr_name in sorted(segments_by_chr.keys()):
        segments_by_chr[chr_name] = sorted(segments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        segments = segments_by_chr[chr_name]
        result.append(segments[0].start_position)
        result.append(segments[-1].end_position)
    return result


def construct_nas_groups(nas_groups_ids, nas):
    from dassp.core.io import EXTERNAL_NA_ID
    nas_by_external_ids = {na.extra[EXTERNAL_NA_ID]: na for na in nas}
    result = []
    for group_ids in nas_groups_ids:
        group = AdjacencyGroup(adjacencies=[nas_by_external_ids[na_id] for na_id in group_ids])
        result.append(group)
    return result


HUMAN_CENTROMERES = [
    ###
    # chr 1
    ###
    Position(chromosome="1", coordinate=125000000, strand=Strand.FORWARD),
    Position(chromosome="chr1", coordinate=125000000, strand=Strand.FORWARD),
    Position(chromosome="1", coordinate=125000001, strand=Strand.REVERSE),
    Position(chromosome="chr1", coordinate=125000001, strand=Strand.REVERSE),
    ###
    # chr 2
    ###
    Position(chromosome="2", coordinate=93300000, strand=Strand.FORWARD),
    Position(chromosome="chr2", coordinate=93300000, strand=Strand.FORWARD),
    Position(chromosome="2", coordinate=93300001, strand=Strand.REVERSE),
    Position(chromosome="chr2", coordinate=93300001, strand=Strand.REVERSE),
    ###
    # chr 3
    ###
    Position(chromosome="3", coordinate=91000000, strand=Strand.FORWARD),
    Position(chromosome="chr3", coordinate=91000000, strand=Strand.FORWARD),
    Position(chromosome="3", coordinate=91000001, strand=Strand.REVERSE),
    Position(chromosome="chr3", coordinate=91000001, strand=Strand.REVERSE),
    ###
    # chr 4
    ###
    Position(chromosome="4", coordinate=50400000, strand=Strand.FORWARD),
    Position(chromosome="chr4", coordinate=50400000, strand=Strand.FORWARD),
    Position(chromosome="4", coordinate=50400001, strand=Strand.REVERSE),
    Position(chromosome="chr4", coordinate=50400001, strand=Strand.REVERSE),
    ###
    # chr 5
    ###
    Position(chromosome="5", coordinate=48400000, strand=Strand.FORWARD),
    Position(chromosome="chr5", coordinate=48400000, strand=Strand.FORWARD),
    Position(chromosome="5", coordinate=48400001, strand=Strand.REVERSE),
    Position(chromosome="chr5", coordinate=48400001, strand=Strand.REVERSE),
    ###
    # chr 6
    ###
    Position(chromosome="6", coordinate=61000000, strand=Strand.FORWARD),
    Position(chromosome="chr6", coordinate=61000000, strand=Strand.FORWARD),
    Position(chromosome="6", coordinate=61000001, strand=Strand.REVERSE),
    Position(chromosome="chr6", coordinate=61000001, strand=Strand.REVERSE),
    ###
    # chr 7
    ###
    Position(chromosome="7", coordinate=59900000, strand=Strand.FORWARD),
    Position(chromosome="chr7", coordinate=59900000, strand=Strand.FORWARD),
    Position(chromosome="7", coordinate=59900001, strand=Strand.REVERSE),
    Position(chromosome="chr7", coordinate=59900001, strand=Strand.REVERSE),
    ###
    # chr 8
    ###
    Position(chromosome="8", coordinate=45600000, strand=Strand.FORWARD),
    Position(chromosome="chr8", coordinate=45600000, strand=Strand.FORWARD),
    Position(chromosome="8", coordinate=45600001, strand=Strand.REVERSE),
    Position(chromosome="chr8", coordinate=45600001, strand=Strand.REVERSE),
    ###
    # chr 9
    ###
    Position(chromosome="9", coordinate=49000000, strand=Strand.FORWARD),
    Position(chromosome="chr9", coordinate=49000000, strand=Strand.FORWARD),
    Position(chromosome="9", coordinate=49000001, strand=Strand.REVERSE),
    Position(chromosome="chr9", coordinate=49000001, strand=Strand.REVERSE),
    ###
    # chr 10
    ###
    Position(chromosome="10", coordinate=40200000, strand=Strand.FORWARD),
    Position(chromosome="chr10", coordinate=40200000, strand=Strand.FORWARD),
    Position(chromosome="10", coordinate=40200001, strand=Strand.REVERSE),
    Position(chromosome="chr10", coordinate=40200001, strand=Strand.REVERSE),
    ###
    # chr 11
    ###
    Position(chromosome="11", coordinate=53700000, strand=Strand.FORWARD),
    Position(chromosome="chr11", coordinate=53700000, strand=Strand.FORWARD),
    Position(chromosome="11", coordinate=53700001, strand=Strand.REVERSE),
    Position(chromosome="chr11", coordinate=53700001, strand=Strand.REVERSE),
    ###
    # chr 12
    ###
    Position(chromosome="12", coordinate=35800000, strand=Strand.FORWARD),
    Position(chromosome="chr12", coordinate=35800000, strand=Strand.FORWARD),
    Position(chromosome="12", coordinate=35800001, strand=Strand.REVERSE),
    Position(chromosome="chr12", coordinate=35800001, strand=Strand.REVERSE),
    ###
    # chr 13
    ###
    Position(chromosome="13", coordinate=17900000, strand=Strand.FORWARD),
    Position(chromosome="chr13", coordinate=17900000, strand=Strand.FORWARD),
    Position(chromosome="13", coordinate=17900001, strand=Strand.REVERSE),
    Position(chromosome="chr13", coordinate=17900001, strand=Strand.REVERSE),
    ###
    # chr 14
    ###
    Position(chromosome="14", coordinate=17600000, strand=Strand.FORWARD),
    Position(chromosome="chr14", coordinate=17600000, strand=Strand.FORWARD),
    Position(chromosome="14", coordinate=17600001, strand=Strand.REVERSE),
    Position(chromosome="chr14", coordinate=17600001, strand=Strand.REVERSE),
    ###
    # chr 15
    ###
    Position(chromosome="15", coordinate=19000000, strand=Strand.FORWARD),
    Position(chromosome="chr15", coordinate=19000000, strand=Strand.FORWARD),
    Position(chromosome="15", coordinate=19000001, strand=Strand.REVERSE),
    Position(chromosome="chr15", coordinate=19000001, strand=Strand.REVERSE),
    ###
    # chr 16
    ###
    Position(chromosome="16", coordinate=36600000, strand=Strand.FORWARD),
    Position(chromosome="chr16", coordinate=36600000, strand=Strand.FORWARD),
    Position(chromosome="16", coordinate=36600001, strand=Strand.REVERSE),
    Position(chromosome="chr16", coordinate=36600001, strand=Strand.REVERSE),
    ###
    # chr 17
    ###
    Position(chromosome="17", coordinate=24000000, strand=Strand.FORWARD),
    Position(chromosome="chr17", coordinate=24000000, strand=Strand.FORWARD),
    Position(chromosome="17", coordinate=24000001, strand=Strand.REVERSE),
    Position(chromosome="chr17", coordinate=24000001, strand=Strand.REVERSE),
    ###
    # chr 18
    ###
    Position(chromosome="18", coordinate=17200000, strand=Strand.FORWARD),
    Position(chromosome="chr18", coordinate=17200000, strand=Strand.FORWARD),
    Position(chromosome="18", coordinate=17200001, strand=Strand.REVERSE),
    Position(chromosome="chr18", coordinate=17200001, strand=Strand.REVERSE),
    ###
    # chr 19
    ###
    Position(chromosome="19", coordinate=26500000, strand=Strand.FORWARD),
    Position(chromosome="chr19", coordinate=26500000, strand=Strand.FORWARD),
    Position(chromosome="19", coordinate=26500001, strand=Strand.REVERSE),
    Position(chromosome="chr19", coordinate=26500001, strand=Strand.REVERSE),
    ###
    # chr 20
    ###
    Position(chromosome="20", coordinate=27500000, strand=Strand.FORWARD),
    Position(chromosome="chr20", coordinate=27500000, strand=Strand.FORWARD),
    Position(chromosome="20", coordinate=27500001, strand=Strand.REVERSE),
    Position(chromosome="chr20", coordinate=27500001, strand=Strand.REVERSE),
    ###
    # chr 21
    ###
    Position(chromosome="21", coordinate=13200000, strand=Strand.FORWARD),
    Position(chromosome="chr21", coordinate=13200000, strand=Strand.FORWARD),
    Position(chromosome="21", coordinate=13200001, strand=Strand.REVERSE),
    Position(chromosome="chr21", coordinate=13200001, strand=Strand.REVERSE),
    ###
    # chr 22
    ###
    Position(chromosome="22", coordinate=14700000, strand=Strand.FORWARD),
    Position(chromosome="chr22", coordinate=14700000, strand=Strand.FORWARD),
    Position(chromosome="22", coordinate=14700001, strand=Strand.REVERSE),
    Position(chromosome="chr22", coordinate=14700001, strand=Strand.REVERSE),
    ###
    # chr X
    ###
    Position(chromosome="X", coordinate=60600000, strand=Strand.FORWARD),
    Position(chromosome="chrX", coordinate=60600000, strand=Strand.FORWARD),
    Position(chromosome="X", coordinate=60600001, strand=Strand.REVERSE),
    Position(chromosome="chrX", coordinate=60600001, strand=Strand.REVERSE),

]