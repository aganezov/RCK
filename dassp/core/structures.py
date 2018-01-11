# -*- coding: utf-8 -*-
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
        if start_position.is_haplotype_specific and end_position.is_haplotype_specific and start_position.haplotype == end_position.haplotype:
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
        if self.position1 > self.position2:
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
    def __init__(self, adjacencies, idx, fp=0.0):
        self.adjacencies = adjacencies
        self.fp = fp
        self.idx = idx


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
