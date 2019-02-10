# -*- coding: utf-8 -*-
import numpy as np
import math
from collections import defaultdict
from copy import deepcopy
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

    def __lt__(self, other):
        if not isinstance(other, Strand):
            return False
        if self.value != other.value:
            return self.value == "-"
        return False


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

    @staticmethod
    def non_hap_distance_strand_specific(pos1, pos2):
        if pos1.chromosome != pos2.chromosome:
            return 3000000000
        if pos1.strand != pos2.strand:
            return 3000000000
        return abs(pos1.coordinate - pos2.coordinate)


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

    def __init__(self, start_position, end_position, extra=None, idx=None, allow_unit_length=True):
        if start_position.strand != Strand.REVERSE or end_position.strand != Strand.FORWARD:
            raise ValueError("START position has to come from RC strand, "
                             "while END position has to come from the F strand.")
        if start_position.coordinate == end_position.coordinate and not allow_unit_length:
            raise ValueError("START position coordinate has to be less than the END position coordinate")
        if start_position.coordinate > end_position.coordinate:
            raise ValueError("START position coordinate has to be less than the END position coordinate")
        if start_position.chromosome != end_position.chromosome:
            raise ValueError("Segment's START and END position has to come from the same chromosome")
        if start_position.is_haplotype_specific and end_position.is_haplotype_specific and start_position.haplotype != end_position.haplotype:
            raise ValueError("Haplotypes for segment do not match")
        self.start_position = start_position
        self.end_position = end_position
        self.extra = extra if extra is not None else {}
        self._idx = idx
        self.allow_unit_length = allow_unit_length

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
    def from_string(cls, string, allow_unit_length=True):
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
        segment = cls(start_position=start_position, end_position=end_position, extra=extra, allow_unit_length=allow_unit_length)
        if segment_is_reversed:
            segment.start_position, segment.end_position = segment.end_position, segment.start_position
        return segment

    @classmethod
    def from_chromosome_coordinates(cls, chromosome, start, end):
        p1 = Position(chromosome=chromosome, coordinate=start, strand=Strand.REVERSE)
        p2 = Position(chromosome=chromosome, coordinate=end, strand=Strand.FORWARD)
        return cls(start_position=p1, end_position=p2)

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
        return int(self.end_position.coordinate - self.start_position.coordinate + 1)

    @property
    def length_1000(self):
        return int(math.ceil(int(self.end_position.coordinate - self.start_position.coordinate) / 1000.0))

    @property
    def length_100(self):
        return int(math.ceil(int(self.end_position.coordinate - self.start_position.coordinate) / 100.0))

    @property
    def length_10(self):
        return int(math.ceil(int(self.end_position.coordinate - self.start_position.coordinate) / 10.0))

    @property
    def start_coordinate(self):
        return self.start_position.coordinate

    @property
    def end_coordinate(self):
        return self.end_position.coordinate


class AdjacencyType(Enum):
    REFERENCE = "R"
    NOVEL = "N"

    @classmethod
    def from_name(cls, name):
        for adj_type in AdjacencyType:
            if adj_type.value == name:
                return adj_type
        raise ValueError("{} is not a valid adjacency type string".format(name))


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

    @property
    def distance_non_hap(self):
        if self.position1.chromosome != self.position2.chromosome:
            return -1
        return abs(self.position1.coordinate - self.position2.coordinate)


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

    @classmethod
    def combined(cls, *scnps):
        result = cls()
        for scnp in scnps:
            for sid in scnp.records:
                for haplotype in [Haplotype.A, Haplotype.B]:
                    current_cn = result.get_cn(sid=sid, haplotype=haplotype)
                    profile_specific_cn = scnp.get_cn(sid=sid, haplotype=haplotype)
                    result.set_cn_record(sid=sid, hap=haplotype, cn=current_cn + profile_specific_cn)
        return result


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

    def haploid_adjacencies_present(self, adjacencies):
        present = []
        for adj in adjacencies:
            if self.get_combined_cn(aid=adj.stable_id_non_phased, default=0) > 0:
                present.append(adj)
        return present

    @classmethod
    def combined(cls, *acnps):
        result = cls()
        for acnp in acnps:
            for aid in acnp.records.keys():
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    current_cn = result.get_cn(aid=aid, phasing=ph)
                    profile_specific = acnp.get_cn(aid=aid, phasing=ph)
                    result._set_cn_record(aid=aid, phasing=ph, cn=current_cn + profile_specific)
        return result


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


class SCNBoundariesStrategies(Enum):
    FIXED = "fix"
    UNIFORM_SPREAD = "uniform-spread"
    LENGTH_SPREAD = "length-spread"
    UNIFORM_MIN_MAX = "uniform-min-max"

    @classmethod
    def from_string(cls, string_value):
        if string_value == "fix":
            return cls.FIXED
        elif string_value == "uniform-spread":
            return cls.UNIFORM_SPREADя
        elif string_value == "length-spread":
            return cls.LENGTH_SPREAD
        elif string_value == "uniform-min-max":
            return cls.UNIFORM_MIN_MAX
        else:
            raise ValueError("{} is not a value SCNBoundariesStrategy".format(string_value))


class SegmentCopyNumberBoundaries(object):
    def __init__(self):
        self._records = defaultdict(lambda: defaultdict(dict))

    def set_cnb_record(self, sid, hap, boundary_type, value):
        self._records[sid][hap][boundary_type] = value

    def get_cnb(self, sid, hap, boundary_type, default=0):
        if not self.has_record(sid=sid, hap=hap, boundary_type=boundary_type):
            return default
        return self._records[sid][hap][boundary_type]

    def has_record(self, sid, hap, boundary_type):
        if sid not in self._records:
            return False
        if hap not in self._records[sid]:
            return False
        if boundary_type not in self._records[sid][hap]:
            return False
        return True

    def update(self, segment_copy_number_boundaries, overwrite=True):
        for sid in segment_copy_number_boundaries:
            for haplotype in segment_copy_number_boundaries[sid]:
                for boundary_type in segment_copy_number_boundaries[sid][haplotype]:
                    if overwrite or self.has_record(sid=sid, hap=haplotype, boundary_type=boundary_type):
                        self._records[sid][haplotype][boundary_type] = segment_copy_number_boundaries[sid][haplotype][boundary_type]

    def fill(self, segments, scnp, missing_only=True, strategy=SCNBoundariesStrategies.UNIFORM_MIN_MAX,
             min_allow_zero_for_positive=-1, max_allow_zero_for_positive=1000000000,
             min_allow_positive_for_zero=-1, max_allow_positive_for_zero=1000000000,
             uniform_spread_size=None, length_spread_relation=None,
             uniform_min=0, uniform_max=10,
             is_female=True):
        for segment in segments:
            sid = segment.stable_id_non_hap
            for haplotype in [Haplotype.A, Haplotype.B]:
                for boundary_type in [CNBoundaries.LOWER, CNBoundaries.UPPER]:
                    if (not missing_only) or (not self.has_record(sid=sid, hap=haplotype, boundary_type=boundary_type)):
                        scn = scnp.get_cn(sid=sid, haplotype=haplotype)
                        s_length = segment.length
                        if boundary_type == CNBoundaries.LOWER:
                            boundary = get_lower_scn_bounds(segment=segment, scn=scn, strategy=strategy,
                                                            uniform_spread_size=uniform_spread_size, length_spread_relationship=length_spread_relation, uniform_min=uniform_min)
                            if boundary == 0 and scn > 0:
                                if s_length < min_allow_zero_for_positive or s_length > max_allow_zero_for_positive:
                                    boundary = 1
                            if boundary > 0 and scn == 0:
                                if s_length < min_allow_positive_for_zero or s_length > max_allow_positive_for_zero:
                                    boundary = 0
                        else:
                            boundary = get_upper_scn_bounds(segment=segment, scn=scn, strategy=strategy,
                                                            uniform_spread_size=uniform_spread_size, length_spread_relationship=length_spread_relation, uniform_max=uniform_max)
                            if boundary == 0 and scn > 0:
                                if s_length < min_allow_zero_for_positive or s_length > max_allow_zero_for_positive:
                                    boundary = 1
                            if boundary > 0 and scn == 0:
                                if s_length < min_allow_positive_for_zero or s_length > max_allow_positive_for_zero:
                                    boundary = 0
                        if is_female and segment.chromosome.lower()[-1] == "y":
                            boundary = 0
                        elif not is_female and segment.chromosome.lower()[-1] == "x" and Haplotype.B:
                            boundary = 0
                        self.set_cnb_record(sid=sid, hap=haplotype, boundary_type=boundary_type, value=boundary)


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


class AdjacencyGroupType(Enum):
    MOLECULE = "M"
    LABELING = "L"

    @classmethod
    def from_string(cls, string):
        for group_type in cls:
            if group_type.value.lower() == string.lower():
                return group_type
        raise ValueError("{} is not a valid Adjacency Group Type".format(string))


class AdjacencyGroup(object):
    def __init__(self, gid, aids, group_type, adjacencies=None, extra=None):
        self.gid = gid
        self.adjacencies_ids = aids
        self.group_type = group_type
        self.adjacencies = adjacencies if adjacencies is not None else []
        self.extra = extra if extra is not None else {}

    @property
    def stable_id(self):
        return "{gid}:{g_type}:<{aids}>".format(gid=self.gid, g_type=str(self.group_type.value), aids=",".join(self.adjacencies_ids))

    def populate_adjacencies_via_ids(self, source, source_by_ids=None):
        from rck.core.io import EXTERNAL_NA_ID
        if source_by_ids is None:
            source_by_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.idx): adj for adj in source}
        for aid in self.adjacencies_ids:
            if aid not in source_by_ids:
                raise ValueError("Trying to populate adjacencies in adjacency group {gid}, but adjacency {aid} (reference in the group) is missing from available adjacencies"
                                 "".format(gid=self.gid, aid=aid))
            self.adjacencies.append(source_by_ids[aid])


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


def cn_distance_inter_scnt(tensor1, tensor2, segments, segments_to_fragments=None, check_clone_ids_match=True):
    if check_clone_ids_match and set(tensor1.keys()) != set(tensor2.keys()):
        raise Exception()
    assert len(sorted(tensor1.keys())) == len(sorted(tensor2.keys()))
    clone_ids_source = sorted(tensor1.keys())
    clone_ids_target = sorted(tensor2.keys())
    sids_to_segments = {segment.stable_id_non_hap: segment for segment in segments}
    segments_to_fragments = check_and_fill_segments_to_fragments(segments=segments, segments_to_fragments=segments_to_fragments)
    fragments_to_segments = segments_to_fragments_reversed(segments_to_fragments=segments_to_fragments)
    result = {clone_id: 0 for clone_id in clone_ids_target}
    for fid, sids in fragments_to_segments.items():
        clone_specific_fid_ab = {clone_id: 0 for clone_id in clone_ids_target}
        clone_specific_fid_ba = {clone_id: 0 for clone_id in clone_ids_target}
        fid_value_ab = 0
        fid_value_ba = 0
        # computing \gamma = 1, \sigma = 2 (i.e., AB for tensor 2)
        for cid_source, cid_target in zip(clone_ids_source, clone_ids_target):
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                cn_dif_B = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                value = ((cn_dif_A + cn_dif_B) * segment_length)
                clone_specific_fid_ab[cid_target] += value
                fid_value_ab += value
        # computing \gamma = 2, \sigma = 1 (i.e, BA for tensor 2)
        for cid_source, cid_target in zip(clone_ids_source, clone_ids_target):
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                cn_dif_B = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                value = ((cn_dif_A + cn_dif_B) * segment_length)
                clone_specific_fid_ba[cid_target] += value
                fid_value_ba += value
        assert fid_value_ab == sum([clone_specific_fid_ab[clone_id] for clone_id in clone_ids_target])
        assert fid_value_ba == sum([clone_specific_fid_ba[clone_id] for clone_id in clone_ids_target])
        if fid_value_ab <= fid_value_ba:
            for clone_id in clone_ids_target:
                result[clone_id] += clone_specific_fid_ab[clone_id]
        else:
            for clone_id in clone_ids_target:
                result[clone_id] += clone_specific_fid_ba[clone_id]
    return result


def cn_distance_intra_scnt(tensor, segments):
    clone_ids = sorted(tensor.keys())
    pairs = list(itertools.combinations(clone_ids, 2))
    result = defaultdict(int)
    for pair in pairs:
        clone1, clone2 = pair
        clone1_scnp = tensor[clone1]
        clone2_scnp = tensor[clone2]
        result[pair] = cn_pairwise_distance(scnp1=clone1_scnp, scnp2=clone2_scnp, segments=segments)
    return result


def cn_pairwise_distance(scnp1, scnp2, segments):
    result = 0
    for s in segments:
        sid = s.stable_id_non_hap
        clone1_cna = scnp1.get_cn(sid=sid, haplotype=Haplotype.A, default=0)
        clone1_cnb = scnp1.get_cn(sid=sid, haplotype=Haplotype.B, default=0)
        clone2_cna = scnp2.get_cn(sid=sid, haplotype=Haplotype.A, default=0)
        clone2_cnb = scnp2.get_cn(sid=sid, haplotype=Haplotype.B, default=0)
        dif_a = abs(clone1_cna - clone2_cna)
        dif_b = abs(clone1_cnb - clone2_cnb)
        length = s.length
        result += (dif_a * length)
        result += (dif_b * length)
    return result


def f_score(precision, recall):
    if precision == 0 and recall == 0:
        return 0
    return 2.0 * (1.0 * precision * recall) / (precision + recall)


def precision_inf(cn_set1, cn_set2):
    return len(cn_set1 & cn_set2) / len(cn_set2)


def recall_inf(cn_set1, cn_set2):
    return len(cn_set1 & cn_set2) / len(cn_set1)


def cn_accuracy_scnt(tensor1, tensor2, segments):
    clone_ids_source = sorted(tensor1.keys())
    clone_ids_target = sorted(tensor2.keys())
    sids_to_segments = {segment.stable_id_non_hap: segment for segment in segments}
    # segments_to_fragments = check_and_fill_segments_to_fragments(segments=segments, segments_to_fragments=segments_to_fragments)
    # fragments_to_segments = segments_to_fragments_reversed(segments_to_fragments=segments_to_fragments)
    sids = [s.stable_id_non_hap for s in segments]
    precision = 0
    recall = 0
    total_length = 0
    for sid in sids:
        segment = sids_to_segments[sid]
        total_length += segment.length
        source_cns = set()
        for cid in clone_ids_source:
            cna = tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0)
            cnb = tensor1[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0)
            source_cns.add((cna, cnb))
        target_ab_cns = set()
        target_ba_cns = set()
        for cid in clone_ids_target:
            cna = tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.A, default=0)
            cnb = tensor2[cid].get_cn(sid=sid, haplotype=Haplotype.B, default=0)
            target_ab_cns.add((cna, cnb))
            target_ba_cns.add((cnb, cna))
        ab_precision = precision_inf(cn_set1=source_cns, cn_set2=target_ab_cns)
        ab_recall = recall_inf(cn_set1=source_cns, cn_set2=target_ab_cns)
        ab_f_score = f_score(precision=ab_precision, recall=ab_recall)
        ba_precision = precision_inf(cn_set1=source_cns, cn_set2=target_ba_cns)
        ba_recall = recall_inf(cn_set1=source_cns, cn_set2=target_ba_cns)
        ba_f_score = f_score(precision=ba_precision, recall=ba_recall)
        if ab_f_score > ba_f_score:
            precision += ab_precision * segment.length
            recall += ab_recall * segment.length
        else:
            precision += ba_precision * segment.length
            recall += ba_recall * segment.length
    return precision / total_length, recall / total_length


def scnt_length(tensor, segments):
    clone_ids = sorted(tensor.keys())
    result = {clone_id: 0 for clone_id in clone_ids}
    for clone_id in clone_ids:
        for segment in segments:
            sid = segment.stable_id_non_hap
            total_cn = 0
            total_cn += tensor[clone_id].get_cn(sid=sid, haplotype=Haplotype.A, default=0)
            total_cn += tensor[clone_id].get_cn(sid=sid, haplotype=Haplotype.B, default=0)
            result[clone_id] += (total_cn * segment.length)
    return result


def scnts_lengths(scnts_by_name, segments):
    result = {}
    for name, scnt in scnts_by_name.items():
        result[name] = scnt_length(tensor=scnt, segments=segments)
    return result


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


def refined_scnt(segments, scnt, merge_fragments=True, max_merge_gap=1000000000, fill_gaps=True, max_fill_gap=1000000000):
    if not merge_fragments and not fill_gaps:
        return deepcopy(segments), deepcopy(scnt)
    clone_ids = sorted(scnt.keys())
    new_fragments_by_chr = defaultdict(list)
    fragments_by_chr = defaultdict(list)
    new_scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    segments_ids_mapping = defaultdict(list)
    for segment in segments:
        chromosome = segment.chromosome
        fragments_by_chr[chromosome].append(segment)
    for chr_name in sorted(fragments_by_chr.keys()):
        fragments_by_chr[chr_name] = sorted(fragments_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
        if not sorted_segments_donot_overlap(segments=fragments_by_chr[chr_name]):
            raise ValueError("Some segments overlap on chromosome {chr_name}.".format(chr_name=chr_name))
    for chr_name in sorted(fragments_by_chr.keys()):
        chr_fragments = fragments_by_chr[chr_name]
        if len(chr_fragments) == 0:
            continue
        current_new_fragment = deepcopy(chr_fragments[0])
        current_new_f_id = current_new_fragment.stable_id_non_hap
        current_f_cnas = {clone_id: scnt[clone_id].get_cn(sid=current_new_f_id, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
        current_f_cnbs = {clone_id: scnt[clone_id].get_cn(sid=current_new_f_id, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
        old_to_new = [current_new_f_id]
        for fragment in chr_fragments[1:]:
            distance = fragment.start_position.coordinate - current_new_fragment.end_position.coordinate
            f_cnas = {clone_id: scnt[clone_id].get_cn(sid=fragment.stable_id_non_hap, haplotype=Haplotype.A, default=0) for clone_id in clone_ids}
            f_cnbs = {clone_id: scnt[clone_id].get_cn(sid=fragment.stable_id_non_hap, haplotype=Haplotype.B, default=0) for clone_id in clone_ids}
            to_merge = merge_fragments and (distance <= max_merge_gap) and \
                       cns_match(cns1=current_f_cnas, cns2=f_cnas, clone_ids=clone_ids) and \
                       cns_match(cns1=current_f_cnbs, cns2=f_cnbs, clone_ids=clone_ids)
            if to_merge:
                old_to_new.append(fragment.stable_id_non_hap)
                current_new_fragment.end_position = fragment.end_position
            else:
                fill = fill_gaps and (2 <= distance <= max_fill_gap)
                mid_coordinate = current_new_fragment.end_position.coordinate + int(distance / 2)

                if fill:
                    current_new_fragment.end_position.coordinate = mid_coordinate

                current_new_f_id = current_new_fragment.stable_id_non_hap
                for old_id in old_to_new:
                    segments_ids_mapping[current_new_f_id].append(old_id)
                    segments_ids_mapping[old_id].append(current_new_f_id)
                old_to_new = [fragment.stable_id_non_hap]

                if fill:
                    fragment.start_position.coordinate = mid_coordinate + 1

                new_fragments_by_chr[chr_name].append(current_new_fragment)
                current_new_f_id = current_new_fragment.stable_id_non_hap
                for clone_id in clone_ids:
                    new_scnt[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.A, cn=current_f_cnas[clone_id])
                    new_scnt[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.B, cn=current_f_cnbs[clone_id])
                current_new_fragment = deepcopy(fragment)
                current_f_cnas = f_cnas
                current_f_cnbs = f_cnbs
        ####
        # Reaching this place only after everything on the chromosome is processed
        ####
        current_new_f_id = current_new_fragment.stable_id_non_hap
        for old_id in old_to_new:
            segments_ids_mapping[current_new_f_id].append(old_id)
            segments_ids_mapping[old_id].append(current_new_f_id)
        new_fragments_by_chr[chr_name].append(current_new_fragment)
        for clone_id in clone_ids:
            new_scnt[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.A, cn=current_f_cnas[clone_id])
            new_scnt[clone_id].set_cn_record(sid=current_new_f_id, hap=Haplotype.B, cn=current_f_cnbs[clone_id])
    new_segments = []
    for chr_name in sorted(new_fragments_by_chr.keys()):
        new_segments.extend(new_fragments_by_chr[chr_name])
    return new_segments, new_scnt, segments_ids_mapping


def refined_scnb(scnb, new_segments, segments_ids_mapping, allow_missing=True):
    clone_ids = sorted(scnb.keys())
    result = {clone_id: SegmentCopyNumberBoundaries() for clone_id in clone_ids}
    for clone_id in clone_ids:
        old_scnb = scnb[clone_id]
        new_scnb = result[clone_id]
        for new_segment in new_segments:
            new_segment_id = new_segment.stable_id_non_hap
            old_segment_ids = segments_ids_mapping[new_segment_id]
            for haplotype in [Haplotype.A, Haplotype.B]:
                for boundary_type in [CNBoundaries.LOWER, CNBoundaries.UPPER]:
                    for old_segment_id in old_segment_ids:
                        if old_scnb.has_record(sid=old_segment_id, hap=haplotype, boundary_type=boundary_type):
                            new_scnb.set_cnb_record(sid=new_segment_id, hap=haplotype, boundary_type=boundary_type,
                                                    value=old_scnb.get_cnb(sid=old_segment_id, hap=haplotype, boundary_type=boundary_type))
                            break
                    else:
                        if not allow_missing:
                            raise ValueError("For new segment {new_sid} there was no record in the Segment Copy Number Boundary for any of the old segments {old_segment_ids}"
                                             " for boundary type {boundary_type}, haplotype {hap} in clone {clone_id}"
                                             "".format(new_sid=new_segment_id, old_segment_ids=",".join(old_segment_ids), boundary_type=boundary_type.value, hap=haplotype.value,
                                                       clone_id=clone_id))
    return result


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


def cns_match(cns1, cns2, clone_ids):
    for clone_id in clone_ids:
        if cns1[clone_id] != cns2[clone_id]:
            return False
    return True


def refined_scnt_with_adjacencies_and_telomeres(segments, scnt, adjacencies=None, telomere_positions=None, allow_unit_segments=True):
    fragments = deepcopy(segments)
    if telomere_positions is None:
        telomere_positions = []
    refined_segments = []
    clone_ids = sorted(scnt.keys())
    refined_scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    fragments_by_chr = defaultdict(list)
    for fragment in fragments:
        fragments_by_chr[fragment.start_position.chromosome].append(fragment)
    positions_by_chr = defaultdict(list)
    segments_ids_mapping = defaultdict(list)
    if adjacencies is None:
        adjacencies = []
    for adj in adjacencies:
        p1 = adj.position1
        p2 = adj.position2
        positions_by_chr[p1.chromosome].append(p1)
        positions_by_chr[p2.chromosome].append(p2)
    for tel_position in telomere_positions:
        positions_by_chr[tel_position.chromosome].append(tel_position)
    if len(list(positions_by_chr.keys())) == 0:
        return fragments, deepcopy(scnt), {f.stable_id_non_hap: [f.stable_id_non_hap] for f in fragments}
    if not positions_within_segments(segments_by_chr=fragments_by_chr, positions_by_chr=positions_by_chr):
        raise ValueError("Either adjacency or telomere positions do not lie within segment")
    for chr_name in sorted(fragments_by_chr.keys()):
        fragments_by_chr[chr_name] = sorted(fragments_by_chr[chr_name], key=lambda s: (s.start_position, s.end_position))
        if not sorted_segments_donot_overlap(segments=fragments_by_chr[chr_name]):
            raise ValueError("Some segments overlap on chromosome {chr_name}.".format(chr_name=chr_name))
    for chr_name in sorted(fragments_by_chr.keys()):
        processed_positions_ids = set()
        chr_fragments = iter(fragments_by_chr[chr_name])
        if chr_name in positions_by_chr:
            chr_nas_positions = iter(sorted(positions_by_chr[chr_name], key=lambda p: (p.coordinate, p.strand)))
        else:
            chr_nas_positions = iter([])
        current_fragment = next(chr_fragments, None)
        current_position = next(chr_nas_positions, None)
        current_segment = deepcopy(current_fragment)
        refined_segments.append(current_segment)
        new_to_old = []
        while current_fragment is not None and current_position is not None:
            if current_position.stable_id_non_hap in processed_positions_ids:
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate < current_fragment.start_position.coordinate:
                raise ValueError("Position {pos_id} lies outside of segments".format(pos_id=current_position.stable_id_non_hap))
            elif current_position.coordinate == current_segment.start_position.coordinate and current_position.strand == Strand.REVERSE:
                processed_positions_ids.add(current_position.stable_id_non_hap)
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate == current_segment.end_position.coordinate and current_position.strand == Strand.FORWARD:
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
                new_to_old.append(current_segment.stable_id_non_hap)
                set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
                current_segment = deepcopy(current_fragment)
                current_segment.start_position = new_start_position
                refined_segments.append(current_segment)
                processed_positions_ids.add(current_position.stable_id_non_hap)
                current_position = next(chr_nas_positions, None)
            elif current_position.coordinate > current_fragment.end_position.coordinate:
                new_to_old.append(current_segment.stable_id_non_hap)
                current_fragment_id = current_fragment.stable_id_non_hap
                for sid in new_to_old:
                    segments_ids_mapping[current_fragment_id].append(sid)
                    segments_ids_mapping[sid].append(current_fragment_id)
                set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
                current_fragment = next(chr_fragments, None)
                current_segment = deepcopy(current_fragment)
                refined_segments.append(current_segment)
                new_to_old = []
        if current_fragment is not None and current_segment is not None:
            new_to_old.append(current_segment.stable_id_non_hap)
            current_fragment_id = current_fragment.stable_id_non_hap
            for sid in new_to_old:
                segments_ids_mapping[current_fragment_id].append(sid)
                segments_ids_mapping[sid].append(current_fragment_id)
            set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
        current_fragment = next(chr_fragments, None)
        current_segment = deepcopy(current_fragment)
        while current_fragment is not None:
            segments_ids_mapping[current_segment].append(current_segment)
            refined_segments.append(current_segment)
            set_cnr(parent_fragment=current_fragment, fcnt=scnt, child_segment=current_segment, scnt=refined_scnt)
            current_fragment = next(chr_fragments, None)
            current_segment = deepcopy(current_fragment)
    return refined_segments, refined_scnt, segments_ids_mapping


def extract_spanned_extremities(source, boundaries):
    result = []
    source_by_chr = defaultdict(list)
    boundaries_by_chr = defaultdict(list)
    for segment in source:
        source_by_chr[segment.chromosome].append(segment)
    for segment in boundaries:
        boundaries_by_chr[segment.chromosome].append(segment)
    for chr_name in source_by_chr.keys():
        if chr_name not in boundaries_by_chr:
            continue

        source_segments = iter(sorted(source_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate)))
        boundaries_segments = iter(sorted(boundaries_by_chr[chr_name], key=lambda s: (s.start_position.coordinate, s.end_position.coordinate)))
        current_source_segment = next(source_segments, None)
        current_boundaries_segment = next(boundaries_segments, None)
        while current_source_segment is not None and current_boundaries_segment is not None:
            if current_source_segment.start_position.coordinate > current_boundaries_segment.end_position.coordinate:
                current_boundaries_segment = next(boundaries_segments, None)
            elif current_source_segment.end_position.coordinate < current_boundaries_segment.start_position.coordinate:
                current_source_segment = next(source_segments, None)
            else:
                if current_boundaries_segment.start_position.coordinate <= current_source_segment.start_position.coordinate <= current_boundaries_segment.end_position.coordinate:
                    result.append(deepcopy(current_source_segment.start_position))
                if current_boundaries_segment.start_position.coordinate <= current_source_segment.end_position.coordinate <= current_boundaries_segment.end_position.coordinate:
                    result.append(deepcopy(current_source_segment.end_position))
                current_source_segment = next(source_segments, None)
    return result


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
        scnt[clone_id].set_cn_record(sid=child_segment.stable_id_non_hap, hap=Haplotype.A,
                                     cn=fcnt[clone_id].get_cn(sid=parent_fragment.stable_id_non_hap, haplotype=Haplotype.A, default=0))
        scnt[clone_id].set_cn_record(sid=child_segment.stable_id_non_hap, hap=Haplotype.B,
                                     cn=fcnt[clone_id].get_cn(sid=parent_fragment.stable_id_non_hap, haplotype=Haplotype.B, default=0))


def aligned_scnts(segments_by_sample_names, scnts_by_sample_names, fill_gaps=True, max_fill_gap=1000000000):
    sample_names = sorted(segments_by_sample_names.keys())
    for sample_name in sample_names:
        if sample_name not in scnts_by_sample_names:
            raise ValueError("Sample {sample} is present w.r.t. segments, and not in SCNTs".format(sample=sample_name))
    result_segments_by_sample_names = {sample_name: [] for sample_name in sample_names}
    result_scnts_by_sample_names = {sample_name: deepcopy(scnts_by_sample_names[sample_name]) for sample_name in sample_names}
    if fill_gaps:
        for sample_name in sample_names:
            segments = segments_by_sample_names[sample_name]
            scnt = result_scnts_by_sample_names[sample_name]
            ref_segments, ref_scnt, _ = refined_scnt(segments=segments, scnt=scnt, merge_fragments=False, max_merge_gap=1000000000, fill_gaps=True, max_fill_gap=max_fill_gap)
            result_segments_by_sample_names[sample_name] = ref_segments
            result_scnts_by_sample_names[sample_name] = ref_scnt
    telomeres_by_sample_names_by_chr = defaultdict(dict)
    segments_by_sample_names_by_chr = defaultdict(lambda: defaultdict(list))
    all_chromosomes = set()
    outer_most_telomeres_by_chr = {}
    for sample_name in sample_names:
        segments = result_segments_by_sample_names[sample_name]
        for segment in segments:
            all_chromosomes.add(segment.chromosome)
            segments_by_sample_names_by_chr[sample_name][segment.chromosome].append(segment)
    all_chromosomes = sorted(all_chromosomes)
    for sample_name in sample_names:
        for chr_name in all_chromosomes:
            if chr_name not in segments_by_sample_names_by_chr[sample_name]:
                continue
            segments = segments_by_sample_names_by_chr[sample_name][chr_name]
            segments = sorted(segments, key=lambda s: (s.start_position.coordinate, s.end_position.coordinate))
            if not sorted_segments_donot_overlap(segments=segments):
                raise ValueError("In {sample_name} on chromosome {chr_name} some segments overlap.".format(sample_name=sample_name, chr_name=chr_name))
            segments_by_sample_names_by_chr[sample_name][chr_name] = segments
            lt = segments[0].start_position
            rt = segments[-1].end_position
            telomeres_by_sample_names_by_chr[sample_name][chr_name] = (lt, rt)
    for chr_name in all_chromosomes:
        left_telomeres = []
        right_telomeres = []
        for sample_name in sample_names:
            if chr_name not in telomeres_by_sample_names_by_chr[sample_name]:
                continue
            lt, rt = telomeres_by_sample_names_by_chr[sample_name][chr_name]
            left_telomeres.append(lt)
            right_telomeres.append(rt)
        outer_most_lt = min(left_telomeres, key=lambda p: p.coordinate)
        outer_most_rt = max(right_telomeres, key=lambda p: p.coordinate)
        outer_most_telomeres_by_chr[chr_name] = (outer_most_lt, outer_most_rt)
    for sample_name in sample_names:
        for chr_name in all_chromosomes:
            if chr_name not in segments_by_sample_names_by_chr[sample_name]:
                continue
            segments = segments_by_sample_names_by_chr[sample_name][chr_name]
            ls = segments[0]
            rs = segments[-1]
            lt, rt = outer_most_telomeres_by_chr[chr_name]
            old_ls = deepcopy(ls)
            old_rs = deepcopy(rs)
            old_ls_id = old_ls.stable_id_non_hap
            old_rs_id = old_rs.stable_id_non_hap
            ls.start_position = lt
            rs.end_position = rt
            new_ls_id = ls.stable_id_non_hap
            new_rs_id = rs.stable_id_non_hap
            scnt = result_scnts_by_sample_names[sample_name]
            for clone_id in sorted(scnt.keys()):
                ls_cna = scnt[clone_id].get_cn(sid=old_ls_id, haplotype=Haplotype.A, default=0)
                ls_cnb = scnt[clone_id].get_cn(sid=old_ls_id, haplotype=Haplotype.B, default=0)
                scnt[clone_id].set_cn_record(sid=new_ls_id, hap=Haplotype.A, cn=ls_cna)
                scnt[clone_id].set_cn_record(sid=new_ls_id, hap=Haplotype.B, cn=ls_cnb)
                rs_cna = scnt[clone_id].get_cn(sid=old_rs_id, haplotype=Haplotype.A, default=0)
                rs_cnb = scnt[clone_id].get_cn(sid=old_rs_id, haplotype=Haplotype.B, default=0)
                scnt[clone_id].set_cn_record(sid=new_rs_id, hap=Haplotype.A, cn=rs_cna)
                scnt[clone_id].set_cn_record(sid=new_rs_id, hap=Haplotype.B, cn=rs_cnb)
                if old_ls_id != new_ls_id and old_ls_id in scnt[clone_id].records:
                    del scnt[clone_id].records[old_ls_id]
                if old_rs_id != new_rs_id and old_rs_id in scnt[clone_id].records:
                    del scnt[clone_id].records[old_rs_id]
    all_positions = {}
    for sample_name in sample_names:
        for chr_name in all_chromosomes:
            if chr_name not in segments_by_sample_names_by_chr[sample_name]:
                continue
            segments = segments_by_sample_names_by_chr[sample_name][chr_name]
            for segment in segments:
                sp_id = segment.start_position.stable_id_non_hap
                ep_id = segment.end_position.stable_id_non_hap
                if sp_id not in all_positions:
                    all_positions[sp_id] = deepcopy(segment.start_position)
                if ep_id not in all_positions:
                    all_positions[ep_id] = deepcopy(segment.end_position)
    all_positions_list = list(all_positions.values())
    for sample_name in sample_names:
        scnt = result_scnts_by_sample_names[sample_name]
        segments = result_segments_by_sample_names[sample_name]
        ref_segments, ref_scnt = refined_scnt_with_adjacencies_and_telomeres(segments=segments, scnt=scnt, adjacencies=[], telomere_positions=all_positions_list)
        result_segments_by_sample_names[sample_name] = ref_segments
        result_scnts_by_sample_names[sample_name] = ref_scnt
    return result_segments_by_sample_names, result_scnts_by_sample_names


class CNBoundaries(Enum):
    LOWER = "l"
    UPPER = "u"

    @classmethod
    def from_string(cls, string):
        for boundary in CNBoundaries:
            if boundary.value == string.lower():
                return boundary
        raise ValueError("{} is not a valid CNBoundary name".format(string))


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

    @classmethod
    def from_string(cls, string):
        for ls_relationship in cls:
            if ls_relationship.value == string.lower():
                return ls_relationship
        raise ValueError("{} is not a valid LengthSpreadRelationship name".format(string))


def get_lower_scn_bounds(segment, scn, strategy,
                         uniform_spread_size=None,
                         length_spread_relationship=None,
                         uniform_min=None):
    sid = segment.stable_id_non_hap
    if strategy == SCNBoundariesStrategies.FIXED:
        return max(scn, 0)
    if strategy == SCNBoundariesStrategies.UNIFORM_MIN_MAX:
        if uniform_min is None:
            raise ValueError("Uniform-min-max strategy is selected for segment {sid}, but no uniform min value is provided".format(sid=sid))
        return int(uniform_min)
    if strategy == SCNBoundariesStrategies.UNIFORM_SPREAD:
        if uniform_spread_size is None:
            raise ValueError("Uniform-spead strategy is selected for segment {sid}, but no uniform spread size value is provided".format(sid=sid))
        return max(scn - uniform_spread_size, 0)
    if strategy == SCNBoundariesStrategies.LENGTH_SPREAD:
        if length_spread_relationship is None:
            raise ValueError("Length-spread strategy is selected for segment {sid}, but no length spread relationship value is provided".format(sid=sid))
        return LengthSpreadRelationships.dummy_lower_boundary(segment=segment, scn=scn)
    raise ValueError("Unknown strategy {strategy}".format(strategy=strategy))


def get_upper_scn_bounds(segment, scn, strategy,
                         uniform_spread_size=None,
                         length_spread_relationship=None,
                         uniform_max=None):
    sid = segment.stable_id_non_hap
    if strategy == SCNBoundariesStrategies.FIXED:
        return max(scn, 0)
    if strategy == SCNBoundariesStrategies.UNIFORM_MIN_MAX:
        if uniform_max is None:
            raise ValueError("Uniform-min-max strategy is selected for segment {sid}, but no uniform max value is provided".format(sid=sid))
        return int(uniform_max)
    if strategy == SCNBoundariesStrategies.UNIFORM_SPREAD:
        if uniform_spread_size is None:
            raise ValueError("Uniform-spead strategy is selected for segment {sid}, but no uniform spread size value is provided".format(sid=sid))
        return max(scn + uniform_spread_size, 0)
    if strategy == SCNBoundariesStrategies.LENGTH_SPREAD:
        if length_spread_relationship is None:
            raise ValueError("Length-spread strategy is selected for segment {sid}, but no length spread relationship value is provided".format(sid=sid))
        return LengthSpreadRelationships.dummy_upper_boundary(segment=segment, scn=scn)
    raise ValueError("Unknown strategy {strategy}".format(strategy=strategy))


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
    from rck.core.io import EXTERNAL_NA_ID
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
            if assign_external_ids and EXTERNAL_NA_ID not in adjacency.extra:
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
        if len(segments) == 0:
            continue
        result.append(segments[0].start_position)
        result.append(segments[-1].end_position)
    return result


def get_aabb_for_ra(haplotype):
    if haplotype == Haplotype.A:
        return Phasing.AA
    elif haplotype == Haplotype.B:
        return Phasing.BB
    else:
        raise Exception()


def get_abba_for_na_and_position(novel_adjacency, position, haplotype):
    left = novel_adjacency.position1.stable_id_non_hap if novel_adjacency.is_sorted_non_phased else novel_adjacency.position2.stable_id_non_hap
    if position.stable_id_non_hap == left:
        if haplotype == Haplotype.A:
            return Phasing.AB
        elif haplotype == Haplotype.B:
            return Phasing.BA
        else:
            raise Exception()
    else:
        if haplotype == Haplotype.A:
            return Phasing.BA
        elif haplotype == Haplotype.B:
            return Phasing.AB
        else:
            raise Exception()