# -*- coding: utf-8 -*-


from enum import Enum


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
        return self.coordinate == other.coordinate and self.strand == other.strand and self.chromosome == other.chromosome

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
        return "Position(chromosome={chrom}, coordinate={coord}, strand={strand})" \
               "".format(chrom=self.chromosome, coord=self.coordinate, strand=str(self.strand))

    def __str__(self):
        return "{chrom}:{strand}{coord}".format(chrom=self.chromosome,
                                                strand=str(self.strand),
                                                coord=self.coordinate)


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


class Segment(object):
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
            return "{chrom}:{start}-{end}".format(chrom=self.chromosome,
                                                  start=self.start_position.coordinate,
                                                  end=self.end_position.coordinate)
        return self._idx

    @property
    def chromosome(self):
        return self.start_position.chromosome

    @idx.setter
    def idx(self, value):
        self._idx = value

    def __str__(self):
        return self.idx

    def __repr__(self):
        return "Segment(id={idx}, start_position={sp}, end_position={ep}, extra={ex})" \
               "".format(sp=repr(self.start_position), ep=repr(self.end_position),
                         ex=repr(self.extra), idx=self._idx)


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
