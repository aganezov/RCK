import random
from collections import Counter
from enum import Enum
from typing import Iterable, Dict, Optional, Set, Union, List, Tuple

from rck.core.structures import Strand, ChromosomeType, InsertedSeq


class RearrangementChrTarget(object):
    def __init__(self, chr_name: str, breakage_coordinates: Iterable[int]):
        self.chr_name = chr_name
        self.breakage_coordinates = breakage_coordinates

    def __str__(self) -> str:
        return f"{self.chr_name}:{','.join(map(str, self.breakage_coordinates))}"

    @classmethod
    def from_str(cls, string: str) -> "RearrangementChrTarget":
        """
        example:
            1:150000,250000,560000
            2:150000,250000,560000
            3:150000,250000,560000
        """
        chr_name, breakage_locations = string.strip().split(":")
        breakage_locations = [int(e) for e in breakage_locations.split(",")]
        return cls(chr_name, breakage_locations)


class RearrangementTarget(object):
    def __init__(self, rearr_chr_targets: Iterable[RearrangementChrTarget]):
        self.rearr_chr_targets = rearr_chr_targets

    def __str__(self) -> str:
        return f"{';'.join(map(str, self.rearr_chr_targets))}"

    @classmethod
    def from_str(cls, string: str) -> "RearrangementTarget":
        """
        example:
            1:150000,250000,560000;2:150000,250000,560000;3:150000,250000,560000
        """
        data = string.strip().split(";")
        rearr_chr_targets = [RearrangementChrTarget.from_str(e) for e in data]
        return cls(rearr_chr_targets)


class OrientedIndexedSegment(object):
    def __init__(self, s_index: int, orientation: Strand):
        """

        "index" determines to "broken" chromosomal part from when the rearrangement is applied
            indexing starts with 0
        "orientation" corresponds to whether or not the broken part will be assembled in the rearrangement target as is, or reversed
        """
        self.s_index = s_index
        self.orientation = orientation

    def __str__(self) -> str:
        return f"{str(self.orientation)}{self.s_index}"

    @classmethod
    def from_str(cls, string: str) -> "OrientedIndexedSegment":
        """
        example:
            +0
            -1
            +2
            -3
            +4
        """
        if len(string.strip()) < 2:
            raise ValueError(f"Unsuccessful attempt to parse '{string}' into an OrientedIndexedSegment entry. Must be at least 2 characters long"
                             f" (one for +- strand sign, one for single-digit index (indexes can be any int))")
        strand, s_index_str = string.strip()[0], string.strip()[1:]
        s_index = int(s_index_str)
        if s_index < 0:
            raise ValueError(f"Unsuccessful attempt to parse '{string}' into an OrientedIndexedSegment entry. Index value (following the +- sign) must be a positive integer,"
                             f" {s_index} was supplied")
        return cls(s_index, Strand.from_pm_string(strand))


class RearrangementChrResult(object):
    def __init__(self, oriented_segments: Iterable[Union[OrientedIndexedSegment, InsertedSeq]], chr_type: ChromosomeType):
        self.oriented_segments = oriented_segments
        self.chr_type = chr_type

    def __str__(self) -> str:
        or_segments = list(self.oriented_segments)
        if len(or_segments) == 0:
            return ""
        return f"{','.join(map(str, or_segments))}{str(self.chr_type)}"

    @classmethod
    def from_str(cls, string: str) -> "RearrangementChrResult":
        """
        example:
             +0,-1,+2$
             -1,-3,+4@
             +0,INS123:A:1-60,+1$
        """
        if len(string.strip()) == 0:
            return cls([], chr_type=ChromosomeType.LINEAR)
        if len(string.strip()) < 3:
            raise ValueError(f"Unsuccessful attempt to parse '{string}' into an OrientedIndexedSegment entry. Must be at least 3 characters long"
                             f" (one, last, for chromosome type ($ or @), and at least 2 string entry for at least one OrientedIndexedSegment entry."
                             f"Can be empty (though not advised), or 3+ symbols")
        data, chr_type_str = string.strip()[:-1], string.strip()[-1]
        oriented_segments = []
        for e in data.split(","):
            try:
                os = OrientedIndexedSegment.from_str(e)
                oriented_segments.append(os)
            except ValueError:
                os = InsertedSeq.from_string(e)
                oriented_segments.append(os)
        return cls(oriented_segments, ChromosomeType.from_str(chr_type_str))


class RearrangementResult(object):
    def __init__(self, chr_results: Iterable[RearrangementChrResult]):
        self.chr_results = chr_results

    def __str__(self) -> str:
        return f"{';'.join(map(str, self.chr_results))}"

    @classmethod
    def from_str(cls, string: str) -> "RearrangementResult":
        """
        example:
            +0,-1,+2$;-1,-3,+4@;+0,INS123:A:1-60,+1$
        """
        data = string.strip().split(";")
        return cls([RearrangementChrResult.from_str(e) for e in data])


class Rearrangement(object):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult, extra: Optional[Dict] = None, validate: bool = True):
        """
        Rearrangement with empty targets are allowed (even with extra info in the 4th column). These usually correspond to whole chromosome deletions.
        """
        self.self_id = self_id
        self.target = target
        self.result = result
        self.extra = extra
        if validate:
            non_unique_target_chrs = self.non_unique_target_chrs(self)
            if len(non_unique_target_chrs) > 0:
                raise ValueError(f"Attempted to create rearrangement {str(self)} in which in {str(non_unique_target_chrs)} targeted chromosomes are mentioned more than once.")

    def __str__(self) -> str:
        extra_str = ";".join(str(k) + "=" + str(v) for k, v in self.extra.items()) if self.extra is not None else ""
        return f"{self.self_id}\t{str(self.target)}\t{str(self.result)}\t{extra_str}"

    @classmethod
    def from_str(cls, string: str) -> "Rearrangement":
        """
        example:
            r1  1:150000,250000,560000;2:150000,250000,560000;3:150000,250000,560000    +0,-1,+2$;-1,-3,+4@
        """
        data = string.strip().split()
        if len(data) < 3:
            raise ValueError("Rearrangement string representation has to have at least 3 tab-separated fields (id, rearr_target, rearr_result), with optional 4th for extra info."
                             " Extra info must be in a 'key1=value1;key1=value1' format")
        self_id = str(data[0])
        rearr_target = RearrangementTarget.from_str(data[1])
        rearr_result = RearrangementResult.from_str(data[2])
        if len(data) > 3:
            extra = {}
            for entry in data[3].split(";"):
                key, value = entry.split("=")
                extra[key] = value
        else:
            extra = None
        return cls(self_id, rearr_target, rearr_result, extra)

    @classmethod
    def constr(cls, *args, **kwargs) -> "Rearrangement":
        raise NotImplemented("classmethod `constr` must only be invoked on Rearrangement's subclasses")

    @classmethod
    def non_unique_target_chrs(cls, rearrangement: "Rearrangement") -> Set[str]:
        chr_counter = Counter((rct.chr_name for rct in rearrangement.target.rearr_chr_targets))
        return {chr_name for chr_name, cnt in chr_counter.items() if cnt > 1}


#########
#
# Specific rearrangements implementations
#
#########

class DupRearr(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult, extra=None):
        super().__init__(self_id, target, result, extra)

    @classmethod
    def constr(cls, self_id, chr_name: str, start: int, end: int, target_location: int, orientation: Strand, target_chr_type: ChromosomeType) -> "DupRearr":
        if end <= start:
            raise ValueError(f"Start coordinate {start} is greater or equal to the end coordinate {end}. Must be distinct with start < end")
        if start < target_location < end:
            raise ValueError("Attempted to insert duplicated fragment in the middle of the fragment being duplicated")
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[start, end, target_location])
        m3_or_segment = None
        left_flanking_or_segment = OrientedIndexedSegment(0, Strand.FORWARD)
        if target_location in [start, end] and target_chr_type == ChromosomeType.LINEAR:
            right_flanking_or_segment = OrientedIndexedSegment(2, Strand.FORWARD)
        elif target_chr_type == ChromosomeType.LINEAR:
            right_flanking_or_segment = OrientedIndexedSegment(3, Strand.FORWARD)
        else:
            right_flanking_or_segment = None
        if target_location < start:
            m1_or_segment = OrientedIndexedSegment(2, orientation)
            m2_or_segment = OrientedIndexedSegment(1, Strand.FORWARD)
            m3_or_segment = OrientedIndexedSegment(2, Strand.FORWARD)
        elif target_location > end:
            m1_or_segment = OrientedIndexedSegment(1, Strand.FORWARD)
            m2_or_segment = OrientedIndexedSegment(2, Strand.FORWARD)
            m3_or_segment = OrientedIndexedSegment(1, orientation)
        else:
            m1_or_segment = OrientedIndexedSegment(1, Strand.FORWARD)
            m2_or_segment = OrientedIndexedSegment(1, orientation)
            if target_location == start:
                m1_or_segment, m2_or_segment = m2_or_segment, m1_or_segment
        oriented_segments = [left_flanking_or_segment, m1_or_segment, m2_or_segment]
        if target_location not in [start, end]:
            oriented_segments.append(m3_or_segment)
        if right_flanking_or_segment is not None:
            oriented_segments.append(right_flanking_or_segment)
        rearr_chr_result = RearrangementChrResult(oriented_segments=oriented_segments,
                                                  chr_type=target_chr_type)
        return cls(self_id=self_id, target=RearrangementTarget([rearr_chr_target]), result=RearrangementResult([rearr_chr_result]))


class TandemDupRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult, extra):
        super().__init__(self_id, target, result, extra)

    @classmethod
    def constr(cls, self_id: str, chr_name: str, start: int, end: int, after: bool = True,
               orientation: Strand = Strand.FORWARD, target_chr_type: ChromosomeType = ChromosomeType.LINEAR) -> "DupRearr":
        return DupRearr.constr(self_id=self_id, chr_name=chr_name, start=start, end=end, target_location=end if after else start, orientation=orientation,
                               target_chr_type=target_chr_type)


class DelRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id, chr_name: str, start: int, end: int, target_chr_type: ChromosomeType):
        if end <= start:
            raise ValueError(f"Start coordinate {start} is greater or equal to the end coordinate {end}. Must be distinct with start < end")
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[start, end])
        oriented_segments = [OrientedIndexedSegment(0, Strand.FORWARD)]
        if target_chr_type == ChromosomeType.LINEAR:
            oriented_segments.append(OrientedIndexedSegment(2, Strand.FORWARD))
        rearr_chr_result = RearrangementChrResult(oriented_segments=oriented_segments, chr_type=target_chr_type)
        return cls(self_id=self_id, target=RearrangementTarget([rearr_chr_target]), result=RearrangementResult([rearr_chr_result]))


class InversionRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id, chr_name: str, start: int, end: int, target_chr_type=ChromosomeType.LINEAR):
        if end <= start:
            raise ValueError(f"Start coordinate {start} is greater or equal to the end coordinate {end}. Must be distinct with start < end")
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[start, end])
        oriented_segments = [OrientedIndexedSegment(0, Strand.FORWARD), OrientedIndexedSegment(1, Strand.REVERSE)]
        if target_chr_type == ChromosomeType.LINEAR:
            oriented_segments.append(OrientedIndexedSegment(2, Strand.FORWARD))
        rearr_chr_result = RearrangementChrResult(oriented_segments=oriented_segments, chr_type=target_chr_type)
        return cls(self_id=self_id, target=RearrangementTarget([rearr_chr_target]), result=RearrangementResult([rearr_chr_result]))


class TranslocationType(Enum):
    CROSSING = "C"
    NON_CROSS = "NC"

    @classmethod
    def from_str(cls, string: str) -> "TranslocationType":
        if string not in {v.value for v in TranslocationType}:
            raise ValueError(f"Can not create a translocation type from supplied string '{string}'. "
                             f"Allowed string values are '{','.join(v.value for v in TranslocationType)}'")
        for v in TranslocationType:
            if v.value == string:
                return v


class TraRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def get_lin2cir0_orient_segments(cls, tra_type: TranslocationType) -> Tuple[List[OrientedIndexedSegment], List[OrientedIndexedSegment]]:
        if tra_type == TranslocationType.CROSSING:
            res_chr1_or_segments = [OrientedIndexedSegment(0, Strand.FORWARD), OrientedIndexedSegment(3, Strand.FORWARD)]
            res_chr2_or_segments = [OrientedIndexedSegment(2, Strand.FORWARD), OrientedIndexedSegment(1, Strand.FORWARD)]
        else:
            res_chr1_or_segments = [OrientedIndexedSegment(0, Strand.FORWARD), OrientedIndexedSegment(2, Strand.REVERSE)]
            res_chr2_or_segments = [OrientedIndexedSegment(1, Strand.REVERSE), OrientedIndexedSegment(3, Strand.FORWARD)]
        return res_chr1_or_segments, res_chr2_or_segments

    @classmethod
    def get_lin1cir1_orient_segments(cls, chr1_type: ChromosomeType, chr2_type: ChromosomeType, tra_type: TranslocationType) -> List[OrientedIndexedSegment]:
        if {chr1_type, chr2_type} != {ChromosomeType.LINEAR, ChromosomeType.CIRCULAR}:
            raise ValueError(f"Attempting to create a translocation rearrangement on 1 linear and 1 circular chromosomes, yet types of both input chromosomes are equal.")
        if chr1_type == ChromosomeType.CIRCULAR:
            tra_type = TranslocationType.CROSSING if tra_type == TranslocationType.NON_CROSS else TranslocationType.NON_CROSS
            return cls.get_lin1cir1_orient_segments(chr1_type=ChromosomeType.LINEAR, chr2_type=ChromosomeType.CIRCULAR, tra_type=tra_type)
        res_chr_segments = [OrientedIndexedSegment(0, Strand.FORWARD)]
        if tra_type.CROSSING:
            res_chr_segments.append(OrientedIndexedSegment(2, Strand.REVERSE))
        else:
            res_chr_segments.append(OrientedIndexedSegment(2, Strand.FORWARD))
        res_chr_segments.append(OrientedIndexedSegment(1, Strand.FORWARD))
        return res_chr_segments

    @classmethod
    def get_lin0cir2_orient_segments(cls, tra_type: TranslocationType) -> List[OrientedIndexedSegment]:
        res_chr_segments = [OrientedIndexedSegment(0, Strand.FORWARD)]
        if tra_type.CROSSING:
            res_chr_segments.append(OrientedIndexedSegment(1, Strand.REVERSE))
        else:
            res_chr_segments.append(OrientedIndexedSegment(1, Strand.FORWARD))
        return res_chr_segments

    @classmethod
    def constr(cls, self_id, chr1_name: str, chr2_name: str, chr1_breakage_location: int, chr2_breakage_location: int,
               chr1_type: ChromosomeType, chr2_type: ChromosomeType,
               translocation_type: TranslocationType = TranslocationType.CROSSING) -> "TraRearrangement":
        rearr_chr1_target = RearrangementChrTarget(chr_name=chr1_name, breakage_coordinates=[chr1_breakage_location])
        rearr_chr2_target = RearrangementChrTarget(chr_name=chr2_name, breakage_coordinates=[chr2_breakage_location])
        rearr_target = RearrangementTarget([rearr_chr1_target, rearr_chr2_target])
        if {chr1_type, chr2_type} == {ChromosomeType.LINEAR, ChromosomeType.LINEAR}:
            rearr_chr1_segments, rearr_chr2_segments = cls.get_lin2cir0_orient_segments(tra_type=translocation_type)
            rearr_result = RearrangementResult([RearrangementChrResult(rearr_chr1_segments, ChromosomeType.LINEAR),
                                                RearrangementChrResult(rearr_chr2_segments, ChromosomeType.LINEAR)])
        elif {chr1_type, chr2_type} == {ChromosomeType.LINEAR, ChromosomeType.CIRCULAR}:
            rearr_result = RearrangementResult([RearrangementChrResult(cls.get_lin1cir1_orient_segments(chr1_type, chr2_type, translocation_type),
                                                                       ChromosomeType.LINEAR)])
        else:
            rearr_result = RearrangementResult([RearrangementChrResult(cls.get_lin0cir2_orient_segments(translocation_type), ChromosomeType.CIRCULAR)])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class ChrDupRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id, chr_name: str, target_chr_type: ChromosomeType = ChromosomeType.LINEAR):
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[])
        rearr_chr1_result = RearrangementChrResult([OrientedIndexedSegment(0, Strand.FORWARD)], chr_type=target_chr_type)
        rearr_chr2_result = RearrangementChrResult([OrientedIndexedSegment(0, Strand.FORWARD)], chr_type=target_chr_type)
        rearr_target = RearrangementTarget([rearr_chr_target])
        rearr_result = RearrangementResult([rearr_chr1_result, rearr_chr2_result])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class MultiChrDupRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id: str, chr_names: Iterable[str], target_chr_types: Iterable[ChromosomeType]):
        chr_names = list(chr_names)
        target_chr_types = list(target_chr_types)
        if len(chr_names) != len(target_chr_types):
            raise ValueError(f"Multi-chromosomal duplication event has a miss-match between input chromosomes names {chr_names} and input chromosome types {target_chr_types}. "
                             f"Must be in length with indexes matching entries in the first list and the second list.")
        rearr_chr_targets = [RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[0]) for chr_name in chr_names]
        rearr_target = RearrangementTarget(rearr_chr_targets)
        rearr_results = []
        for cnt, (chr_name, target_chr_types) in enumerate(zip(chr_names, target_chr_types)):
            rearr_results.append(RearrangementChrResult([OrientedIndexedSegment(cnt, Strand.FORWARD)], chr_type=target_chr_types))
            rearr_results.append(RearrangementChrResult([OrientedIndexedSegment(cnt, Strand.FORWARD)], chr_type=target_chr_types))
        rearr_result = RearrangementResult(rearr_results)
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class ChrDelRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult, extra=None):
        super().__init__(self_id, target, result, extra)

    @classmethod
    def constr(cls, self_id: str, chr_name: str):
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[])
        rearr_target = RearrangementTarget([rearr_chr_target])
        rearr_result = RearrangementResult([])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class MultiChrDelRearrangement(Rearrangement):
    def __init__(self, self_id: str, target: RearrangementTarget, result: RearrangementResult, extra=None):
        super().__init__(self_id, target, result, extra)

    @classmethod
    def constr(cls, self_id: str, chr_names: Iterable[str]):
        rearr_chr_targets = [RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[]) for chr_name in chr_names]
        rearr_target = RearrangementTarget(rearr_chr_targets)
        rearr_result = RearrangementResult([])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class SNPRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def const(cls, self_id: str, chr_name: str, location: int, target_chr_type: ChromosomeType, snp_id: Optional[str] = None, snp_id_prefix: str = "SNP"):
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[location, location + 1])
        rearr_target = RearrangementTarget([rearr_chr_target])
        oriented_segments: List[Union[OrientedIndexedSegment, InsertedSeq]] = [OrientedIndexedSegment(0, Strand.FORWARD)]
        if snp_id is None:
            snp_id = snp_id_prefix + str(random.randint(0, 100000))
        snp = InsertedSeq.from_string(f"{snp_id}:A:1-2")
        oriented_segments.append(snp)
        oriented_segments.append(OrientedIndexedSegment(1, Strand.FORWARD))
        if target_chr_type == ChromosomeType.LINEAR:
            oriented_segments.append(OrientedIndexedSegment(2, Strand.FORWARD))
        rearr_chr_result = RearrangementChrResult(oriented_segments, chr_type=target_chr_type)
        rearr_result = RearrangementResult([rearr_chr_result])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class InsRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id: str, chr_name: str, location: int, target_chr_type: ChromosomeType,
               ins_length: int, ins_id: Optional[str] = None, ins_id_prefix: str = "INS"):
        if ins_length < 1:
            raise ValueError(f"Attempted to create an insertion rearrangement with insertion sequence length < 1. Insertion length Must be positive (i.e., >= 1).")
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[location])
        rearr_target = RearrangementTarget([rearr_chr_target])
        oriented_segments: List[Union[OrientedIndexedSegment, InsertedSeq]] = [OrientedIndexedSegment(0, Strand.FORWARD)]
        if ins_id is None:
            ins_id = ins_id_prefix + str(random.randint(0, 100000))
        insertion = InsertedSeq.from_string(f"{ins_id}:A:1-{ins_length}")
        oriented_segments.append(insertion)
        if target_chr_type == ChromosomeType.LINEAR:
            oriented_segments.append(OrientedIndexedSegment(1, Strand.FORWARD))
        rearr_chr_result = RearrangementChrResult(oriented_segments, chr_type=target_chr_type)
        rearr_result = RearrangementResult([rearr_chr_result])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class BFBRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)

    @classmethod
    def constr(cls, self_id: str, chr_name: str, breakage_location: int, target_chr_type: ChromosomeType = ChromosomeType.LINEAR) -> "BFBRearrangement":
        if target_chr_type != ChromosomeType.LINEAR:
            raise ValueError(f"Breakage-Fusion-Bridge (BFB) rearrangement can only be constructed for a liner chromosome. Circular chr {chr_name} supplied.")
        rearr_chr_target = RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=[breakage_location])
        rearr_target = RearrangementTarget([rearr_chr_target])
        oriented_segments: List[OrientedIndexedSegment] = [OrientedIndexedSegment(0, Strand.FORWARD), OrientedIndexedSegment(0, Strand.REVERSE)]
        rearr_chr_result = RearrangementChrResult(oriented_segments, target_chr_type)
        rearr_result = RearrangementResult([rearr_chr_result])
        return cls(self_id=self_id, target=rearr_target, result=rearr_result)


class ChromothripsisRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)


class ChromoplexyRearrangement(Rearrangement):
    def __init__(self, self_id, target: RearrangementTarget, result: RearrangementResult):
        super().__init__(self_id, target, result)
