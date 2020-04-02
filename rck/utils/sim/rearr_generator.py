import logging
import random
import re
from enum import Enum
from typing import Tuple, Type, Optional, Set, Dict, List

import numpy

from rck.core.structures import Genome, Segment, Strand, ChromosomeType
from rck.utils.sim.rearrangement import Rearrangement, InsRearrangement, DelRearrangement, InversionRearrangement, DupRearr, TandemDupRearrangement, TraRearrangement, \
    MultiChrDupRearrangement, MultiChrDelRearrangement, TranslocationType, BFBRearrangement, OrientedIndexedSegment, RearrangementChrTarget, RearrangementChrResult, \
    RearrangementTarget, RearrangementResult, ChromothripsisRearrangement, \
    ChromoplexyRearrangement, SNPRearrangement

REARRANGEMENT_STRING = "rearr_str"
CNT = "cnt"
ENTRIES = "entries"
TYPE = "type"
REGIONS = "regions"
LENGTH_RANGES = "length_ranges"
LENGTH_DISTRS = "length_distrs"
PROBABILITY = "prob"
ID_PREFIX = "id_prefix"
INS_ORIENTATIONS = "ins_orientations"
CHR1 = "chr1"
CHR2 = "chr2"
CHR_CNT = "chr_cnt"
BREAKAGE_CNT = "breakage_cnt"
DELETION_PROBABILITY = "del_prob"
AMPLIFICATION_PROBABILITY = "ampl_prob"
TRA_TYPE = "tra_type"


class ProbabilitiesDistributionsTypes(Enum):
    UNIFORM = "u"
    NORMAL = "norm"

    @classmethod
    def from_str(cls, string: str) -> "ProbabilitiesDistributionsTypes":
        if string == cls.UNIFORM.value:
            return cls.UNIFORM
        elif string == cls.NORMAL:
            return cls.NORMAL
        else:
            raise ValueError(f"Attempted to create a probability distribution from string '{string}', that does not match possible values "
                             f"'{','.join(e.value for e in cls)}'")


class ProbabilitiesDistributionsRange(object):
    def __init__(self, prob_type: ProbabilitiesDistributionsTypes, values_range: Tuple[int, int]):
        self.prob_type: ProbabilitiesDistributionsTypes = prob_type
        self.values_range: Tuple[int, int] = values_range

    def generate_value(self) -> int:
        if self.prob_type == ProbabilitiesDistributionsTypes.UNIFORM:
            random_func = numpy.random.randint
        else:
            raise ValueError(f"Error trying to generate value")
        return random_func(self.values_range[0], self.values_range[1])


class Constraint(object):

    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def rearrangement_allowed(cls, genome: Genome, rearrangement: Rearrangement, *args, **kwargs) -> bool:
        raise NotImplemented(f"rearrangement_allowed method shall be called on particular Constraint subclass")


class ForbiddenTelomereCreationConstraint(Constraint):
    @classmethod
    def rearrangement_allowed(cls, genome: Genome, rearrangement: Rearrangement, *args, **kwargs) -> bool:
        if isinstance(rearrangement, (InsRearrangement, DelRearrangement, InversionRearrangement)):
            for rearr_chr_target in rearrangement.target.rearr_chr_targets:
                target_chr = genome.chromosomes[rearr_chr_target.chr_name]
                target_locations = set(min(loc, target_chr.length) for loc in rearr_chr_target.breakage_coordinates)
                if len(target_locations.intersection({1, target_chr.length})) > 0:
                    return False
            return True
        if isinstance(rearrangement, DupRearr):
            for rearr_chr_target in rearrangement.target.rearr_chr_targets:
                target_chr = genome.chromosomes[rearr_chr_target.chr_name]
                target_locations = set(min(loc, target_chr.length) for loc in rearr_chr_target.breakage_coordinates)
                if len(target_locations.intersection({1, target_chr.length})) > 0:
                    return False
            return True
        return True


class HaploidISConstraint(Constraint):
    @classmethod
    def rearrangement_allowed(cls, genome: Genome, rearrangement: Rearrangement, *args, **kwargs) -> bool:
        return True


class DiploidISConstraint(Constraint):
    @classmethod
    def rearrangement_allowed(cls, genome: Genome, rearrangement: Rearrangement, *args, **kwargs) -> bool:
        return True


class RearrangementGenerator(object):
    str_rearr_type_mapping = {
        "ins": InsRearrangement,
        "del": DelRearrangement,
        "tandemdup": TandemDupRearrangement,
        "tradup": DupRearr,
        "inv": InversionRearrangement,
        "tra": TraRearrangement,
        "chrdup": MultiChrDupRearrangement,
        "chrdel": MultiChrDelRearrangement,
        "wgd": MultiChrDupRearrangement,
    }

    @classmethod
    def rearr_type_from_str(cls, string: str) -> Type[Rearrangement]:
        if string not in cls.str_rearr_type_mapping:
            raise ValueError(f"Unable to create a rearrangement generator from `{string}`. "
                             f"Allowed rearrangement generators {','.join(cls.str_rearr_type_mapping.keys())}")
        return cls.str_rearr_type_mapping[string]

    def __init__(self, rearr_type: Type[Rearrangement], id_prefix: str = "REAR", *args, **kwargs):
        self.rearr_type: Type[Rearrangement] = rearr_type
        self.id_prefix = id_prefix

    def generate_rearrangement(self, *args, **kwargs) -> Rearrangement:
        raise NotImplemented(f"generate_rearrangement method can not be called on RearrangementGenerator class instance, rather ir must be invoked on specializing subclasses")

    def _generate_rearrangement_attempt(self, *args, **kwargs) -> Rearrangement:
        raise NotImplemented(f"generate_rearrangement method can not be called on RearrangementGenerator class instance, rather ir must be invoked on specializing subclasses")

    @classmethod
    def suitable(cls, rearrangement: Optional[Rearrangement], genome: Genome, constraints: Optional[Set[Type[Constraint]]] = None) -> bool:
        if rearrangement is None:
            return False
        if constraints is None:
            return True
        for constr in constraints:
            if not constr.rearrangement_allowed(genome, rearrangement):
                return False
        return True

    @classmethod
    def constr(cls, data: Dict) -> "RearrangementGenerator":
        raise NotImplemented(f"Constr method can not be called on RearrangementGenerator class instance, rather ir must be invoked on specializing subclasses")

    @classmethod
    def get_random_target_span(cls, regions: List[Segment], genome: Genome, exclude_chromosomes: Optional[Set[str]] = None) -> Segment:
        if exclude_chromosomes is None:
            exclude_chromosomes = set()
        possible_spans: List[Segment] = []
        for region in regions:
            for chr_name, chromosome in genome.chromosomes.items():
                if chr_name in exclude_chromosomes:
                    continue
                if re.match(region.chromosome, chr_name):
                    if region.start_coordinate > chromosome.length:
                        continue
                    end_coordinate = min(chromosome.length, region.end_coordinate)
                    possible_spans.append(Segment.from_chromosome_coordinates(chromosome=chr_name,
                                                                              start=region.start_coordinate, end=end_coordinate))
        if len(possible_spans) < 1:
            raise ValueError(f"Could not find a random target region")
        return list(numpy.random.choice(possible_spans, size=1))[0]

    @classmethod
    def get_random_location_in_span(cls, span: Segment) -> int:
        return numpy.random.randint(span.start_coordinate, span.end_coordinate + 1)

    @classmethod
    def get_random_length_distribution(cls, lengths_distributions: List[Tuple[float, ProbabilitiesDistributionsRange]]) -> ProbabilitiesDistributionsRange:
        return list(numpy.random.choice([e[1] for e in lengths_distributions], size=1, p=[e[0] for e in lengths_distributions]))[0]

    def get_random_rearr_id(self, forbidden_ids: Optional[Set[str]] = None) -> str:
        if forbidden_ids is None:
            forbidden_ids = set()
        result = self.id_prefix + str(random.randint(1, 100000))
        while result in forbidden_ids:
            result = self.id_prefix + str(random.randint(1, 100000))
        return result

    def get_rearr_id(self, rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None) -> str:
        if rearr_id is None:
            return self.get_random_rearr_id(forbidden_ids)
        return rearr_id

    @classmethod
    def parse_range(cls, range_str: str) -> Tuple[int, int]:
        range_data = range_str.split("-")
        start = int(range_data[0])
        if len(range_data) == 1:
            end = start
        elif len(range_data) == 2:
            end = int(range_data[1])
        else:
            raise ValueError(f"Could not parse a region entry '{range_str}'")
        return start, end

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "RearrangementGenerator":
        raise NotImplemented(f"Abstract rearrangement generator can not be created from a yaml dict. "
                             f"Specific subclasses of Rearrangement Generator can be created from yaml dict, though.")

    @classmethod
    def parse_regions(cls, length_range_str: str) -> List[Segment]:
        if len(length_range_str) == 0:
            return [Segment.from_chromosome_coordinates(chromosome=".*", start=1, end=3000000000)]
        chr_specific_ranges = length_range_str.split(";")
        result = []
        for chr_specific_range in chr_specific_ranges:
            if len(chr_specific_range) == 0:
                result.append(Segment.from_chromosome_coordinates(chromosome=".*", start=1, end=3000000000))
                return result
            data = chr_specific_range.split(":")
            chr_name = data[0]
            if len(data) == 1:
                result.append(Segment.from_chromosome_coordinates(chromosome=chr_name, start=1, end=3000000000))
            else:
                ranges_str = data[1].split(",")
                for range_str in ranges_str:
                    start, end = cls.parse_range(range_str)
                    result.append(Segment.from_chromosome_coordinates(chromosome=chr_name, start=start, end=end))
        return result

    @classmethod
    def parse_lengths(cls, length_ranges_str: str, length_distrs_str: Optional[str]) -> List[Tuple[float, ProbabilitiesDistributionsRange]]:
        length_ranges_data = length_ranges_str.split(",")
        if length_distrs_str is not None and len(length_distrs_str) > 0:
            length_dirtrs_probs = [float(e) for e in length_distrs_str.split(",") if len(e) > 0]
        else:
            length_dirtrs_probs = [1.0 / len(length_ranges_data) for _ in length_ranges_data]
        if len(length_ranges_data) != len(length_dirtrs_probs):
            raise ValueError(f"Lengths ranges '{length_ranges_str}' has a different number of entries than lengths distrs entries '{length_distrs_str}'")
        result = []
        for range_str, prob in zip(length_ranges_data, length_dirtrs_probs):
            start, end = cls.parse_range(range_str)
            result.append((prob, ProbabilitiesDistributionsRange(prob_type=ProbabilitiesDistributionsTypes.UNIFORM,
                                                                 values_range=(start, end))))
        return result

    @classmethod
    def get_suitable_chr_names(cls, regions: List[Segment], genome: Genome, exclude_chromosomes: Optional[Set[str]] = None) -> Set[str]:
        if exclude_chromosomes is None:
            exclude_chromosomes = set()
        suitable_chr_names = set()
        for region in regions:
            for chr_name in genome.chromosomes.keys():
                if chr_name in exclude_chromosomes:
                    continue
                if re.match(region.chromosome, chr_name):
                    suitable_chr_names.add(chr_name)
        return suitable_chr_names


class ExactRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[Rearrangement], rearr_string: str, *args, **kwargs):
        super().__init__(rearr_type, *args, **kwargs)
        self.rearr_string = rearr_string

    @classmethod
    def constr(cls, data: Dict) -> "ExactRearrangementGenerator":
        return cls(rearr_type=Rearrangement.__class__, rearr_string=data[REARRANGEMENT_STRING])

    def generate_rearrangement(self, *args, **kwargs) -> Rearrangement:
        return Rearrangement.from_str(string=self.rearr_string)

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "ExactRearrangementGenerator":
        return cls.constr(data=yaml_dict)


class InsRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[InsRearrangement],
                 regions: List[Segment],
                 lengths: List[Tuple[float, ProbabilitiesDistributionsRange]],
                 id_prefix: str = "INS", *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions
        self.lengths = lengths

    def _generate_rearrangement_attempt(self, genome: Genome, ins_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, retry_cnt: int = 100,
                                        constraints: Optional[List[Constraint]] = None, *args, **kwarg) -> Rearrangement:
        target_span: Segment = self.get_random_target_span(regions=self.regions, genome=genome)
        target_location: int = self.get_random_location_in_span(span=target_span)
        target_length_distribution: ProbabilitiesDistributionsRange = self.get_random_length_distribution(self.lengths)
        ins_id = self.get_rearr_id(ins_rearr_id, forbidden_ids)
        target_length: int = target_length_distribution.generate_value()
        target_ins_rearrangement = InsRearrangement.constr(self_id=ins_id, chr_name=target_span.chromosome, location=target_location,
                                                           target_chr_type=genome.chromosomes[target_span.chromosome].chr_type, ins_length=target_length,
                                                           ins_id=ins_id)
        return target_ins_rearrangement

    def generate_rearrangement(self, genome: Genome, ins_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, retry_cnt: int = 100,
                               constraints: Optional[Set[Type[Constraint]]] = None, *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate an insertion. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, ins_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "InsRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "INS")
        regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        lengths = cls.parse_lengths(yaml_dict.get(LENGTH_RANGES, "50-20000"), yaml_dict.get(LENGTH_DISTRS, ""))
        return cls(InsRearrangement, regions=regions, lengths=lengths, id_prefix=id_prefix)


class SNPRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[SNPRearrangement],
                 regions: List[Segment],
                 id_prefix: "SNP", *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions

    def _generate_rearrangement_attempt(self, genome: Genome, ins_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, retry_cnt: int = 100,
                                        constraints: Optional[List[Constraint]] = None, *args, **kwarg) -> Rearrangement:
        pass

    def generate_rearrangement(self, genome: Genome, snp_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, retry_cnt: int = 100,
                               constraints: Optional[Set[Type[Constraint]]] = None, *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a SNP. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, snp_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr


class DelRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[DelRearrangement],
                 regions: List[Segment],
                 lengths: List[Tuple[float, ProbabilitiesDistributionsRange]],
                 id_prefix: str = "DEL", *args, **kwargs):
        super().__init__(rearr_type, *args, **kwargs)
        self.regions = regions
        self.lengths = lengths
        self.id_prefix = id_prefix

    def _generate_rearrangement_attempt(self, genome: Genome, del_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                                        *args, **kwargs) -> Rearrangement:
        target_span: Segment = self.get_random_target_span(regions=self.regions, genome=genome)
        target_location: int = self.get_random_location_in_span(span=target_span)
        target_length_distribution: ProbabilitiesDistributionsRange = self.get_random_length_distribution(self.lengths)
        target_length: int = target_length_distribution.generate_value()
        del_id = self.get_rearr_id(del_rearr_id, forbidden_ids)
        target_del_rearrangement = DelRearrangement.constr(self_id=del_id, chr_name=target_span.chromosome, start=target_location, end=target_location + target_length,
                                                           target_chr_type=genome.chromosomes[target_span.chromosome].chr_type)
        return target_del_rearrangement

    def generate_rearrangement(self, genome: Genome, del_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a deletion. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, del_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "DelRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "DEL")
        regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        lengths = cls.parse_lengths(yaml_dict.get(LENGTH_RANGES, "50-2000000"), yaml_dict.get(LENGTH_DISTRS, ""))
        return cls(DelRearrangement, regions=regions, lengths=lengths, id_prefix=id_prefix)


class DupRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[DupRearr],
                 target_regions: List[Segment],
                 ins_regions: Optional[List[Segment]],
                 lengths: List[Tuple[float, ProbabilitiesDistributionsRange]],
                 inserted_orientations: List[Tuple[float, Strand]],
                 id_prefix: str = "DUP",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.target_regions = target_regions
        self.ins_regions = ins_regions
        self.lengths = lengths
        self.inserted_orientations = inserted_orientations

    @classmethod
    def get_random_suitable_insertion_span_for_dup(cls, dup_region: Segment, genome: Genome, ins_regions: Optional[List[Segment]] = None, retry_cnt: int = 100) -> Segment:
        try_cnt = 0
        if ins_regions is not None:
            suitable_ins_regions = [r for r in ins_regions if re.match(r.chromosome, dup_region.chromosome)]
        else:
            chr_name = f"^{dup_region.chromosome}$"
            suitable_ins_regions = [Segment.from_chromosome_coordinates(chromosome=chr_name, start=dup_region.start_coordinate, end=dup_region.start_coordinate),
                                    Segment.from_chromosome_coordinates(chromosome=chr_name, start=dup_region.end_coordinate, end=dup_region.end_coordinate)]
        if len(suitable_ins_regions) == 0:
            raise ValueError(f"Can not generate a duplication rearrangement, as non of possible insertion regions match the chromosome of a duplicated segment")
        while try_cnt <= retry_cnt:
            try_cnt += 1
            ins_span = cls.get_random_target_span(regions=suitable_ins_regions, genome=genome)
            assert ins_span.chromosome == dup_region.chromosome
            if dup_region.start_coordinate < ins_span.start_coordinate < ins_span.end_coordinate < dup_region.end_coordinate:
                continue
            left_side = None
            if ins_span.start_coordinate <= dup_region.start_coordinate:
                left_side = Segment.from_chromosome_coordinates(chromosome=ins_span.chromosome, start=ins_span.start_coordinate,
                                                                end=min(ins_span.end_coordinate, dup_region.start_coordinate))
            right_side = None
            if ins_span.end_coordinate >= dup_region.end_coordinate:
                right_side = Segment.from_chromosome_coordinates(chromosome=ins_span.chromosome, start=max(dup_region.end_coordinate, ins_span.start_coordinate),
                                                                 end=ins_span.end_coordinate)
            choices = [s for s in (left_side, right_side) if s is not None]
            if len(choices) == 0:
                continue
            return list(numpy.random.choice(choices, size=1))[0]

    def _generate_rearrangement_attempt(self, genome: Genome, dup_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                                        *args, **kwargs) -> Rearrangement:
        dup_segment: Segment = self.get_random_dup_segment(genome=genome)
        target_ins_span: Segment = self.get_random_suitable_insertion_span_for_dup(dup_region=dup_segment, ins_regions=self.ins_regions, genome=genome)
        target_ins_location: int = self.get_random_location_in_span(span=target_ins_span)
        insertion_orientation = list(numpy.random.choice([s[1] for s in self.inserted_orientations], size=1, p=[s[0] for s in self.inserted_orientations]))[0]
        assert dup_segment.chromosome == target_ins_span.chromosome
        assert target_ins_location <= dup_segment.start_coordinate or target_ins_location >= dup_segment.end_coordinate
        dup_id = self.get_rearr_id(dup_rearr_id, forbidden_ids)
        return DupRearr.constr(self_id=dup_id, chr_name=dup_segment.chromosome, start=dup_segment.start_coordinate, end=dup_segment.end_coordinate,
                               target_location=target_ins_location, orientation=insertion_orientation, target_chr_type=genome.chromosomes[dup_segment.chromosome].chr_type)

    def generate_rearrangement(self, genome: Genome, dup_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a deletion. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, dup_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    def get_random_dup_segment(self, genome: Genome) -> Segment:
        target_span: Segment = self.get_random_target_span(regions=self.target_regions, genome=genome)
        dup_start_location: int = self.get_random_location_in_span(span=target_span)
        target_length_distr = self.get_random_length_distribution(self.lengths)
        target_length: int = target_length_distr.generate_value()
        return Segment.from_chromosome_coordinates(chromosome=target_span.chromosome, start=dup_start_location, end=dup_start_location + target_length)

    @classmethod
    def parse_insertion_orientation(cls, ins_orientation_str: str) -> List[Tuple[float, Strand]]:
        if len(ins_orientation_str) == 0:
            return [(1.0, Strand.FORWARD)]
        result = []
        for e in ins_orientation_str.split(";"):
            data = e.split(",")
            result.append((float(data[0]), Strand.from_pm_string(data[1])))
        return result

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "DupRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "DUP")
        target_regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        lengths = cls.parse_lengths(yaml_dict.get(LENGTH_RANGES, "50-20000000"), yaml_dict.get(LENGTH_DISTRS, ""))
        ins_regions = cls.parse_regions(yaml_dict[REGIONS]) if REGIONS in yaml_dict else None
        ins_orientations = cls.parse_insertion_orientation(yaml_dict.get(INS_ORIENTATIONS, ""))
        return cls(DupRearr, target_regions, ins_regions, lengths, ins_orientations, id_prefix)


class TandemDupRearrangementGenerator(DupRearrangementGenerator):
    def __init__(self, rearr_type: Type[DupRearr],
                 target_regions: List[Segment],
                 lengths: List[Tuple[float, ProbabilitiesDistributionsRange]],
                 ins_orientations: List[Tuple[float, Strand]],
                 id_prefix: str = "TDUP",
                 *args, **kwargs):
        super().__init__(rearr_type, target_regions, None, lengths, ins_orientations, id_prefix, *args, **kwargs)

    def generate_rearrangement(self, genome: Genome, dup_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> Rearrangement:
        return super().generate_rearrangement(genome, dup_rearr_id, forbidden_ids, constraints, retry_cnt, *args, **kwargs)

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "TandemDupRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "TDUP")
        target_regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        lengths = cls.parse_lengths(yaml_dict.get(LENGTH_RANGES, "50-20000000"), yaml_dict.get(LENGTH_DISTRS, ""))
        ins_orientations = cls.parse_insertion_orientation(yaml_dict.get(INS_ORIENTATIONS, ""))
        return cls(DupRearr, target_regions, lengths, ins_orientations, id_prefix)


class InvRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[InversionRearrangement],
                 regions: List[Segment],
                 lengths: List[Tuple[float, ProbabilitiesDistributionsRange]],
                 id_prefix: str = "INV",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions
        self.lengths = lengths

    def _generate_rearrangement_attempt(self, genome: Genome, inv_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                                        *args, **kwargs) -> Rearrangement:
        target_span: Segment = self.get_random_target_span(regions=self.regions, genome=genome)
        target_location: int = self.get_random_location_in_span(span=target_span)
        target_length_distribution: ProbabilitiesDistributionsRange = self.get_random_length_distribution(self.lengths)
        target_length: int = target_length_distribution.generate_value()
        inv_id = self.get_rearr_id(inv_rearr_id, forbidden_ids)
        return InversionRearrangement.constr(self_id=inv_id, chr_name=target_span.chromosome, start=target_location, end=target_location + target_length,
                                             target_chr_type=genome.chromosomes[target_span.chromosome].chr_type)

    def generate_rearrangement(self, genome: Genome, inv_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a deletion. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, inv_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "InvRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "INV")
        target_regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        lengths = cls.parse_lengths(yaml_dict.get(LENGTH_RANGES, "50-20000000"), yaml_dict.get(LENGTH_DISTRS, ""))
        return cls(InversionRearrangement, target_regions, lengths, id_prefix)


class WGDRearrGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[MultiChrDupRearrangement], id_prefix: str = "WGD", *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)

    def generate_rearrangement(self, genome: Genome, wgd_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, *args, **kwargs) -> Rearrangement:
        wgd_id = self.get_rearr_id(wgd_rearr_id, forbidden_ids)
        return MultiChrDupRearrangement.constr(self_id=wgd_id, chr_names=[chromosome.name for chromosome in genome.chromosomes.values()],
                                               target_chr_types=[chromosome.chr_type for chromosome in genome.chromosomes.values()])

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "WGDRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "WGD")
        return cls(MultiChrDupRearrangement, id_prefix)


class MultiChrDupRearrGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[MultiChrDupRearrangement],
                 regions: List[Segment], cnt: Tuple[int, int] = (1, 1),
                 id_prefix: str = "MCHRDUP",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions
        self.cnt = cnt

    def generate_rearrangement(self, genome: Genome, mchr_dup_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, *args, **kwargs) -> Rearrangement:
        suitable_chr_names = self.get_suitable_chr_names(self.regions, genome)
        if len(suitable_chr_names) < self.cnt[0]:
            raise ValueError(f"Can not generate multi-chr duplication rearrangement. Suitable chromosomes cnt ={len(suitable_chr_names)},"
                             f" specified range of affected chromosomes {self.cnt[0]}-{self.cnt[1]}")
        suitable_chr_names = list(suitable_chr_names)
        target_chr_cnt = numpy.random.randint(self.cnt[0], min(len(suitable_chr_names), self.cnt[1]) + 1)
        target_chr_names = list(numpy.random.choice(suitable_chr_names, size=target_chr_cnt, replace=False))
        mchr_dup_id = self.get_rearr_id(mchr_dup_rearr_id, forbidden_ids)
        return MultiChrDupRearrangement.constr(self_id=mchr_dup_id, chr_names=[chr_name for chr_name in target_chr_names],
                                               target_chr_types=[genome.chromosomes[chr_name].chr_type for chr_name in target_chr_names])

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "MultiChrDupRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "MCHRDUP")
        target_regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        cnt = cls.parse_range(str(yaml_dict.get(CNT, "1")))
        return cls(MultiChrDupRearrangement, target_regions, cnt, id_prefix)


class MultiChrDelRearrGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[MultiChrDelRearrangement],
                 regions: List[Segment],
                 cnt: Tuple[int, int] = (1, 1),
                 id_prefix: str = "MCHRDEL",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions
        self.cnt = cnt

    def generate_rearrangement(self, genome: Genome, mchr_del_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, *args, **kwargs) -> Rearrangement:
        suitable_chr_names = self.get_suitable_chr_names(self.regions, genome)
        if len(suitable_chr_names) < self.cnt[0]:
            raise ValueError(f"Can not generate multi-chr deletion rearrangement. Suitable chromosomes cnt ={len(suitable_chr_names)},"
                             f" specified range of affected chromosomes {self.cnt[0]}-{self.cnt[1]}")
        suitable_chr_names = list(suitable_chr_names)
        target_chr_cnt = numpy.random.randint(self.cnt[0], min(len(suitable_chr_names), self.cnt[1]) + 1)
        target_chr_names = list(numpy.random.choice(suitable_chr_names, size=target_chr_cnt, replace=False))
        mchr_del_id = self.get_rearr_id(mchr_del_rearr_id, forbidden_ids)
        return MultiChrDelRearrangement.constr(self_id=mchr_del_id, chr_names=[chr_name for chr_name in target_chr_names])

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "MultiChrDelRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "MCHRDEL")
        target_regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        cnt = cls.parse_range(str(yaml_dict.get(CNT, "1")))
        return cls(MultiChrDelRearrangement, target_regions, cnt, id_prefix)


class TraRearrangementGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[TraRearrangement],
                 chr1_regions: List[Segment],
                 chr2_regions: List[Segment],
                 tra_crossings: List[Tuple[float, TranslocationType]],
                 allow_circ_chromosomes: bool,
                 id_prefix: str = "TRA",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.chr1_regions: List[Segment] = chr1_regions
        self.chr2_regions: List[Segment] = chr2_regions
        self.tra_crossings: List[Tuple[float, TranslocationType]] = tra_crossings
        self.allow_circ_chromosomes = allow_circ_chromosomes

    def _generate_rearrangement_attempt(self, genome: Genome, tra_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                                        *args, **kwargs) -> Optional[TraRearrangement]:
        circ_chromosomes_names: Set[str] = set()
        if not self.allow_circ_chromosomes:
            circ_chromosomes_names = {chromosome.name for chromosome in genome.chromosomes.values() if chromosome.chr_type == ChromosomeType.CIRCULAR}
        chr1_span = self.get_random_target_span(regions=self.chr1_regions, genome=genome, exclude_chromosomes=circ_chromosomes_names)
        chr2_suitable_chrs = [chrom for chrom in genome.chromosomes.values() if chrom.name != chr1_span.chromosome]
        if not self.allow_circ_chromosomes:
            chr2_suitable_chrs = [chrom for chrom in chr2_suitable_chrs if chrom.chr_type != ChromosomeType.CIRCULAR]
        if len(chr2_suitable_chrs) == 0:
            return None
        exclude = {chr1_span.chromosome}
        if not self.allow_circ_chromosomes:
            exclude.update(circ_chromosomes_names)
        chr2_span = self.get_random_target_span(regions=self.chr2_regions, genome=genome, exclude_chromosomes=exclude)
        chr1_location: int = self.get_random_location_in_span(span=chr1_span)
        chr2_location: int = self.get_random_location_in_span(span=chr2_span)
        tra_crossing_type = list(numpy.random.choice([e[1] for e in self.tra_crossings], size=1, p=[e[0] for e in self.tra_crossings]))[0]
        tra_id = self.get_rearr_id(tra_rearr_id, forbidden_ids)
        return TraRearrangement.constr(self_id=tra_id, chr1_name=chr1_span.chromosome, chr2_name=chr2_span.chromosome,
                                       chr1_breakage_location=chr1_location, chr2_breakage_location=chr2_location,
                                       chr1_type=genome.chromosomes[chr1_span.chromosome].chr_type, chr2_type=genome.chromosomes[chr2_span.chromosome].chr_type,
                                       translocation_type=tra_crossing_type)

    def generate_rearrangement(self, genome: Genome, tra_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> TraRearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a translocation. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, tra_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def parse_tra_types(cls, tra_types_string: str) -> List[Tuple[float, TranslocationType]]:
        if len(tra_types_string) == 0:
            return [(1.0, TranslocationType.CROSSING)]
        entries = tra_types_string.split(",")
        result = []
        for entry in entries:
            data = entry.split(";")
            if len(data) == 1:
                prob = None
                tra_type = TranslocationType.from_str(data[0])
            else:
                type_str = data[0]
                if len(type_str) == 0:
                    tra_type = TranslocationType.CROSSING
                else:
                    tra_type = TranslocationType.from_str(type_str)
                prob_str = data[1]
                if len(prob_str) == 0:
                    prob = None
                else:
                    prob = float(prob_str)
            result.append((prob, tra_type))
        remaining_prob = 1.0
        specified_prob = [e for e in result if e[0] is not None]
        non_specified_prob = [e for e in result if e[0] is None]
        for entry in result:
            if entry[0] is not None:
                remaining_prob -= entry[0]
        if remaining_prob < 0:
            raise ValueError(f"Total supplied probabilities for Translocation type are > 1.0")
        if len(non_specified_prob) > 0:
            non_specified_prob = [(remaining_prob / len(non_specified_prob), e[1]) for e in non_specified_prob]
        result = specified_prob
        result += non_specified_prob
        return result

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "TraRearrangementGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "TRA")
        chr1_regions = cls.parse_regions(yaml_dict.get(CHR1, ""))
        chr2_regions = cls.parse_regions(yaml_dict.get(CHR2, ""))
        tra_types = cls.parse_tra_types(yaml_dict.get(TRA_TYPE, ""))
        return cls(TraRearrangement, chr1_regions, chr2_regions, tra_types, id_prefix)


class BFBRRearrGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[BFBRearrangement],
                 regions: List[Segment],
                 id_prefix: str = "BFB",
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions = regions

    def _generate_rearrangement_attempt(self, genome: Genome, inv_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None, *args, **kwargs) -> BFBRearrangement:
        circular_chr_names = {chr_name for chr_name, chromosome in genome.chromosomes.items() if chromosome.chr_type == ChromosomeType.CIRCULAR}
        target_span = self.get_random_target_span(self.regions, genome, exclude_chromosomes=circular_chr_names)
        target_location = self.get_random_location_in_span(target_span)
        bfb_id = self.get_rearr_id(inv_rearr_id, forbidden_ids)
        return BFBRearrangement.constr(self_id=bfb_id, chr_name=target_span.chromosome, breakage_location=target_location,
                                       target_chr_type=genome.chromosomes[target_span.chromosome].chr_type)

    def generate_rearrangement(self, genome: Genome, bfb_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100,
                               *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a BFB. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, bfb_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "BFBRRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "BFB")
        regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        return cls(BFBRearrangement, regions, id_prefix)


class KBreakGenerator(RearrangementGenerator):
    def __init__(self, rearr_type: Type[Rearrangement],
                 regions: List[Segment],
                 chr_cnt: Tuple[int, int], breakage_cnt: Tuple[int, int],
                 deletion_prob: float, amplification_prob: float,
                 id_prefix: str = "KBREAK",
                 allow_telomere_translocation: bool = False,
                 *args, **kwargs):
        super().__init__(rearr_type, id_prefix, *args, **kwargs)
        self.regions: List[Segment] = regions
        self.chr_cnt: Tuple[int, int] = chr_cnt
        self.breakage_cnt: Tuple[int, int] = breakage_cnt
        self.deletion_prob: float = deletion_prob
        self.amplification_prob: float = amplification_prob
        self.allow_telomere_translocation = allow_telomere_translocation

    def _generate_rearrangement_attempt(self, genome: Genome, kbreak_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                                        *args, **kwargs) -> Optional[Rearrangement]:
        suitable_chr_names: Set[str] = set()
        for region in self.regions:
            for chromosome in genome.chromosomes.values():
                if chromosome.chr_type == ChromosomeType.CIRCULAR:
                    continue
                if re.match(region.chromosome, chromosome.name):
                    suitable_chr_names.add(chromosome.name)
        if len(suitable_chr_names) < self.chr_cnt[0]:
            raise ValueError(f"Can not generate multi-chr duplication rearrangement. Suitable chromosomes cnt ={len(suitable_chr_names)},"
                             f" specified range of affected chromosomes {self.chr_cnt[0]}-{self.chr_cnt[1]}")
        suitable_chr_names: List[str] = list(suitable_chr_names)
        target_chr_cnt = numpy.random.randint(self.chr_cnt[0], min(len(suitable_chr_names), self.chr_cnt[1]))
        breakage_cnt = numpy.random.randint(self.breakage_cnt[0], self.breakage_cnt[1] + 1)
        target_chr_names = list(numpy.random.choice(suitable_chr_names, size=min(target_chr_cnt, breakage_cnt), replace=False))
        target_chr_breakage_cnt_distr = self.get_chr_breakage_cnt_assignment(chr_names=target_chr_names, total_breakage_cnt=breakage_cnt)
        target_chr_breakage_locations = [(chr_name, self.get_chr_breakage_locations(chr_name, breakage_cnt, genome))
                                         for chr_name, breakage_cnt in target_chr_breakage_cnt_distr.items()]
        for _, locations in target_chr_breakage_locations:
            if locations is None:
                return None
        target_chr_blocks = []
        chr_telomere_blocks = []
        non_telomere_blocks = []
        cnt = 0
        for chr_name, locations in target_chr_breakage_locations:
            chr_blocks = [OrientedIndexedSegment(cnt + x, Strand.FORWARD) for x in range(len(locations) + 1)]
            cnt += len(chr_blocks)
            target_chr_blocks.append((chr_name, chr_blocks))
            chr_telomere_blocks.append((chr_name, (chr_blocks[0], chr_blocks[-1])))
            for block in chr_blocks[1:-1]:
                to_delete = list(numpy.random.choice([True, False], size=1, p=[self.deletion_prob, 1.0 - self.deletion_prob]))[0]
                if to_delete:
                    continue
                to_duplicate = list(numpy.random.choice([True, False], size=1, p=[self.amplification_prob, 1.0 - self.amplification_prob]))[0]
                if to_duplicate:
                    non_telomere_blocks.append(block)
                non_telomere_blocks.append(block)
        new_telomere_pairs = self.get_random_telomere_matching(chr_telomere_blocks)
        new_chromosomes = []
        for cnt, new_telomere_pair in enumerate(new_telomere_pairs):
            new_chromosome = [new_telomere_pair[0]]
            if len(non_telomere_blocks) > 0:
                filling = self.get_random_block_filling(non_telomere_blocks, has_left=cnt != (len(new_telomere_pairs) - 1))
                for block in filling:
                    non_telomere_blocks.remove(block)
                for block in filling:
                    block.orientation = list(numpy.random.choice([Strand.FORWARD, Strand.REVERSE], size=1))[0]
                    new_chromosome.append(block)
            new_chromosome.append(new_telomere_pair[1])
            new_chromosomes.append(new_chromosome)
        rearr_id = self.get_rearr_id(kbreak_rearr_id, forbidden_ids)
        rearr_chr_targets = [RearrangementChrTarget(chr_name=chr_name, breakage_coordinates=coordinates) for chr_name, coordinates in target_chr_breakage_locations]
        rearr_chr_results = [RearrangementChrResult(or_segments, chr_type=ChromosomeType.LINEAR) for or_segments in new_chromosomes]
        return Rearrangement(self_id=rearr_id, target=RearrangementTarget(rearr_chr_targets), result=RearrangementResult(rearr_chr_results))

    def generate_rearrangement(self, genome: Genome, bfb_rearr_id: Optional[str] = None, forbidden_ids: Optional[Set[str]] = None,
                               constraints: Optional[Set[Type[Constraint]]] = None, retry_cnt: int = 100, *args, **kwargs) -> Rearrangement:
        try_cnt = 0
        while True:
            try_cnt += 1
            if try_cnt >= retry_cnt:
                raise ValueError(f"Tried {try_cnt} number of times to generate a k-break. Failed to do so with constraints.")
            rearr = self._generate_rearrangement_attempt(genome, bfb_rearr_id, forbidden_ids, retry_cnt, constraints)
            if self.suitable(rearr, genome, constraints):
                return rearr

    @classmethod
    def get_random_block_filling(cls, remaining_blocks: List[OrientedIndexedSegment], has_left=True):
        filling = list(
            numpy.random.choice(remaining_blocks, size=numpy.random.randint(low=1 if has_left else len(remaining_blocks), high=len(remaining_blocks) + 1), replace=False))
        return filling

    @classmethod
    def get_random_telomere_matching(cls, chr_telomere_blocks: List[Tuple[str, Tuple[OrientedIndexedSegment, OrientedIndexedSegment]]]
                                     ) -> List[Tuple[OrientedIndexedSegment, OrientedIndexedSegment]]:
        left_telomeres: List[OrientedIndexedSegment] = []
        right_telomeres: List[OrientedIndexedSegment] = []
        for chr_name, (left, right) in chr_telomere_blocks:
            left_telomeres.append(left)
            right_telomeres.append(right)
        result = []
        sources = [left_telomeres, right_telomeres]
        while len(left_telomeres) > 0 or len(right_telomeres) > 0:
            low = 0 if len(left_telomeres) > 0 else 1
            high = 2 if len(right_telomeres) > 0 else 1
            new_left_source_index = numpy.random.randint(low, high)
            new_left_source = sources[new_left_source_index]
            new_left: OrientedIndexedSegment = list(numpy.random.choice(new_left_source, size=1))[0]
            new_left_source.remove(new_left)
            if new_left_source_index == 1:
                new_left.orientation = Strand.REVERSE
            low = 0 if len(left_telomeres) > 0 else 1
            high = 2 if len(right_telomeres) > 0 else 1
            new_right_source_index = numpy.random.randint(low, high)
            new_right_source = sources[new_right_source_index]
            new_right: OrientedIndexedSegment = list(numpy.random.choice(new_right_source, size=1))[0]
            if new_right_source_index == 0:
                new_right.orientation = Strand.REVERSE
            new_right_source.remove(new_right)
            result.append((new_left, new_right))
        return result

    def get_chr_breakage_locations(self, target_chr: str, target_breakage_cnt: int, genome: Genome) -> Optional[List[int]]:
        result = set()
        suitable_regions = [s for s in self.regions if re.match(s.chromosome, target_chr) is not None]
        seen_spans = set()
        while len(result) < target_breakage_cnt:
            if len(seen_spans) == len(suitable_regions):
                return None
            for region in suitable_regions:
                if len(result) >= target_breakage_cnt:
                    break
                chromosome = genome.chromosomes[target_chr]
                if region.start_coordinate > chromosome.length:
                    continue
                end_coordinate = min(chromosome.length, region.end_coordinate)
                span = Segment.from_chromosome_coordinates(chromosome=chromosome.name,
                                                           start=region.start_coordinate, end=end_coordinate)
                seen_spans.add(span)
                can_extract_cnt = span.length_100 + 1
                for value in numpy.random.randint(low=span.start_coordinate, high=span.end_coordinate + 1, size=can_extract_cnt):
                    if value not in result:
                        result.add(value)
                        if len(result) == target_breakage_cnt:
                            break
        return sorted(result)

    @classmethod
    def get_chr_breakage_cnt_assignment(cls, chr_names: List[set], total_breakage_cnt: int) -> Dict[str, int]:
        if total_breakage_cnt < len(chr_names):
            raise ValueError(f"Total number of breaks {total_breakage_cnt} is less than the number {len(chr_names)} of chromosomes on which those breaks must occur")
        result = {chr_name: 1 for chr_name in chr_names}
        total_breakage_cnt -= len(chr_names)
        while total_breakage_cnt > 0:
            total_breakage_cnt -= 1
            chr_name = list(numpy.random.choice(chr_names, size=1))[0]
            result[chr_name] += 1
        return result


class ChromothripsisRearrGenerator(KBreakGenerator):
    def __init__(self, rearr_type: Type[ChromothripsisRearrangement],
                 regions: List[Segment],
                 chr_cnt: Tuple[int, int], breakage_cnt: Tuple[int, int],
                 deletion_prob: float = 0.1, amplification_prob: float = 0.0,
                 id_prefix: str = "CHROMOTHRIPSIS",
                 *args, **kwargs):
        super().__init__(rearr_type, regions, chr_cnt, breakage_cnt, deletion_prob, amplification_prob, id_prefix, *args, **kwargs)

    @classmethod
    def get_random_telomere_matching(cls, chr_telomere_blocks: List[Tuple[str, Tuple[OrientedIndexedSegment, OrientedIndexedSegment]]]
                                     ) -> List[Tuple[OrientedIndexedSegment, OrientedIndexedSegment]]:
        result = []
        for chr_name, (left, right) in chr_telomere_blocks:
            result.append((left, right))
        return result

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "ChromothripsisRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "CHROMOTHRIPSIS")
        regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        chr_cnt_range = cls.parse_range(str(yaml_dict.get(CHR_CNT, "1-2")))
        breakage_cnt_range = cls.parse_range(str(yaml_dict.get(BREAKAGE_CNT, "20-100")))
        del_prob = float(yaml_dict.get(DELETION_PROBABILITY, "0.1"))
        ampl_prob = float(yaml_dict.get(AMPLIFICATION_PROBABILITY, "0.0"))
        return cls(ChromothripsisRearrangement, regions, chr_cnt_range, breakage_cnt_range, del_prob, ampl_prob, id_prefix)


class ChromoplexyRearrGenerator(KBreakGenerator):
    def __init__(self, rearr_type: Type[ChromoplexyRearrangement],
                 regions: List[Segment],
                 chr_cnt: Tuple[int, int], breakage_cnt: Tuple[int, int],
                 deletion_prob: float = 0.1, amplification_prob: float = 0.0,
                 id_prefix: str = "CHROMOPLEXY",
                 *args, **kwargs):
        super().__init__(rearr_type, regions, chr_cnt, breakage_cnt, deletion_prob, amplification_prob, id_prefix, *args, **kwargs)

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict) -> "ChromoplexyRearrGenerator":
        id_prefix = yaml_dict.get(ID_PREFIX, "CHROMOPLEXY")
        regions = cls.parse_regions(yaml_dict.get(REGIONS, ""))
        chr_cnt_range = cls.parse_range(str(yaml_dict.get(CHR_CNT, "5-10")))
        breakage_cnt_range = cls.parse_range(str(yaml_dict.get(BREAKAGE_CNT, "5-10")))
        del_prob = float(yaml_dict.get(DELETION_PROBABILITY, "0.0"))
        ampl_prob = float(yaml_dict.get(AMPLIFICATION_PROBABILITY, "0.0"))
        return cls(ChromoplexyRearrangement, regions, chr_cnt_range, breakage_cnt_range, del_prob, ampl_prob, id_prefix)


class RearrangementGeneratorGroup(object):
    str_rearr_gen_type_mapping = {
        "ins": InsRearrangementGenerator,
        "del": DelRearrangementGenerator,
        "tandemdup": TandemDupRearrangementGenerator,
        "tradup": DupRearr,
        "inv": InvRearrangementGenerator,
        "tra": TraRearrangementGenerator,
        "chrdup": MultiChrDupRearrGenerator,
        "chrdel": MultiChrDelRearrGenerator,
        "wgd": WGDRearrGenerator,
        "bfb": BFBRRearrGenerator,
        "chromoplexy": ChromoplexyRearrGenerator,
        "chromothripsis": ChromothripsisRearrGenerator,
        "kbreak": KBreakGenerator,
    }

    def __init__(self, cnt: int, rearrangements: List[RearrangementGenerator], probabilities: List[float]):
        self.cnt = cnt
        self.generated_cnt = 0
        self.rearr_generators: List[RearrangementGenerator] = rearrangements
        self.probabilities: List[float] = probabilities
        if sum(self.probabilities) != 1.0:
            raise ValueError(f"Attempted to create a rearrangement group, but combined probabilities don't add up to 1, but rather to {sum(self.probabilities)}")
        if len(self.rearr_generators) != len(self.probabilities):
            raise ValueError(f"Attempted to create a rearrangement group, but the number {len(self.rearr_generators)} of possible rearrangements "
                             f"does not equal to the number {len(self.probabilities)} probabilities")

    def pull_next_rearr_generator(self) -> RearrangementGenerator:
        return list(numpy.random.choice(self.rearr_generators, size=1, p=self.probabilities))[0]

    def get_next_random_rearrangement(self, genome: Genome, forbidden_ids: Optional[Set[str]] = None,
                                      constraints: Optional[Set[Type[Constraint]]] = None) -> Rearrangement:
        if forbidden_ids is None:
            forbidden_ids = set()
        rearr_generator = self.pull_next_rearr_generator()
        self.generated_cnt += 1
        return rearr_generator.generate_rearrangement(genome=genome, forbidden_ids=forbidden_ids, constraints=constraints)

    @classmethod
    def from_yaml_dict(cls, yaml_dict: dict, logger: Optional[logging.Logger] = None) -> "RearrangementGeneratorGroup":
        logger = logger or logging.getLogger('dummy')
        rearr_generators: List[Tuple[RearrangementGenerator, Optional[float]]] = []
        probability_left = 1.0
        cnt_range = RearrangementGenerator.parse_range(str(yaml_dict.get(CNT, 1)))
        cnt = numpy.random.randint(cnt_range[0], cnt_range[1] + 1)
        logger.info(f"Number of rearrangements for the generator = {cnt} sampled from range {cnt_range[0]}-{cnt_range[1]}")
        if REARRANGEMENT_STRING in yaml_dict:
            logger.info(f"Generators = ExactRearrangementGenerator (1.0) as there is a {REARRANGEMENT_STRING} entry")
            return cls(1, [ExactRearrangementGenerator(Rearrangement, rearr_string=yaml_dict[REARRANGEMENT_STRING])], [1.0])
        elif ENTRIES in yaml_dict:
            for entry in yaml_dict[ENTRIES]:
                rearr_generator = cls.create_rearr_generator_from_yaml_dict(entry)
                probability = float(entry[PROBABILITY]) if PROBABILITY in entry else None
                rearr_generators.append((rearr_generator, probability))
            assigned_probability = [r for r in rearr_generators if r[1] is not None]
            unassigned_probability = [r for r in rearr_generators if r[1] is None]
            used_probability = sum(r[1] for r in assigned_probability)
            probability_left -= used_probability
            if probability_left < 0:
                raise ValueError(f"Can not create a rearrangement generators group as explicitly specified probabilities add up to > 1")
            if len(unassigned_probability) > 0:
                unassigned_individual_probability = probability_left / len(unassigned_probability)
            else:
                unassigned_individual_probability = 0.0
            result_rearr_generators = []
            result_probabilities = []
            for rg, prob in rearr_generators:
                result_rearr_generators.append(rg)
                result_probabilities.append(prob if prob is not None else unassigned_individual_probability)
            assert sum(result_probabilities) == 1.0
            gen_prob_strings = []
            for rearr_g, rearr_g_prob in zip(result_rearr_generators, result_probabilities):
                gen_prob_strings.append(f"{type(rearr_g).__name__} ({rearr_g_prob:0.2f})")
            logger.info(f"Generators = {', '.join(gen_prob_strings)}")
            return cls(cnt, result_rearr_generators, result_probabilities)
        else:
            raise ValueError(f"Can not parse a rearrangement generators group. No `{ENTRIES}` entry and no '{REARRANGEMENT_STRING}' entry")

    @classmethod
    def create_rearr_generator_from_yaml_dict(cls, yaml_dict: dict) -> RearrangementGenerator:
        rearr_gen_type = yaml_dict.get(TYPE, None)
        if rearr_gen_type is None:
            raise ValueError(f"Rearrangement generator type not specified")
        if rearr_gen_type not in cls.str_rearr_gen_type_mapping:
            raise ValueError(f"Unknown rearrangement generator type '{rearr_gen_type}")
        return cls.str_rearr_gen_type_mapping[rearr_gen_type].from_yaml_dict(yaml_dict)
