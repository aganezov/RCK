import bisect
import itertools
import logging
from collections import defaultdict
from copy import deepcopy
from enum import Enum
from io import StringIO
from typing import Set, Union, Iterable, Tuple, List, Optional, Iterator, Type, Dict

import Bio
from Bio import Phylo

import networkx as nx

from rck.core.structures import Position, Segment, Strand, reverse_segment, ChromosomeType, Chromosome, Genome
from rck.utils.sim.rearrangement import Rearrangement, RearrangementChrTarget, OrientedIndexedSegment, RearrangementChrResult
from rck.utils.sim.rearr_generator import Constraint, ForbiddenTelomereCreationConstraint, HaploidISConstraint, DiploidISConstraint, RearrangementGeneratorGroup, \
    REARRANGEMENT_STRING

STARTING_GENOME = "starting_genome"
INFINITE_SITES = "infinite_sites"
REARRANGEMENTS = "rearrangements"
BRANCH = "branch"
TREE = "tree"
CONSTRAINTS = "constraints"
ROOT = "root"
EVOLUTION = "evolution"
REFERENCE = "ref"
REFERENCE_FILE = "ref_file"


class InfiniteSitesConstrain(Enum):
    DIPLOID = "dIS"
    HAPLOID = "hIS"

    def __str__(self) -> str:
        return self.value

    @classmethod
    def from_str(cls, string: str) -> "InfiniteSitesConstrain":
        if string not in ["dIS", "hIS"]:
            raise ValueError(f"Attempted to create Infinite Sites Constraint from string {string}. Must be 'dIS' or 'hIS'")
        if string == "dIS":
            return cls.DIPLOID
        return cls.HAPLOID


class RearrManager(object):

    @classmethod
    def rearrange_genome(cls, rearrangement: Rearrangement, target: Genome, isc: Union[InfiniteSitesConstrain, None] = InfiniteSitesConstrain.DIPLOID,
                         rearr_genome_chr_prefix: str = "") -> Genome:
        """
        Does not affect the target genome, creates a rearranged one, if rearrangement is possible.
        """
        target_genome = target
        #
        # pre application check can be done without applying the rearrangement
        #
        unsuitable_chrs = cls.unsuitable_rearr_chr_for_genome(target_genome, rearrangement)
        if len(unsuitable_chrs) > 0:
            raise ValueError(f"Attempted to apply rearrangement {str(rearrangement)} to a genome. "
                             f"Rearrangement target chromosomes are not a subset of genome's chromosomes. "
                             f"Unsuitable chromosomes {str(unsuitable_chrs)}.")
        for rearr_chr_target in rearrangement.target.rearr_chr_targets:
            target_chromosome = target_genome.chromosomes[rearr_chr_target.chr_name]
            unsuitable_breakage_positions = cls.unsuitable_rearr_loc_for_chromosome(target_chromosome, rearr_chr_target, allow_adj_positions=isc is not None)
            if len(unsuitable_breakage_positions) > 0:
                raise ValueError(f"Attempted to apply rearrangement {rearrangement} to a genome. "
                                 f"For target chromosome {rearr_chr_target.chr_name} breakage positions are not suitable. "
                                 f"Unsuitable chromosomes {str(unsuitable_breakage_positions)}.")
        #
        # breaking target genome
        #
        broken_chromosomes: List[Iterable[Chromosome]] = []
        all_breakage_positions = set()
        target_chr_names = set()
        for chr_target in rearrangement.target.rearr_chr_targets:
            target_chr_names.add(chr_target.chr_name)
            breakage_positions: Iterable[Position] = cls.get_breakage_positions(target=target_genome.chromosomes[chr_target.chr_name],
                                                                                locations=chr_target.breakage_coordinates)
            all_breakage_positions.update(breakage_positions)
            broken_chromosomes.append(cls.break_chromosome(target=target_genome.chromosomes[chr_target.chr_name], breakage_positions=breakage_positions,
                                                           locations=chr_target.breakage_coordinates))
        broken_chromosomes: List[Chromosome] = list(itertools.chain(*broken_chromosomes))
        #
        # check that the result of the rearrangement can be constructed from the broken parts
        #
        unsuitable_rearr_result_indexes = cls.unsuitable_rearr_result_indexes(broken_chromosomes=broken_chromosomes, rearrangement=rearrangement)
        if len(unsuitable_rearr_result_indexes) > 0:
            raise ValueError(f"Attempted to apply rearrangement {rearrangement} to a genome. "
                             f"Following results oriented indexes {','.join(map(str, unsuitable_rearr_result_indexes))} are not suitable for current rearrangement, "
                             f"as integer indexes in them are too large (i.e., refer to broken parts numbers that are not produced by this rearrangement.")
        #
        # arranging and assembling new chromosomes in place the ones that were targeted and broken
        #
        rearranged_chromosomes = []
        for rearr_chr_result in rearrangement.result.chr_results:
            rearranged_chromosomes.append(cls.arrange_result_chromosome(broken_chromosomes=broken_chromosomes, rearr_chr_result=rearr_chr_result))

        #
        # breaking non-target chromosomes with obtained breakage locations and then reassembling them with the same order and orientation
        #   (i.e., segment fragmentation without any structural rearrangements)
        #
        unaffected_refined_chromosomes = []
        for chr_name, chromosome in target_genome.chromosomes.items():
            if chr_name in target_chr_names:
                continue
            unaffected_refined_chromosomes.append(cls.refine_chromosome(target=chromosome, breakage_positions=all_breakage_positions))

        #
        # refining targeted chromosomes, as breakage positions may affect multiple segments,
        #   but only breakages that correspond to absolute coordinate locations were broken in the previous step
        #
        rearranged_chromosomes = [cls.refine_chromosome(target=chromosome, breakage_positions=all_breakage_positions) for chromosome in rearranged_chromosomes]

        #
        # constructing the overall rearranged genome
        #
        result = Genome(chromosomes=[])
        for cnt, chromosome in enumerate(unaffected_refined_chromosomes + rearranged_chromosomes):
            chromosome.name = cls.constr_new_chr_name(step=rearr_genome_chr_prefix, chr_name=cnt)
            result.chromosomes[chromosome.name] = chromosome

        #
        # recording the genome that underwent rearrangement into the history
        #
        return result

    @classmethod
    def constr_new_chr_name(cls, step: Union[str, int], chr_name: Union[str, int], separator: str = "_") -> str:
        return str(step) + separator + str(chr_name)

    @classmethod
    def arrange_result_chromosome(cls, broken_chromosomes: List[Chromosome], rearr_chr_result: RearrangementChrResult) -> Chromosome:
        result_chr_fragments: List[Chromosome] = []
        for oi in rearr_chr_result.oriented_segments:
            if isinstance(oi, OrientedIndexedSegment):
                fragment = deepcopy(broken_chromosomes[oi.s_index])
                if oi.orientation == Strand.REVERSE:
                    fragment.segments = cls.reverse_segments(fragment.segments)
                result_chr_fragments.append(fragment)
            else:
                result_chr_fragments.append(Chromosome(chr_name=oi.chromosome, segments=[oi], chr_type=ChromosomeType.LINEAR))
        result = cls.assemble_broken_chromosomes(chr_name="", chr_type=rearr_chr_result.chr_type, broken_chromosomes=result_chr_fragments)
        return result

    @classmethod
    def refine_genome(cls, target: Genome, breakage_positions: Iterable[Position]) -> Genome:
        if not isinstance(breakage_positions, list):
            breakage_positions = list(breakage_positions)
        result = Genome(chromosomes=[])
        for chromosome in target.chromosomes.values():
            refined_chromosome = cls.refine_chromosome(target=chromosome, breakage_positions=breakage_positions)
            result.chromosomes[refined_chromosome.name] = refined_chromosome
        return result

    @classmethod
    def refine_chromosome(cls, target: Chromosome, breakage_positions: Iterable[Position]) -> Chromosome:
        broken_chromosomes = cls.break_chromosome(target=target, breakage_positions=breakage_positions)
        return cls.assemble_broken_chromosomes(chr_name=target.name, chr_type=target.chr_type, broken_chromosomes=broken_chromosomes)

    @classmethod
    def assemble_broken_chromosomes(cls, chr_name: str, chr_type: ChromosomeType, broken_chromosomes: Iterable[Chromosome]) -> Chromosome:
        """
        Does not alter the broken chromosomes. References to segments input chromosomes are in the output.
        """
        result = Chromosome(chr_name=chr_name, chr_type=chr_type, segments=list(itertools.chain(*[chromosome.segments for chromosome in broken_chromosomes])))
        return result

    @classmethod
    def unsuitable_rearr_chr_for_genome(cls, genome: Genome, rearrangement: Rearrangement) -> Set[str]:
        """
        Does not alter the genome
        """
        available_chromosomes = set(genome.chromosomes.keys())
        target_chromosomes = set([rct.chr_name for rct in rearrangement.target.rearr_chr_targets])
        return target_chromosomes - available_chromosomes

    @classmethod
    def unsuitable_rearr_result_indexes(cls, broken_chromosomes: List[Chromosome], rearrangement: Rearrangement) -> List[OrientedIndexedSegment]:
        """
        Does not alter the broken chromosomes list
        """
        result = []
        for result_chr in rearrangement.result.chr_results:
            for oi in result_chr.oriented_segments:
                if isinstance(oi, OrientedIndexedSegment) and oi.s_index >= len(broken_chromosomes):
                    result.append(oi)
        return result

    @classmethod
    def unsuitable_rearr_loc_for_chromosome(cls, chromosome: Chromosome, rearr_chr_target: RearrangementChrTarget, allow_adj_positions: bool = False) -> Set[int]:
        """
        Does not alter the chromosome
        """
        chr_length = sum(s.length for s in chromosome.segments)
        bad_coordinates = set()
        for coordinate in set(rearr_chr_target.breakage_coordinates):
            if coordinate < 0 or coordinate > chr_length:
                bad_coordinates.add(coordinate)
        if not allow_adj_positions:
            current_coordinate = 0
            adj_positions = {current_coordinate}
            for segment in chromosome.segments:
                current_coordinate += segment.length
                adj_positions.add(current_coordinate)
            adj_coordinates = set(rearr_chr_target.breakage_coordinates).intersection(adj_positions)
            bad_coordinates.update(adj_coordinates)
        return bad_coordinates

    @classmethod
    def get_breakage_positions(cls, target: Chromosome, locations: Iterable[int]) -> Iterable[Position]:
        """
        Does not alter target chromosomes
        """
        result = set()
        if target.chr_type == ChromosomeType.CIRCULAR:
            # for circular chromosomes breakage coordinates that are longed, than the chromosomes length must "wrap-around" via modular division
            #   coordinate 0 and chromosome.length are thus identical
            chr_length = target.length  # this is a property attribute computed on the fly, thus the memoization
            locations = [location % chr_length for location in locations]
        locations = sorted(set(locations))  # if multiple breakage location are pointing to the same underlying coordinate, they are treated as identical
        current_offset = 0  # current length offset from the beginning of the target chromosome
        segments = iter(target.segments)
        current_segment = next(segments, None)
        if len(locations) == 0:
            return result
        if locations[0] == 0:
            result.add(current_segment.start_position.get_non_hap_copy())
            locations = locations[1:]
        locations = iter(locations)
        current_location = next(locations, None)
        while current_location is not None and current_segment is not None:
            if current_offset < current_location <= current_offset + current_segment.length:
                result.add(cls.get_breakage_position(current_segment, current_location - current_offset))
                current_location = next(locations, None)
            else:
                current_offset += current_segment.length
                current_segment = next(segments, None)
        if current_location is not None and current_location > current_offset:
            # can only happen for linear chromosomes, as for circular locations are % against the total chromosomal length
            result.add(target.segments[-1].end_position.get_non_hap_copy())
        return result

    @classmethod
    def get_breakage_position(cls, segment: Segment, relative_location: int) -> Position:
        """
        Does not alter segment
        """
        if relative_location > segment.length:
            raise ValueError(f"Attempted to retrieve breakage position for segment {str(segment)} with relative location {relative_location}. "
                             f"Relative location must be less or equal to segment's length.")
        result = segment.start_position.get_non_hap_copy()
        if segment.is_reversed:
            result.coordinate -= (relative_location - 1)
            result.strand = Strand.REVERSE
        else:
            result.strand = Strand.FORWARD
            result.coordinate += (relative_location - 1)
        return result

    @classmethod
    def break_chromosome(cls, target: Chromosome, breakage_positions: Iterable[Position], locations: Optional[List[int]] = None) -> Iterable[Chromosome]:
        """
        Does not alter the target chromosome, nor segments in target chromosome, nor references to segments are part of the result
        """
        result = []
        positions_by_chr = defaultdict(list)
        for p in breakage_positions:
            positions_by_chr[p.chromosome].append(p)
        segments = iter(target.segments)
        current_segment = deepcopy(next(segments, None))
        current_new_chr_segments = []
        current_offset = 0
        while current_segment is not None:
            related_positions = sorted(positions_by_chr[current_segment.chromosome], key=lambda e: e.coordinate)
            # positions are sorted with by coordinates and strand taken into account, with strands order (for same coordinate) determined by segment orientation
            positions_in_current_segment = cls.positions_in_segment(current_segment, related_positions)
            at_the_end_of_segment = False
            for p in positions_in_current_segment:
                assert not at_the_end_of_segment  # sanity check
                if locations is not None and current_offset + abs(p.coordinate - current_segment.start_coordinate) + 1 not in locations:
                    continue
                broken_segments = cls.break_segment(current_segment, p)  # produces a pair of segments, with possible None values (should not have two None values)
                assert broken_segments != (None, None)
                if broken_segments[0] is not None:  # if broken_segment[0] is None: breakage was by the current segment's left extremity, no new "left" segment is produced
                    current_new_chr_segments.append(broken_segments[0])
                    current_offset += broken_segments[0].length
                result.append(current_new_chr_segments)  # once the breakage happens, new sub-chromosome is produced
                current_new_chr_segments = []  # starting new sub-chromosome
                if broken_segments[1] is not None:  # if broken_segment[1] is None, we have broken segment by its right-most extremity
                    current_segment = broken_segments[1]
                else:
                    # ended the current segment fragmentation, shall stop iterations over positions at the top, flag below work with check assertion in the positions iter loop
                    # current_segment = deepcopy(next(segments, None))
                    at_the_end_of_segment = True
            if current_segment is not None:
                # current_segment here can be non only if we've exhausted the iteration over all segment in the target chromosome
                if not at_the_end_of_segment:
                    current_new_chr_segments.append(current_segment)
                    current_offset += current_segment.length
                current_segment = deepcopy(next(segments, None))
        result.append(current_new_chr_segments)
        if target.chr_type == ChromosomeType.CIRCULAR:
            # for circular chromosome, left-most sub-chromosome, and right-most sub-chromosome must be merged together
            #   when breakage happens at the extremities of the linear version of the circular chromosome, will create empty sub-chromosome lists
            #   and such merge will correctly merge empty list(s) with (non)empty lists, correctly producing linear sub-chromosomes
            result[0] = cls.reverse_segments(result[-1]) + result[0]
        result = [Chromosome(f"_tmp_{cnt}", n_segments, chr_type=ChromosomeType.LINEAR) for cnt, n_segments in enumerate(result) if len(n_segments) > 0]
        return result

    @classmethod
    def positions_in_segment(cls, segment: Segment, positions: Iterable[Position]) -> Iterable[Position]:
        positions = list(positions)
        positions = sorted(positions, key=lambda p: (p.coordinate, 0 if p.strand == Strand.REVERSE else 1))
        left_index = bisect.bisect_left([p.coordinate for p in positions], min(segment.start_coordinate, segment.end_coordinate))
        right_index = bisect.bisect_right([p.coordinate for p in positions], max(segment.start_coordinate, segment.end_coordinate))
        if 0 < len(positions) == left_index and max(segment.start_coordinate, segment.end_coordinate) == max([p.coordinate for p in positions]):
            left_index -= 1
        positions = positions[left_index:right_index]
        return sorted(positions, key=lambda p: (p.coordinate, 0 if p.strand == Strand.REVERSE else 1), reverse=segment.is_reversed)

    @classmethod
    def break_segment(cls, segment: Segment, position: Position) -> Tuple[Optional[Segment], Optional[Segment]]:
        """
        Does not alter the segment
        """
        if segment.chromosome != position.chromosome:
            raise ValueError(f"Attempting to break a segment {str(segment)} with position {str(position)} and they disagree on chromosome of origin.")
        if position.coordinate < min(segment.start_coordinate, segment.end_coordinate) or position.coordinate > max(segment.start_coordinate, segment.end_coordinate):
            raise ValueError(f"Attempting to break segment {str(segment)} with position {str(position)} and position coordinate lies outside of segment.")
        if position.non_hap_eq(segment.start_position):
            return None, deepcopy(segment)
        elif position.non_hap_eq(segment.end_position):
            return deepcopy(segment), None
        nl_segment = deepcopy(segment)
        nr_segment = deepcopy(segment)
        nl_extra_coord_padding = 0
        nr_extra_coord_padding = 0
        if position.strand == Strand.FORWARD:
            if segment.is_reversed:
                nl_extra_coord_padding = 1
            else:
                nr_extra_coord_padding = 1
        else:
            if segment.is_reversed:
                nr_extra_coord_padding = -1
            else:
                nl_extra_coord_padding = -1
        nl_segment.end_position.coordinate = position.coordinate + nl_extra_coord_padding
        nr_segment.start_position.coordinate = position.coordinate + nr_extra_coord_padding
        return nl_segment, nr_segment

    @classmethod
    def reverse_segments(cls, segments: Union[List[Optional[Segment]], Tuple[Optional[Segment], ...]]) -> List[Optional[Segment]]:
        """
        Does not alter the segments iterable
        """
        return [reverse_segment(s) if s is not None else s for s in reversed(segments)]


class SimBranchManager(object):
    def __init__(self, starting_genome: Optional[Genome], rearr_history: List[RearrangementGeneratorGroup],
                 branch: str = "R-A", external_forbidden_ids: Optional[Set[str]] = None,
                 constraints: Optional[Set[Type[Constraint]]] = None,
                 logger: Optional[logging.Logger] = None):
        self.forbidden_ids = external_forbidden_ids if external_forbidden_ids is not None else set()
        self.constraints = constraints if constraints is not None else self.default_constraints()
        self.branch = branch
        self.starting_genome: Genome = starting_genome
        self.current_genome: Genome = self.starting_genome
        self.input_history: List[RearrangementGeneratorGroup] = rearr_history
        self.applied_history: List[Rearrangement] = []
        self.rearranged_history: List[Genome] = []
        self.current_rearr_generator_group = self.input_history[0]
        self.current_step = 0
        self.logger = logger or logging.getLogger('dummy')

    @classmethod
    def default_constraints(cls):
        return {ForbiddenTelomereCreationConstraint}

    @classmethod
    def from_yaml_dict(cls, starting_genome: Optional[Genome], yaml_dict: dict, logger: Optional[logging.Logger] = None) -> "SimBranchManager":
        logger = logger or logging.getLogger('dummy')
        rearrangements_entries = yaml_dict.get(REARRANGEMENTS, [])
        branch = yaml_dict.get(BRANCH, "R-A")
        logger.info(f"Generating simulation manager for the branch {branch}")
        if len(rearrangements_entries) == 0:
            raise ValueError(f"Can not create a branch simulator manager, as supplied list of rearrangements is empty or does not exist")
        rearr_history = [RearrangementGeneratorGroup.from_yaml_dict(e, logger=logger) for e in rearrangements_entries]
        return cls(starting_genome, rearr_history, branch, logger=logger)

    def generate_next_rearrangement(self) -> Iterator[Rearrangement]:
        for rearr_gen_group in self.input_history:
            for _ in range(rearr_gen_group.cnt):
                self.current_step += 1
                yield rearr_gen_group.get_next_random_rearrangement(self.current_genome, self.forbidden_ids, self.constraints)

    def apply_rearrangement(self, rearrangement: Rearrangement, store_rearr_history: bool = True, store_intermediate_genomes: bool = False):
        if store_rearr_history:
            self.applied_history.append(rearrangement)
        rearranged_genome = RearrManager.rearrange_genome(rearrangement, self.current_genome, rearr_genome_chr_prefix=str(self.current_step))
        if store_intermediate_genomes:
            self.rearranged_history.append(self.current_genome)
        self.forbidden_ids.add(rearrangement.self_id)
        self.current_genome = rearranged_genome

    def simulate(self, store_rearr_history: bool = True, store_intermediate_genomes: bool = False):
        if self.starting_genome is None:
            raise ValueError(f"Can not simulate evolution along branch '{self.branch}' because the starting genome is not set")
        for rearr in self.generate_next_rearrangement():
            self.logger.info(f"Applying generated rearrangement {str(rearr)}")
            self.apply_rearrangement(rearr, store_rearr_history=store_rearr_history, store_intermediate_genomes=store_intermediate_genomes)

    def applied_rearr_history_as_yaml(self) -> List[dict]:
        return [{REARRANGEMENT_STRING: str(entry)} for entry in self.applied_history]


class SimManager(object):
    CONSTR_STR_TO_TYPE = {
        "NoNewTelomeres": ForbiddenTelomereCreationConstraint,
        "DiploidIS": DiploidISConstraint,
        "HaploidIS": HaploidISConstraint,
    }

    def __init__(self,
                 reference: Genome,
                 tree: nx.DiGraph,
                 branches_setup: Dict[str, SimBranchManager],
                 root_name: str = "R",
                 tree_constraints: Optional[Set[Type[CONSTRAINTS]]] = None,
                 biopython_phylo_tree: Optional = None,
                 logger: Optional[logging.Logger] = None,
                 ):
        self.reference = reference
        self.tree_constraints: Optional[Set[Type[CONSTRAINTS]]] = tree_constraints
        self.forbidden_ids: Set[str] = set()
        self.bio_python_phylo_tree = biopython_phylo_tree
        self.tree: nx.DiGraph = tree
        self.root_name: str = root_name
        self.branches_setup: Dict[str, SimBranchManager] = branches_setup
        self.genomes: Dict[str, Genome] = {self.root_name: deepcopy(reference)}
        self.refined_genomes: Dict[str, Genome] = {}
        self.logger = logger or logging.getLogger('dummy')

    @classmethod
    def parse_newick_str_into_biopython_phylo(cls, tree: str, root_name: str = "R"):
        tree = Phylo.read(StringIO(tree), format="newick")
        tree.root_with_outgroup({"name": root_name})
        tree.rooted = True
        return tree

    @classmethod
    def parse_tree(cls, tree: str, root_name: str = "R"):
        tree = cls.parse_newick_str_into_biopython_phylo(tree, root_name)
        tree = Phylo.to_networkx(tree)
        return tree

    @classmethod
    def constraint_from_yaml_list(cls, yaml_list: list) -> Set[Type[Constraint]]:
        result = set()
        for e in yaml_list:
            if e not in cls.CONSTR_STR_TO_TYPE:
                raise ValueError(f"Unsupported constraint str value '{e}'")
            result.add(cls.CONSTR_STR_TO_TYPE[e])
        return result

    @classmethod
    def get_branches_setup_from_setup_list(cls, setup_list: List, logger: Optional[logging.Logger] = None) -> Dict[str, SimBranchManager]:
        logger = logger or logging.getLogger('dummy')
        result: Dict[str, SimBranchManager] = {}
        for entry in setup_list:
            branch = entry.get(BRANCH, "R-A")
            result[branch] = SimBranchManager.from_yaml_dict(starting_genome=None, yaml_dict=entry, logger=logger)
        return result

    @classmethod
    def phylogeny_sources(cls, tree: nx.DiGraph) -> Set:
        result = {node.name for node in tree.nodes}
        for node in tree:
            for _ in tree.predecessors(node):
                result.remove(node.name)
                break
        return result

    @classmethod
    def validate_tree(cls, tree, root_name: str):
        if not isinstance(tree, nx.DiGraph):
            raise ValueError(f"Phylogeny tree failed to convert to the networkx.DiGraph object. This is mandatory")
        if not nx.is_directed_acyclic_graph(tree):
            raise ValueError(f"Supplied phylogeny tree is not a DAG. This is mandatory")
        tree_roots = cls.phylogeny_sources(tree)
        if len(tree_roots) != 1:
            raise ValueError(f"The number of roots in supplied phylogeny does not equal to 1. Should not be happening")
        inferred_root = tree_roots.pop()
        if inferred_root != root_name:
            ValueError(f"Phylogeny tree root {inferred_root} does not equal specified/default root {root_name}")

    @classmethod
    def from_yaml_dict(cls, yaml_dict: Dict, reference: Optional[Genome] = None, logger: Optional[logging.Logger] = None) -> "SimManager":
        logger = logger or logging.getLogger('dummy')
        constraints = cls.constraint_from_yaml_list(yaml_dict.get(CONSTRAINTS, []))
        root_name = yaml_dict.get(ROOT, "R")
        logger.info(f"Phylogeny root set to {root_name}")
        tree = cls.parse_tree(tree=yaml_dict.get(TREE, "R(A);"), root_name=root_name)
        biopython_tree = cls.parse_newick_str_into_biopython_phylo(tree=yaml_dict.get(TREE, "R(A);"), root_name=root_name)
        cls.validate_tree(tree, root_name)
        logger.info(f"Validated phylogeny tree successfully")
        if reference is None:
            reference_str = yaml_dict.get(REFERENCE, None)
            if reference_str is not None:
                reference = Genome.from_strs(reference_str.split("\n"))
            else:
                ref_file = yaml_dict.get(REFERENCE_FILE, None)
                if ref_file is not None:
                    lines = open(ref_file, "rt").readlines()
                    reference = Genome.from_strs(lines)
                else:
                    raise ValueError(f"No reference was set. Neither directly, nor through the 'ref' or 'ref_file' fields in the settings")
        logger.info(f"Reference genome for root node in the phylogeny was set")
        evolution: List[Dict] = yaml_dict[EVOLUTION]
        logger.info(f"Generating phylogeny branch-specific individual evolution managers")
        branches_setup = cls.get_branches_setup_from_setup_list(setup_list=evolution, logger=logger)
        return cls(reference, tree, branches_setup, root_name, constraints, biopython_tree, logger=logger)

    @classmethod
    def get_tree_node_by_name(cls, tree: nx.DiGraph, node_name: str):
        for node in tree:
            if node.name == node_name:
                return node
        raise ValueError(f"Could not find a node with name {node_name} in phylo tree")

    def simulate(self, refine_after: bool = True, store_intermediate_genomes: bool = False):
        for branch in nx.dfs_edges(self.tree, source=self.get_tree_node_by_name(self.tree, self.root_name)):
            branch_str = f"{branch[0].name}-{branch[1].name}"
            self.logger.info(f"Starting simulation for branch {branch_str}")
            self.simulate_branch(branch=branch_str, store_intermediate_genomes=store_intermediate_genomes)
            self.logger.info(f"Success. Completed simulation for branch {branch_str}")
        if refine_after:
            self.refine_genomes()

    def refine_genomes(self):
        breakage_positions = set()
        for genome in self.genomes.values():
            for position in genome.iter_positions():
                breakage_positions.add(position.get_non_hap_copy())
        for genome_name, genome in self.genomes.items():
            ref_genome = RearrManager.refine_genome(target=genome, breakage_positions=breakage_positions)
            self.refined_genomes[genome_name] = ref_genome

    @classmethod
    def flip_branch(cls, branch_str: str) -> str:
        return "-".join(reversed(branch_str.split("-")))

    def get_order_source_descendant(self, branch_str: str) -> Tuple[str, str]:
        node1_str, node2_str = branch_str.split("-")
        node1, node2 = self.get_tree_node_by_name(self.tree, node1_str), self.get_tree_node_by_name(self.tree, node2_str)
        if self.tree.has_edge(node1, node2):
            return node1_str, node2_str
        if self.tree.has_edge(node2, node1):
            return node2_str, node1_str
        raise ValueError(f"Can not determine which node from branch {branch_str} is a parental, as no such edge is present in the phylo tree")

    def simulate_branch(self, branch: str, store_rearrangements: bool = True, store_intermediate_genomes: bool = False):
        source_str, descendant_str = self.get_order_source_descendant(branch)
        source_node, descendant_node = self.get_tree_node_by_name(self.tree, source_str), self.get_tree_node_by_name(self.tree, descendant_str)
        branch_setup: Optional[SimBranchManager] = self.branches_setup.get(branch, self.branches_setup.get(self.flip_branch(branch), None))
        if branch_setup is None:
            raise ValueError(f"Can not simulate evolution for branch {branch}: no evolutionary setup provided")
        if source_str not in self.genomes:
            raise ValueError(f"Can not simulate evolution for branch {branch}: not genome for source node {source_node} is available to start from")
        starting_genome = deepcopy(self.genomes[source_str])
        branch_setup.starting_genome = starting_genome
        branch_setup.current_genome = starting_genome
        branch_setup.simulate(store_rearr_history=store_rearrangements, store_intermediate_genomes=store_intermediate_genomes)
        self.genomes[descendant_str] = deepcopy(branch_setup.current_genome)

    def as_yaml(self) -> Dict:
        tree_str_io = StringIO()
        Phylo.write(self.bio_python_phylo_tree, tree_str_io, "newick", plain=True)
        tree_str_io.seek(0)
        result = {
            TREE: tree_str_io.readline().strip(),
            ROOT: self.root_name,
            EVOLUTION: []
        }
        for branch_name, branch_setup in self.branches_setup.items():
            branch_entry = {BRANCH: branch_name, REARRANGEMENTS: branch_setup.applied_rearr_history_as_yaml()}
            result[EVOLUTION].append(branch_entry)
        return result




