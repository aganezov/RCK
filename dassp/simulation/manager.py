from collections import Counter, defaultdict, namedtuple
from copy import deepcopy
import itertools
import numpy as np
import networkx as nx
import math
from functools import lru_cache

from dassp.core.structures import Adjacency, AdjacencyType, Haplotype, Phasing, HAPLOTYPE_PAIRS_TO_PHASING
from dassp.simulation.parts import ChromosomeGenerator, MutationType, Mutation, reverse_segments, duplicate_segments, delete_segments, HAPLOTYPE, \
    translocate_segments, ChrPart, split_and_reassemble_chromosomes
from core.structures import reverse_segment, propagate_haplotype_genome_to_positions
from simulation.parts import duplicate_chromosome, duplicate_genome

CHROMOSOMES_SIZE = 100
AB = True
PRESERVE_TELOMERES = "preserve_telomeres"
HIIS = "haplotype_independent_infinite_sites"
HSIS = "haplotype_specific_infinite_sites"

MUTATED_EXTREMITIES = "mutated_extremities"
REVERSAL_TELOMERE = "telomere_reversal"
DELETION_TELOMERE = "telomere_deletion"
DUPLICATION_TELOMERE = "telomere_duplication"
TRANSLOCATION_CHR1_TELOMERE = "telomere_chr1_translocation"
TRANSLOCATION_CHR2_TELOMERE = "telomere_chr2_translocation"
HOMOZYGOUS = "homozygous"
MUTATIONS_TYPES_SPECIFICATIONS = "mutations_types_specification"
MUTATIONS_TYPES = "mutations_types"
MUTATIONS_PROBABILITIES = "mutations_probabilities"
MUTATIONS_PROBABILITIES_BY_MUTATIONS_TYPES = "mutations_probabilities_by_mutations_types"
DUPLICATION_TANDEM_PROBABILITY = "tandem_duplication_probability"
DUPLICATION_REVERSED_PROBABILITY = "reversed_duplication_probability"

CHROMOSOME_INDEX = "chromosome_index"
CHROMOSOME_INDEXES = "chromosome_indexes"
CHROMOSOME_1_INDEX = "chromosome_1_index"
CHROMOSOME_2_INDEX = "chromosome_2_index"
REVERSAL_START_INDEX = "reversal_start_index"
REVERSAL_END_INDEX = "reversal_end_index"
REVERSAL_TRY_LIMIT = "reversal_try_limit"

POSITIONS = "positions"
HAPLOTYPES = "haplotypes"
LOCATIONS = "locations"

DUPLICATION_START_INDEX = "duplication_start_index"
DUPLICATION_END_INDEX = "duplication_end_index"
DUPLICATION_INSERT_SHIFT = "duplication_insert_shift"
DUPLICATION_INSERT_INDEX = "duplication_insert_index"
DUPLICATION_REVERSED = "duplication_reversed"
DUPLICATION_TRY_LIMIT = "duplication_try_limit"
DUPLICATION_ALLOW_DOWNSTREAM = "duplication_allow_downstream"

DELETION_START_INDEX = "deletion_start_index"
DELETION_END_INDEX = "deletion_end_index"
DELETION_TRY_LIMIT = "deletion_try_limit"

TRANSLOCATION_CHR1_INDEX = "translocation_chromosome_1_index"
TRANSLOCATION_CHR2_INDEX = "translocation_chromosome_2_index"
TRANSLOCATION_CC = "translocation_cc"
TRANSLOCATION_TRY_LIMIT = "translocation_try_limit"

GENOMES = "genomes"
GENOME = "genome"
MUTATIONS = "mutations"
MUTATION_CNT = "mutation_count"
HISTORY = "history"

CHROMOTHRIPSIS_MIN_CHR_CNT = "chromothripsis_min_chrs_cnt"
CHROMOTHRIPSIS_MAX_CHR_CNT = "chromothripsis_max_chrs_cnt"
CHROMOTHRIPSIS_MIN_BR_CNT = "chromothrpsis_min_breakage_cnt"
CHROMOTHRIPSIS_MAX_BR_CNT = "chromothripsis_max_breakage_cnt"
CHROMOTHRIPSIS_DUPLICATIONS = "chromothripsis_allow_duplications"
CHROMOTHRIPSIS_DELETIONS = "chromothripsis_allow_deletions"
CHROMOTHRIPSIS_TELOMERE = "telomere_chromothripsis"
CHROMOTHRIPSIS_TRY_LIMIT = "chromothripsis_try_limit"
CHROMOTHRIPSIS_DUPLICATIONS_PROB = "chromothripsis_duplications_probability"
CHROMOTHRIPSIS_DELETIONS_PROB = "chromothripsis_deletions_probability"
CHROMOTHRIPSIS_REVERSALS = "chromothripsis_allow_reversals"
CHROMOTHRIPSIS_REVERSALS_PROB = "chromothripsis_reversals_probability"
CHROMOTHRIPSIS_TELOMERE_PAIRING_PRESERVE = "chromothripsis_preserve_telomere_pairs"

CHROMOTHRIPSIS_PER_CHRS_SEGMENT_INDEXES = "chromothripsis_per_chrs_segment_indexes"
CHROMOTHRIPSIS_TARGET_STRUCTURE = "chromothripsis_target_structure"
MUT_CONFIG = "mutation_config"

SAVE_INTERMEDIATE_GENOMES = "save_intermediate_genomes"

MUTATION_CONFIG = {
    MUTATIONS_TYPES:
        [
            MutationType.REVERSAL,
            MutationType.DUPLICATION,
            MutationType.DELETION,
            MutationType.TRANSLOCATION,
            MutationType.CHROMOTHRIPSIS,
            MutationType.CHROMOSOME_DUP,
            MutationType.CHROMOSOME_DEL,
            MutationType.WGD,
        ],
    MUTATIONS_PROBABILITIES:
        [
            0.145,
            0.429,
            0.340,
            0.054,
            0.020,
            0.010,
            0.000,
            0.002,
        ],
    MUTATIONS_PROBABILITIES_BY_MUTATIONS_TYPES:
        {
            MutationType.REVERSAL: 0.145,
            MutationType.DUPLICATION: 0.429,
            MutationType.DELETION: 0.340,
            MutationType.TRANSLOCATION: 0.054,
            MutationType.CHROMOTHRIPSIS: 0.020,
            MutationType.CHROMOSOME_DUP: 0.010,
            MutationType.CHROMOSOME_DEL: 0.000,
            MutationType.WGD: 0.002,
        },
    PRESERVE_TELOMERES: True,
    HIIS: False,
    HSIS: True,
    MUTATIONS_TYPES_SPECIFICATIONS: {
        MutationType.DELETION: {
            DELETION_TELOMERE: False,
            DELETION_TRY_LIMIT: 100,
        },
        MutationType.REVERSAL: {
            REVERSAL_TELOMERE: False,
            REVERSAL_TRY_LIMIT: 100,
        },
        MutationType.DUPLICATION: {
            DUPLICATION_TELOMERE: False,
            DUPLICATION_TANDEM_PROBABILITY: 0.9,
            DUPLICATION_REVERSED_PROBABILITY: 0.25,
            DUPLICATION_TRY_LIMIT: 100,
            DUPLICATION_ALLOW_DOWNSTREAM: True,
        },
        MutationType.TRANSLOCATION: {
            TRANSLOCATION_CHR1_TELOMERE: False,
            TRANSLOCATION_CHR2_TELOMERE: False,
            HOMOZYGOUS: False,
            TRANSLOCATION_TRY_LIMIT: 100,
        },
        MutationType.CHROMOTHRIPSIS: {
            CHROMOTHRIPSIS_MIN_CHR_CNT: 2,
            CHROMOTHRIPSIS_MAX_CHR_CNT: 4,
            CHROMOTHRIPSIS_MIN_BR_CNT: 4,
            CHROMOTHRIPSIS_MAX_BR_CNT: 10,
            CHROMOTHRIPSIS_DUPLICATIONS: False,
            CHROMOTHRIPSIS_DELETIONS: False,
            CHROMOTHRIPSIS_TELOMERE: False,
            CHROMOTHRIPSIS_TRY_LIMIT: 100,
            CHROMOTHRIPSIS_DELETIONS_PROB: 0.2,
            CHROMOTHRIPSIS_DUPLICATIONS_PROB: 0.4,
            CHROMOTHRIPSIS_REVERSALS: True,
            CHROMOTHRIPSIS_REVERSALS_PROB: .5,
            CHROMOTHRIPSIS_TELOMERE_PAIRING_PRESERVE: False,
            HOMOZYGOUS: False,
        }
    }
}


class SimulationManager(object):
    def __init__(self, chrs_cnt, clones_cnt, chrs_size=CHROMOSOMES_SIZE, ab=AB, mutation_config=None):
        self.config = mutation_config if mutation_config is not None else MUTATION_CONFIG
        self.chromosomes_cnt = chrs_cnt
        self.clones_cnt = clones_cnt
        self.ab = ab
        self.chromosomes_size = chrs_size
        self.initial_genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=self.chromosomes_size,
                                                                       chromosomes_cnt=self.chromosomes_cnt,
                                                                       ab=self.ab)


class ClonalSubClonalSimulationManager(SimulationManager):
    def __init__(self, chrs_cnt,
                 clonal_mut_cnt, subclonal_mut_cnt,
                 chrs_size=CHROMOSOMES_SIZE, ab=AB):
        super(ClonalSubClonalSimulationManager, self).__init__(chrs_cnt=chrs_cnt,
                                                               chrs_size=chrs_size,
                                                               clones_cnt=2,
                                                               ab=ab)
        self.clonal_mutation_cnt = clonal_mut_cnt
        self.subclonal_mutation_cnt = subclonal_mut_cnt

    def generate_mutated_genomes(self):
        self.majority_clone_history = generate_mutated_genome(starting_genome=self.initial_genome, mutation_cnt=self.clonal_mutation_cnt, config=self.config)
        self.majority_clone = self.majority_clone_history[GENOMES][-1]
        self.minority_clone_history = generate_mutated_genome(starting_genome=self.majority_clone, mutation_cnt=self.subclonal_mutation_cnt, config=self.config)
        self.minority_clone = self.minority_clone_history[GENOMES][-1]


class Manager(object):
    def __init__(self, phylo_tree, config=None, assign_artificial_root=True):
        self.tree = phylo_tree
        self.root = next(nx.topological_sort(G=self.tree))
        if assign_artificial_root:
            root = "R"
            self.tree.add_edge(u=root, v=self.root)
            self.root = root
        self.config = config if config is not None else {}

    def generate_mutated_genomes(self, processed_nodes=None, processed_edges=None):
        if processed_nodes is None:
            processed_nodes = set()
        if processed_edges is None:
            processed_edges = set()
        processed_nodes.add(self.root)
        mutated_extremities = set()
        for edge in nx.bfs_edges(G=self.tree, source=self.root):
            u, v = edge
            mutation_config = self.tree.edges[u, v][MUT_CONFIG]
            save_inter = mutation_config.get(SAVE_INTERMEDIATE_GENOMES, None)
            if save_inter is None:
                save_inter = self.config.get(SAVE_INTERMEDIATE_GENOMES, False)
            mut_cnt = mutation_config.get(MUTATION_CNT, None)
            if mut_cnt is None:
                mut_cnt = self.tree.edges[u, v][MUTATION_CNT]
            starting_genome = self.tree.nodes[u][GENOME]
            mutation_config[MUTATED_EXTREMITIES] = mutated_extremities
            history = generate_mutated_genome(starting_genome=starting_genome,
                                              mutation_cnt=mut_cnt,
                                              config=mutation_config,
                                              save_intermediate_genome=save_inter)
            self.tree.edges[u, v][HISTORY] = history
            processed_edges.add(edge)
            self.tree.nodes[v][GENOME] = history[GENOMES][-1]
            processed_nodes.add(v)

    def propagate_haplotypes_to_positions(self):
        for node in self.tree:
            genome = self.tree.nodes[node][GENOME]
            propagate_haplotype_genome_to_positions(genome=genome, inplace=True)


def assign_uniform_mutation_configs(manager, mut_cnt, mut_config):
    for edge in nx.bfs_edges(G=manager.tree, source=manager.root):
        u, v = edge
        manager.tree.edges[u, v][MUT_CONFIG] = deepcopy(mut_config)
        manager.tree.edges[u, v][MUTATION_CNT] = mut_cnt


def assign_reference_genome(manager, ref_genome, copy=True):
    manager.tree.nodes[manager.root][GENOME] = ref_genome if not copy else deepcopy(ref_genome)


def generate_mutation(mutation_type, genome, config):
    if mutation_type not in MUTATION_TYPES_to_GENERATING_FUNCTIONS:
        raise ValueError("Unsupported mutation type {mt}".format(mt=repr(mutation_type)))
    generating_function = MUTATION_TYPES_to_GENERATING_FUNCTIONS[mutation_type]
    return generating_function(genome=genome, config=config)


BreakageExtremity = namedtuple('BreakageExtremity', ['p', 'h'])
BreakageLocation = namedtuple('BreakageLocation', ['be1', 'be2', 'bi'])


def breakage_location_is_on_telomere(breakage_location):
    return breakage_location_is_on_left_telomere(breakage_location=breakage_location) or breakage_location_is_on_right_telomere(breakage_location=breakage_location)


def breakage_location_is_on_left_telomere(breakage_location):
    return breakage_location.be1.p is None


def breakage_location_is_on_right_telomere(breakage_location):
    return breakage_location.be2.p is None


def suitable_for_breakage(breakage_location, config):
    bl = breakage_location
    if config[HIIS]:
        if bl.be1 is not None and bl.be1.p in config[MUTATED_EXTREMITIES]:
            return False
        if bl.be2 is not None and bl.be2.p in config[MUTATED_EXTREMITIES]:
            return False
    elif config[HSIS]:
        if bl.be1 is not None and (bl.be1.p, bl.be1.h) in config[MUTATED_EXTREMITIES]:
            return False
        if bl.be2 is not None and (bl.be2.p, bl.be2.h) in config[MUTATED_EXTREMITIES]:
            return False
    if config[PRESERVE_TELOMERES] and breakage_location_is_on_telomere(breakage_location=breakage_location):
        return False
    return True


def breakage_location_hashable_representation(breakage_location, config):
    bl = breakage_location
    if bl.be1.p is None and bl.be2.p is None:
        return None, None
    if config[HIIS]:
        if bl.be1.p is None:
            return None, bl.be2.p
        elif bl.be2.p is None:
            return bl.be1.p, None
        else:
            return (bl.be1.p, bl.be2.p) if bl.be1.p < bl.be2.p else (bl.be2.p, bl.be1.p)
    if config[HSIS]:
        if bl.be1.p is None:
            return None, bl.be2
        elif bl.be2.p is None:
            return bl.be2, None
        else:
            return (bl.be1, bl.be2) if bl.be1.p < bl.be2.p else (bl.be2, bl.be1)
    return bl.be1, bl.be2


def get_unique_breakage_locations(br_locations_by_chr_indexes, config):
    if not config[HIIS] and not config[HSIS]:
        return br_locations_by_chr_indexes
    result = defaultdict(list)
    entries = defaultdict(lambda: defaultdict(list))
    for chr_index, br_locations in br_locations_by_chr_indexes.items():
        for bl in br_locations:
            entries[breakage_location_hashable_representation(breakage_location=bl, config=config)][chr_index].append(bl)
    for key, values in entries.items():
        chr_choices = list(values.keys())
        chr_cnts = []
        for chr_index in chr_choices:
            chr_cnts.append(len(values[chr_index]))
        total_cnt = sum(chr_cnts)
        chr_probs = [cnt / total_cnt for cnt in chr_cnts]
        chr_index = np.random.choice(a=chr_choices, p=chr_probs)
        bl = values[chr_index][np.random.choice(len(values[chr_index]))]
        result[chr_index].append(bl)
    for chr_index in sorted(result.keys()):
        br_locations = result[chr_index]
        s_br_locations = sorted(br_locations, key=lambda l: l.bi)
        result[chr_index] = s_br_locations
    return result


def get_available_breakage_locations(chromosome, config):
    breakage_locations = []
    if len(chromosome) == 0:
        return breakage_locations
    ###
    #
    # processing left telomere
    #
    ###
    start_position, haplotype = chromosome[0].start_position, chromosome[0].extra[HAPLOTYPE]
    be1 = BreakageExtremity(p=None, h=None)
    be2 = BreakageExtremity(p=start_position, h=haplotype)
    b_location = BreakageLocation(be1=be1, be2=be2, bi=0)
    if suitable_for_breakage(breakage_location=b_location, config=config):
        breakage_locations.append(b_location)
    ###
    #
    # processing all but telomeres in a chromosome
    #
    ###
    for i in range(len(chromosome) - 1):
        p1 = chromosome[i].end_position
        p1_haplotype = chromosome[i].extra[HAPLOTYPE]
        be1 = BreakageExtremity(p=p1, h=p1_haplotype)
        p2 = chromosome[i + 1].start_position
        p2_haplotype = chromosome[i].extra[HAPLOTYPE]
        be2 = BreakageExtremity(p=p2, h=p2_haplotype)
        b_location = BreakageLocation(be1=be1, be2=be2, bi=i + 1)
        if suitable_for_breakage(breakage_location=b_location, config=config):
            breakage_locations.append(b_location)
    ###
    #
    # processing right telomere
    #
    ###
    end_position, haplotype = chromosome[-1].end_position, chromosome[-1].extra[HAPLOTYPE]
    be1 = BreakageExtremity(p=end_position, h=haplotype)
    be2 = BreakageExtremity(p=None, h=None)
    b_location = BreakageLocation(be1=be1, be2=be2, bi=len(chromosome))
    if suitable_for_breakage(breakage_location=b_location, config=config):
        breakage_locations.append(b_location)
    return breakage_locations


def generate_reversal_mutation(genome, config):
    try_cnt = 0
    while True:
        try_cnt += 1
        if try_cnt >= config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.REVERSAL][REVERSAL_TRY_LIMIT]:
            return Mutation(mutation_type=MutationType.REVERSAL, mutation_data=None)
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_locations = get_available_breakage_locations(chromosome=chromosome, config=config)
        br_locations = get_unique_breakage_locations(br_locations_by_chr_indexes={chromosome_index: br_locations}, config=config)[chromosome_index]
        if len(br_locations) == 0:
            continue
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.REVERSAL][REVERSAL_TELOMERE]:
            if breakage_location_is_on_telomere(breakage_location=br_locations[0]):
                br_locations = br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=br_locations[-1]):
                br_locations = br_locations[:-1]
        if len(br_locations) >= 2:
            break
    reversal_start_index_i = np.random.randint(low=0, high=len(br_locations) - 1)
    rev_start_bl = br_locations[reversal_start_index_i]
    if reversal_start_index_i + 1 == len(br_locations) - 1:
        reversal_end_index_i = len(br_locations) - 1
    else:
        reversal_end_index_i = np.random.randint(low=reversal_start_index_i + 1, high=len(br_locations))
    rev_end_bl = br_locations[reversal_end_index_i]
    return Mutation(mutation_type=MutationType.REVERSAL,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        REVERSAL_START_INDEX: rev_start_bl.bi,
                        REVERSAL_END_INDEX: rev_end_bl.bi,
                        LOCATIONS: [rev_start_bl, rev_end_bl],
                        POSITIONS: [rev_start_bl.be1.p, rev_start_bl.be2.p, rev_end_bl.be1.p, rev_end_bl.be2.p],
                        HAPLOTYPES: [rev_start_bl.be1.h, rev_start_bl.be2.h, rev_end_bl.be1.h, rev_end_bl.be2.h]
                    })


def generate_duplication_mutation(genome, config):
    try_cnt = 0
    while True:
        try_cnt += 1
        if try_cnt >= config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_TRY_LIMIT]:
            return Mutation(mutation_type=MutationType.DUPLICATION, mutation_data=None)
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_locations = get_available_breakage_locations(chromosome=chromosome, config=config)
        br_locations = get_unique_breakage_locations(br_locations_by_chr_indexes={chromosome_index: br_locations}, config=config)[chromosome_index]
        if len(br_locations) == 0:
            continue
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_TELOMERE]:
            if breakage_location_is_on_telomere(breakage_location=br_locations[0]):
                br_locations = br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=br_locations[-1]):
                br_locations = br_locations[:-1]
        if len(br_locations) >= 2:
            break
    duplication_start_index_i = np.random.randint(low=0, high=len(br_locations) - 1)
    dup_start_bl = br_locations[duplication_start_index_i]
    if duplication_start_index_i + 1 == len(br_locations) - 1:
        duplication_end_index_i = len(br_locations) - 1
    else:
        duplication_end_index_i = np.random.randint(low=duplication_start_index_i + 1, high=len(br_locations))
    dup_end_bl = br_locations[duplication_end_index_i]
    tandem_dup_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_TANDEM_PROBABILITY]
    tandem_duplication = np.random.choice(a=[True, False], p=[tandem_dup_prob, 1 - tandem_dup_prob])
    if tandem_duplication:
        dup_insert_bl = dup_end_bl
    else:
        downstream_allowed = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_ALLOW_DOWNSTREAM]
        downstream_end_index = max(duplication_start_index_i - 1, 0)
        down_indexes = list(range(duplication_start_index_i))[:downstream_end_index]
        upstream_start_index = min(duplication_end_index_i, len(br_locations) - 1)
        up_indexes = list(range(len(br_locations)))[upstream_start_index:]
        indexes = []
        if downstream_allowed:
            indexes.extend(down_indexes)
        indexes.extend(up_indexes)
        if len(indexes) == 0:
            dup_insert_index_i = duplication_end_index_i
        else:
            dup_insert_index_i = np.random.choice(indexes)
        dup_insert_bl = br_locations[dup_insert_index_i]
    reversed_dup_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_REVERSED_PROBABILITY]
    reversed_duplication = np.random.choice(a=[True, False], p=[reversed_dup_prob, 1 - reversed_dup_prob])
    return Mutation(mutation_type=MutationType.DUPLICATION,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        DUPLICATION_START_INDEX: dup_start_bl.bi,
                        DUPLICATION_END_INDEX: dup_end_bl.bi,
                        DUPLICATION_INSERT_INDEX: dup_insert_bl.bi,
                        LOCATIONS: [dup_start_bl, dup_end_bl, dup_insert_bl],
                        DUPLICATION_REVERSED: reversed_duplication,
                        POSITIONS: [dup_start_bl.be1.p, dup_start_bl.be2.p, dup_end_bl.be1.p, dup_end_bl.be2.p, dup_insert_bl.be1.p, dup_insert_bl.be2.p],
                        HAPLOTYPES: [dup_start_bl.be1.h, dup_start_bl.be2.h, dup_end_bl.be1.h, dup_end_bl.be2.h, dup_insert_bl.be1.h, dup_insert_bl.be2.h]
                    })


def generate_deletion_mutation(genome, config):
    try_cnt = 0
    while True:
        try_cnt += 1
        if try_cnt >= config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DELETION][DELETION_TRY_LIMIT]:
            return Mutation(mutation_type=MutationType.DELETION, mutation_data=None)
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_locations = get_available_breakage_locations(chromosome=chromosome, config=config)
        br_locations = get_unique_breakage_locations(br_locations_by_chr_indexes={chromosome_index: br_locations}, config=config)[chromosome_index]
        if len(br_locations) == 0:
            continue
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DELETION][DELETION_TELOMERE]:
            if breakage_location_is_on_telomere(breakage_location=br_locations[0]):
                br_locations = br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=br_locations[-1]):
                br_locations = br_locations[:-1]
        if len(br_locations) >= 2:
            break
    deletion_start_index_i = np.random.randint(low=0, high=len(br_locations) - 1)
    del_start_bl = br_locations[deletion_start_index_i]
    if deletion_start_index_i + 1 == len(br_locations) - 1:
        deletion_end_index_i = len(br_locations) - 1
    else:
        deletion_end_index_i = np.random.randint(low=deletion_start_index_i + 1, high=len(br_locations))
    del_end_bl = br_locations[deletion_end_index_i]
    return Mutation(mutation_type=MutationType.DELETION,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        DELETION_START_INDEX: del_start_bl.bi,
                        DELETION_END_INDEX: del_end_bl.bi,
                        LOCATIONS: [del_start_bl, del_end_bl],
                        POSITIONS: [del_start_bl.be1.p, del_start_bl.be2.p, del_end_bl.be1.p, del_end_bl.be2.p],
                        HAPLOTYPES: [del_start_bl.be1.h, del_start_bl.be2.h, del_end_bl.be1.h, del_end_bl.be2.h]
                    })


def chromosomes_are_mates(chromosome1, chromosome2):
    chr1_chr = [s.chromosome for s in chromosome1]
    chr2_chr = [s.chromosome for s in chromosome2]
    chr1_counter = Counter(chr1_chr)
    chr2_counter = Counter(chr2_chr)
    chr1_most_common_chr = chr1_counter.most_common(n=1)[0][0]
    chr2_most_common_chr = chr2_counter.most_common(n=1)[0][0]
    if chr1_most_common_chr != chr2_most_common_chr:
        return False
    return True


def generate_translocation_mutation(genome, config):
    try_cnt = 0
    while True:
        try_cnt += 1
        if try_cnt >= config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.TRANSLOCATION][TRANSLOCATION_TRY_LIMIT]:
            return Mutation(mutation_type=MutationType.TRANSLOCATION, mutation_data=None)
        chromosome_1_index, chromosome_2_index = np.random.randint(low=0, high=len(genome), size=2)
        if chromosome_1_index == chromosome_2_index:
            continue
        chromosome_1 = genome[chromosome_1_index]
        chromosome_2 = genome[chromosome_2_index]
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.TRANSLOCATION][HOMOZYGOUS] and \
                chromosomes_are_mates(chromosome1=chromosome_1, chromosome2=chromosome_2):
            continue
        chr1_br_locations = get_available_breakage_locations(chromosome=chromosome_1, config=config)
        if len(chr1_br_locations) == 0:
            continue
        chr2_br_locations = get_available_breakage_locations(chromosome=chromosome_2, config=config)
        if len(chr2_br_locations) == 0:
            continue
        bls_by_chr_indexes = {
            chromosome_1_index: chr1_br_locations,
            chromosome_2_index: chr2_br_locations,
        }
        u_bls_by_chr_indexes = get_unique_breakage_locations(br_locations_by_chr_indexes=bls_by_chr_indexes, config=config)
        chr1_br_locations = u_bls_by_chr_indexes[chromosome_1_index]
        chr2_br_locations = u_bls_by_chr_indexes[chromosome_2_index]
        if len(chr1_br_locations) == 0 or len(chr2_br_locations) == 0:
            continue
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.TRANSLOCATION][TRANSLOCATION_CHR1_TELOMERE]:
            if breakage_location_is_on_telomere(breakage_location=chr1_br_locations[0]):
                chr1_br_locations = chr1_br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=chr1_br_locations[-1]):
                chr1_br_locations = chr1_br_locations[:-1]
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.TRANSLOCATION][TRANSLOCATION_CHR2_TELOMERE]:
            if breakage_location_is_on_telomere(breakage_location=chr2_br_locations[0]):
                chr2_br_locations = chr2_br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=chr2_br_locations[-1]):
                chr2_br_locations = chr2_br_locations[:-1]
        if len(chr1_br_locations) >= 1 and len(chr2_br_locations) >= 1:
            break
    chr1_index = np.random.randint(low=0, high=len(chr1_br_locations))
    chr1_bl = chr1_br_locations[chr1_index]
    chr2_index = np.random.randint(low=0, high=len(chr2_br_locations))
    chr2_bl = chr2_br_locations[chr2_index]
    cc = np.random.choice([True, False])
    return Mutation(mutation_type=MutationType.TRANSLOCATION,
                    mutation_data={
                        CHROMOSOME_1_INDEX: chromosome_1_index,
                        CHROMOSOME_2_INDEX: chromosome_2_index,
                        TRANSLOCATION_CHR1_INDEX: chr1_bl.bi,
                        TRANSLOCATION_CHR2_INDEX: chr2_bl.bi,
                        LOCATIONS: [chr1_bl, chr2_bl],
                        TRANSLOCATION_CC: cc,
                        POSITIONS: [chr1_bl.be1.p, chr1_bl.be2.p, chr2_bl.be1.p, chr2_bl.be2.p],
                        HAPLOTYPES: [chr1_bl.be1.h, chr1_bl.be2.h, chr2_bl.be1.h, chr2_bl.be2.h]
                    })


def reversed_part(chr_part):
    return ChrPart(chr_index=chr_part.chr_index, pi=chr_part.pi, reverse=not chr_part.reverse)


def generate_chromothripsis_mutation(genome, config):
    try_cnt = 0
    min_chr_cnt = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_MIN_CHR_CNT]
    if len(genome) < min_chr_cnt:
        return Mutation(mutation_type=MutationType.CHROMOTHRIPSIS, mutation_data=None)
    max_chr_cnt = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_MAX_CHR_CNT]
    while True:
        try_cnt += 1
        if try_cnt >= config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_TRY_LIMIT]:
            return Mutation(mutation_type=MutationType.CHROMOTHRIPSIS, mutation_data=None)
        chr_cnt = np.random.randint(low=min_chr_cnt, high=min(len(genome) + 1, max_chr_cnt + 1))
        chromosomes_indexes = sorted(np.random.choice(a=list(range(len(genome))), size=chr_cnt, replace=False))
        chr_by_indexes = {chr_index: genome[chr_index] for chr_index in chromosomes_indexes}
        chromosomes_operational = [genome[chr_index] for chr_index in chromosomes_indexes]
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][HOMOZYGOUS]:
            all_non_mates = True
            for chr1, chr2 in itertools.combinations(chromosomes_operational, r=2):
                if chromosomes_are_mates(chr1, chr2):
                    all_non_mates = False
                    break
            if not all_non_mates:
                continue
        all_br_locations_by_chr_indexes = {}
        all_chr_have_br_locations = True
        for chr_index, chromosome in chr_by_indexes.items():
            chromosome_br_locations = get_available_breakage_locations(chromosome=chromosome, config=config)
            ###
            #
            # Telomeres must always be excluded as potential "cut" positions, as in chromothripsis the idea is to "cut chromosome(s) into a set of non-empty "spans" of segments
            #   and cutting on the telomere does not comply with that
            #
            ###
            if len(chromosome_br_locations) == 0:
                all_chr_have_br_locations = False
                break
            if breakage_location_is_on_telomere(breakage_location=chromosome_br_locations[0]):
                chromosome_br_locations = chromosome_br_locations[1:]
            if breakage_location_is_on_telomere(breakage_location=chromosome_br_locations[-1]):
                chromosome_br_locations = chromosome_br_locations[:-1]
            if len(chromosome_br_locations) == 0:
                all_chr_have_br_locations = False
                break
            all_br_locations_by_chr_indexes[chr_index] = chromosome_br_locations
        if not all_chr_have_br_locations:
            continue
        all_br_locations_by_chr_indexes = get_unique_breakage_locations(br_locations_by_chr_indexes=all_br_locations_by_chr_indexes, config=config)
        for chr_index in chromosomes_indexes:
            if chr_index not in all_br_locations_by_chr_indexes:
                all_chr_have_br_locations = False
                break
            br_locations = all_br_locations_by_chr_indexes[chr_index]
            if len(br_locations) == 0:
                all_chr_have_br_locations = False
                break
        if not all_chr_have_br_locations:
            continue
        total_br_locations_cnt = sum(map(lambda l: len(l), all_br_locations_by_chr_indexes.values()))
        if total_br_locations_cnt < config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_MIN_BR_CNT]:
            continue
        min_br_cnt = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_MIN_BR_CNT]
        max_br_cnt = min(config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_MAX_BR_CNT], total_br_locations_cnt)
        br_cnt = np.random.randint(low=min_br_cnt, high=max_br_cnt + 1)
        br_cnt_by_chr_indexes = defaultdict(int)
        ###
        #
        # producing the number of breakages per chromosome
        #
        ###
        # every chromosomes needs at least 1 breakage to happen on it
        for chr_index, _ in zip(chromosomes_indexes, range(br_cnt)):
            br_cnt_by_chr_indexes[chr_index] = 1
        if br_cnt <= len(chromosomes_indexes):
            # if the number of operational chromosomes were less or equal to the number of breakages -- each chromosome stays with 1/0 breakages to its name
            #   indexes of chromosomes with 0 breakages are cut from the index list
            chromosomes_indexes = chromosomes_indexes[:br_cnt]
        else:
            # br_cnt now has the number of breakages we need to distribute among the chromosomes, keeping in mind that each chromosomes already got 1 breakage
            br_cnt -= len(chromosomes_indexes)
            br_capacity_by_chr_index = {chr_index: len(all_br_locations_by_chr_indexes[chr_index]) - 1 for chr_index in chromosomes_indexes}
            while br_cnt > 0:
                available_chromosomes_indexes = get_indexes_chromosomes_available_for_breakage(chromosomes_indexes=chromosomes_indexes,
                                                                                               chr_br_capacity_by_chr_indexes=br_capacity_by_chr_index)
                chr_index = np.random.choice(a=available_chromosomes_indexes)
                br_cnt_by_chr_indexes[chr_index] += 1
                br_capacity_by_chr_index[chr_index] -= 1
                br_cnt -= 1
        ###
        #
        # producing the actual breakage locations per chromosome
        #
        ###
        br_locations_by_chr_indexes = {}
        for chr_index in chromosomes_indexes:
            chr_br_cnt = br_cnt_by_chr_indexes[chr_index]
            chr_br_locations = all_br_locations_by_chr_indexes[chr_index]
            selected_br_locations_index = sorted(np.random.choice(a=list(range(len(chr_br_locations))), size=chr_br_cnt, replace=False))
            selected_br_locations = [chr_br_locations[index] for index in selected_br_locations_index]
            br_locations_by_chr_indexes[chr_index] = selected_br_locations
        ###
        #
        # producing parts that can be randomly shuffled around in the further target structure
        #   only concerned about chromosomes, that are chopped up (i.e., chromosomes which indexes are present in `chromosomes_indexes`)
        #
        ###
        non_telomere_parts_by_chr_index = {}
        telomere_parts_by_chr_index = {}
        for chr_index in chromosomes_indexes:
            chr_non_telomere_parts = [ChrPart(chr_index=chr_index, pi=part_index, reverse=False) for part_index in range(1, br_cnt_by_chr_indexes[chr_index])]
            chr_telomere_parts = [ChrPart(chr_index=chr_index, pi=0, reverse=False), ChrPart(chr_index=chr_index, pi=br_cnt_by_chr_indexes[chr_index], reverse=False)]
            non_telomere_parts_by_chr_index[chr_index] = chr_non_telomere_parts
            telomere_parts_by_chr_index[chr_index] = chr_telomere_parts
        ###
        #
        #
        #
        ###
        parts_to_shuffle = []
        for part in non_telomere_parts_by_chr_index.values():
            parts_to_shuffle.extend(part)
        if config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_TELOMERE]:
            for part in telomere_parts_by_chr_index.values():
                parts_to_shuffle.extend(part)
        ###
        #
        # introducing duplications/deletions for parts
        #
        ###
        result_parts_to_shuffle = []
        for part in parts_to_shuffle:
            possible_to_duplicate = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_DUPLICATIONS]
            possible_to_delete = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_DELETIONS]
            dup_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_DUPLICATIONS_PROB]
            del_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_DELETIONS_PROB]
            duplicate = False if not possible_to_duplicate else np.random.choice(a=[True, False], p=[dup_prob, 1 - dup_prob])
            delete = False if not possible_to_delete else np.random.choice(a=[True, False], p=[del_prob, 1 - del_prob])
            if duplicate and delete:
                rel_dup_probability = dup_prob / (dup_prob + del_prob)
                duplicate = np.random.choice(a=[True, False], p=[rel_dup_probability, 1 - rel_dup_probability])
                delete = not duplicate
            if duplicate:
                result_parts_to_shuffle.append(deepcopy(part))
            if delete:
                continue
            result_parts_to_shuffle.append(part)
        parts = result_parts_to_shuffle
        ###
        #
        # introducing possible reversals for resulting parts, that will be shuffled
        #
        ###
        result_parts_to_shuffle = []
        for part in parts:
            possible_to_reverse = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_REVERSALS]
            rev_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_REVERSALS_PROB]
            reverse = False if not possible_to_reverse else np.random.choice(a=[True, False], p=[rev_prob, 1 - rev_prob])
            part = part if not reverse else reversed_part(chr_part=part)
            result_parts_to_shuffle.append(part)
        parts = result_parts_to_shuffle
        ###
        #
        #
        #
        ###
        np.random.shuffle(parts)
        ###
        #
        # constructing target chromosome structures
        #
        ###
        target_structures_by_chr_indexes = defaultdict(list)
        for part in parts:
            target_chr_index = np.random.choice(chromosomes_indexes)
            target_structures_by_chr_indexes[target_chr_index].append(part)
        ###
        #
        # if telomeres were allowed to be used as parts -- they would have been in the parts list to begin with and thus no further work would have been needed.
        #   Otherwise, telomeres needs to be ap/pre-pended to target chromosomal list structures
        #
        ###
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_TELOMERE]:
            if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.CHROMOTHRIPSIS][CHROMOTHRIPSIS_TELOMERE_PAIRING_PRESERVE]:

                for chr_index in chromosomes_indexes:
                    flip = np.random.choice(a=[True, False], p=[0.5, 0.5])
                    if flip:
                        l_telomere, r_telomere = telomere_parts_by_chr_index[chr_index]
                        r_l_telomere = reversed_part(chr_part=l_telomere)
                        r_r_telomere = reversed_part(chr_part=r_telomere)
                        telomere_parts_by_chr_index[chr_index] = [r_r_telomere, r_l_telomere]
                l_telomeres = [telomere_parts_by_chr_index[chr_index][0] for chr_index in chromosomes_indexes]
                r_telomeres = [telomere_parts_by_chr_index[chr_index][1] for chr_index in chromosomes_indexes]
                np.random.shuffle(r_telomeres)
                for chr_index, l_telomere, r_telomere in zip(chromosomes_indexes, l_telomeres, r_telomeres):
                    telomere_parts_by_chr_index[chr_index] = [l_telomere, r_telomere]
            for chr_index in chromosomes_indexes:
                l_telomere, r_telomere = telomere_parts_by_chr_index[chr_index]
                target_structures_by_chr_indexes[chr_index].insert(0, l_telomere)
                target_structures_by_chr_indexes[chr_index].append(r_telomere)
        ###
        #
        # creating per chromosome cut indexes lists. Chromosomes that were not part of chromothripsis event have empty segment cutting indexes lists.
        #
        ###
        per_chrs_segment_indexes = []
        for chr_index, _ in enumerate(genome):
            if chr_index in chromosomes_indexes:
                br_locations = br_locations_by_chr_indexes[chr_index]
                segment_indexes = [brl.bi for brl in br_locations]
                per_chrs_segment_indexes.append(segment_indexes)
            else:
                per_chrs_segment_indexes.append([])
        ###
        #
        # explicitly adding chromosomes that were not selected for "cutting" into the target structure
        #
        ###
        target_structure = []
        for chr_index, _ in enumerate(genome):
            if chr_index in chromosomes_indexes:
                target_structure.append(target_structures_by_chr_indexes[chr_index])
            else:
                # dummy target for chromosomes that were not part of the "cutting" spree.
                target_structure.append([ChrPart(chr_index=chr_index, pi=0, reverse=False)])
        ###
        #
        #
        #
        ###
        all_br_locations = []
        for chr_index in chromosomes_indexes:
            all_br_locations.extend(br_locations_by_chr_indexes[chr_index])
        all_positions = []
        all_haplotypes = []
        for brl in all_br_locations:
            all_positions.append(brl.be1.p)
            all_positions.append(brl.be2.p)
            all_haplotypes.append(brl.be1.h)
            all_haplotypes.append(brl.be2.h)
        ###
        #
        # have to exist the True loop, since we've completed all we've wanted and everything for mutation to be created exists
        #
        ###
        break
    return Mutation(mutation_type=MutationType.CHROMOTHRIPSIS,
                    mutation_data={
                        CHROMOSOME_INDEXES: chromosomes_indexes,
                        LOCATIONS: all_br_locations,
                        POSITIONS: all_positions,
                        CHROMOTHRIPSIS_PER_CHRS_SEGMENT_INDEXES: per_chrs_segment_indexes,
                        CHROMOTHRIPSIS_TARGET_STRUCTURE: target_structure,
                    })


def get_indexes_chromosomes_available_for_breakage(chromosomes_indexes, chr_br_capacity_by_chr_indexes):
    return [chr_index for chr_index in chromosomes_indexes if chr_br_capacity_by_chr_indexes[chr_index] > 0]


def generate_chromosome_duplication_mutation(genome, config):
    chromosome_index = np.random.randint(low=0, high=len(genome))
    return Mutation(mutation_type=MutationType.CHROMOSOME_DUP,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        LOCATIONS: [],
                    })


def generate_chromosome_deletion_mutation(genome, config):
    chromosome_index = np.random.randint(low=0, high=len(genome))
    return Mutation(mutation_type=MutationType.CHROMOSOME_DEL,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        LOCATIONS: [],
                    })


def generate_wgd_mutation(genome, config):
    return Mutation(mutation_type=MutationType.WGD,
                    mutation_data={
                        LOCATIONS: [],
                    })


def apply_mutation(genome, mutation, inplace=True):
    if mutation.mutation_type not in MUTATION_TYPES_to_MUTATING_FUNCTIONS:
        raise ValueError("Unsupported mutation type {mt}".format(mt=repr(mutation.mutation_type)))
    mutating_function = MUTATION_TYPES_to_MUTATING_FUNCTIONS[mutation.mutation_type]
    return mutating_function(genome=genome, mutation=mutation, inplace=inplace)


def apply_reversal_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chr_index = mutation.mutation_data[CHROMOSOME_INDEX]
    chromosome = new_genome[chr_index]
    rev_start_index = mutation.mutation_data[REVERSAL_START_INDEX]
    rev_end_index = mutation.mutation_data[REVERSAL_END_INDEX]
    new_chromosome = reverse_segments(chromosome=chromosome,
                                      start_segment_index=rev_start_index,
                                      end_segment_index=rev_end_index,
                                      copy=False)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_duplication_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chr_index = mutation.mutation_data[CHROMOSOME_INDEX]
    chromosome = new_genome[chr_index]
    dup_start_index = mutation.mutation_data[DUPLICATION_START_INDEX]
    dup_end_index = mutation.mutation_data[DUPLICATION_END_INDEX]
    dup_insert_index = mutation.mutation_data[DUPLICATION_INSERT_INDEX]
    duplication_reverse = mutation.mutation_data[DUPLICATION_REVERSED]
    new_chromosome = duplicate_segments(chromosome=chromosome,
                                        start_segment_index=dup_start_index,
                                        end_segment_index=dup_end_index,
                                        insert_segment_index=dup_insert_index,
                                        reverse_copy=duplication_reverse,
                                        copy=False)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_deletion_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chr_index = mutation.mutation_data[CHROMOSOME_INDEX]
    chromosome = new_genome[chr_index]
    del_start_index = mutation.mutation_data[DELETION_START_INDEX]
    del_end_index = mutation.mutation_data[DELETION_END_INDEX]
    new_chromosome = delete_segments(chromosome=chromosome,
                                     start_segment_index=del_start_index,
                                     end_segment_index=del_end_index,
                                     copy=False)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_translocation_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chr1_index = mutation.mutation_data[CHROMOSOME_1_INDEX]
    chr2_index = mutation.mutation_data[CHROMOSOME_2_INDEX]
    chromosome_1 = new_genome[chr1_index]
    chromosome_2 = new_genome[chr2_index]
    chr1_trans_index = mutation.mutation_data[TRANSLOCATION_CHR1_INDEX]
    chr2_trans_index = mutation.mutation_data[TRANSLOCATION_CHR2_INDEX]
    cc = mutation.mutation_data[TRANSLOCATION_CC]
    new_chr1, new_chr2 = translocate_segments(chromosome1=chromosome_1, chromosome2=chromosome_2,
                                              chromosome1_segment_index=chr1_trans_index,
                                              chromosome2_segment_index=chr2_trans_index,
                                              cc=cc,
                                              copy=False)
    new_genome[chr1_index] = new_chr1
    new_genome[chr2_index] = new_chr2
    return new_genome


def apply_chromothripsis_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    per_chrs_segment_indexes = mutation.mutation_data[CHROMOTHRIPSIS_PER_CHRS_SEGMENT_INDEXES]
    target_structure = mutation.mutation_data[CHROMOTHRIPSIS_TARGET_STRUCTURE]
    result_genome = split_and_reassemble_chromosomes(chromosomes=new_genome,
                                                     per_chrs_segment_indexes=per_chrs_segment_indexes,
                                                     target_structures=target_structure)
    return result_genome


def apply_chromosome_duplication_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chromosome_index = mutation.mutation_data[CHROMOSOME_INDEX]
    duplicated_chromosome = duplicate_chromosome(chromosome=genome[chromosome_index])
    new_genome.append(duplicated_chromosome)
    return new_genome


def apply_chromosome_deletion_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    chromosome_index = mutation.mutation_data[CHROMOSOME_INDEX]
    del new_genome[chromosome_index]
    return new_genome


def apply_wgd_mutation(genome, mutation, inplace=True):
    new_genome = genome if inplace else deepcopy(genome)
    genome_copy = duplicate_genome(chromosomes=new_genome)
    new_genome.extend(genome_copy)
    return new_genome


MUTATION_TYPES_to_GENERATING_FUNCTIONS = {
    MutationType.REVERSAL: generate_reversal_mutation,
    MutationType.DELETION: generate_deletion_mutation,
    MutationType.DUPLICATION: generate_duplication_mutation,
    MutationType.TRANSLOCATION: generate_translocation_mutation,
    MutationType.CHROMOTHRIPSIS: generate_chromothripsis_mutation,
    MutationType.CHROMOSOME_DUP: generate_chromosome_duplication_mutation,
    MutationType.CHROMOSOME_DEL: generate_chromosome_deletion_mutation,
    MutationType.WGD: generate_wgd_mutation,
}

MUTATION_TYPES_to_MUTATING_FUNCTIONS = {
    MutationType.REVERSAL: apply_reversal_mutation,
    MutationType.DELETION: apply_deletion_mutation,
    MutationType.DUPLICATION: apply_duplication_mutation,
    MutationType.TRANSLOCATION: apply_translocation_mutation,
    MutationType.CHROMOTHRIPSIS: apply_chromothripsis_mutation,
    MutationType.CHROMOSOME_DUP: apply_chromosome_duplication_mutation,
    MutationType.CHROMOSOME_DEL: apply_chromosome_deletion_mutation,
    MutationType.WGD: apply_wgd_mutation,
}


def generate_mutated_genome(starting_genome, mutation_cnt, config=MUTATION_CONFIG, save_intermediate_genome=False):
    if MUTATED_EXTREMITIES not in config:
        config[MUTATED_EXTREMITIES] = set()
    current_genome = deepcopy(starting_genome)
    history = {
        GENOMES: [],
        MUTATIONS: [],
    }
    for mut_cnt in range(mutation_cnt):
        mut_types = []
        mut_prob = []
        for mut, prob in config[MUTATIONS_PROBABILITIES_BY_MUTATIONS_TYPES].items():
            mut_types.append(mut)
            mut_prob.append(prob)
        mutation_type = np.random.choice(a=mut_types,
                                         p=mut_prob)
        mutation = generate_mutation(mutation_type=mutation_type,
                                     genome=current_genome,
                                     config=config)
        history[MUTATIONS].append(mutation)
        if mutation.mutation_data is not None:
            current_genome = apply_mutation(genome=current_genome,
                                            mutation=mutation,
                                            inplace=True)
            for bl in mutation.mutation_data[LOCATIONS]:
                for ep in [bl.be1, bl.be2]:
                    if ep.p is not None:
                        if config[HIIS]:
                            config[MUTATED_EXTREMITIES].add(ep.p)
                        elif config[HSIS]:
                            config[MUTATED_EXTREMITIES].add((ep.p, ep.h))
        if save_intermediate_genome or mut_cnt == (mutation_cnt - 1):
            history[GENOMES].append(current_genome)
            current_genome = deepcopy(current_genome)
    return history


def generate_random_phylogeny_tree(n, return_root=True):
    nodes = list(range(n))
    tree, root = rec_generate_random_phylogeny_tree(nodes=nodes)
    if return_root:
        return tree, root
    return tree


@lru_cache(maxsize=1000)
def tree_cnt(i):
    @lru_cache(maxsize=500)
    def f(n):
        odd = n % 2 == 1
        half = n // 2
        if odd:
            half_cnt = tree_cnt(half)
            return ((half_cnt + 1) * half_cnt) // 2
        else:
            return tree_cnt(half - 1) * tree_cnt(half)
    if i == 0:
        return 1
    if i == 1:
        return 1
    if i == 2:
        return 1
    upper_bound = int(math.ceil(i / 2)) - 1
    result = 0
    for t_i in range(upper_bound):
        result += (tree_cnt(i - t_i - 1) * tree_cnt(t_i))
    result += f(i)
    return result


@lru_cache(maxsize=1000)
def get_root_selection_probabilities(nodes_cnt):
    cnts = []
    middle = int(math.ceil(nodes_cnt / 2))
    for i in range(middle):
        cnts.append(tree_cnt(nodes_cnt - i - 1) * tree_cnt(i))
    overall_cnt = sum(cnts)
    return [cnt / overall_cnt for cnt in cnts]


def rec_generate_random_phylogeny_tree(nodes):
    result = nx.DiGraph()
    if len(nodes) == 1:
        root = nodes[0]
        result.add_node(root)
        return result, root
    middle = int(math.ceil(len(nodes) / 2))
    possible_roots = nodes[:middle]
    root_selection_prob = get_root_selection_probabilities(nodes_cnt=len(nodes))
    root = np.random.choice(a=possible_roots, p=root_selection_prob)
    root_index = nodes.index(root)
    rt, rt_root = rec_generate_random_phylogeny_tree(nodes=nodes[root_index + 1:])
    if root_index > 0:
        lt, lt_root = rec_generate_random_phylogeny_tree(nodes=nodes[:root_index])
        result = nx.compose(G=lt, H=rt)
        result.add_node(root)
        result.add_edge(u=root, v=lt_root)
        result.add_edge(u=root, v=rt_root)
    else:
        result = rt
        result.add_node(root)
        result.add_edge(u=root, v=rt_root)
    return result, root


def get_adjacencies_from_genome_old(genome, is_reference=False, ref_adjacencies=None):
    if ref_adjacencies is None:
        ref_adjacencies = set()
    result = {}
    for chromosome in genome:
        for s1, s2 in zip(chromosome[:-1], chromosome[1:]):
            adjacency = Adjacency(position1=s1.end_position, position2=s2.start_position, adjacency_type=AdjacencyType.NOVEL)
            if adjacency.idx in ref_adjacencies or is_reference:
                adjacency.adjacency_type = AdjacencyType.REFERENCE
                ref_adjacencies.add(adjacency.idx)
            if s1.end_position == adjacency.position2:
                phasing = (s1.extra.get(HAPLOTYPE, Haplotype.UNKNOWN), s2.extra.get(HAPLOTYPE, Haplotype.UNKNOWN))
            else:
                phasing = (s2.extra.get(HAPLOTYPE, Haplotype.UNKNOWN), s1.extra.get(HAPLOTYPE, Haplotype.UNKNOWN))
            phasing = HAPLOTYPE_PAIRS_TO_PHASING[phasing]
            if adjacency.idx not in result:
                result[adjacency.idx] = defaultdict(list)
            result[adjacency.idx][phasing].append(adjacency)
    return result


def get_scn_profile_from_genome(genome):
    result = defaultdict(lambda: defaultdict(list))
    for chromosome in genome:
        for segment in chromosome:
            haplotype = segment.extra.get(HAPLOTYPE, Haplotype.UNKNOWN)
            if segment.is_reversed:
                segment_idx = reverse_segment(segment=segment).idx  # internal changes are made to segment on `reverse_segment` function call, have to undo them
                reverse_segment(segment=segment)
            else:
                segment_idx = segment.idx
            result[segment_idx][haplotype].append(segment)
    return result


def get_unphased_adjacency_cn(adjacency_id, acnp, default=0):
    if adjacency_id not in acnp:
        return default
    total_cnt = 0
    for phasing in Phasing:
        if phasing not in acnp[adjacency_id]:
            continue
        value = acnp[adjacency_id][phasing]
        if isinstance(value, list):
            value = len(value)
        total_cnt += value
    return total_cnt


def get_unphased_segment_cn(segment_id, scnp, default=0):
    if segment_id not in scnp:
        return default
    total_cnt = 0
    for haplotype in Haplotype:
        if haplotype not in scnp[segment_id]:
            continue
        value = scnp[segment_id][haplotype]
        if isinstance(value, list):
            value = len(value)
        total_cnt += value
    return total_cnt


def get_correctly_inferred_present_absent_unphased_adjacencies(ref_acnp, inf_acnp):
    result = get_correctly_inferred_present_absent_unphased_adjacencies(ref_acnp=ref_acnp, inf_acnp=inf_acnp)
    return result.union(get_correctly_inferred_present_unphased_adjacencies(ref_acnp=ref_acnp, inf_acnp=inf_acnp))


def get_correctly_inferred_present_unphased_adjacencies(ref_acnp, inf_acnp):
    result = set()
    adjacency_ids = set(ref_acnp.keys()).union(inf_acnp.keys())
    for a_id in adjacency_ids:
        ref_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=ref_acnp)
        inf_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=inf_acnp)
        if ref_total_cnt > 0 and inf_total_cnt > 0:
            result.add(a_id)
    return result


def get_correctly_inferred_absent_unphased_adjacencies(ref_acnp, inf_acnp):
    result = set()
    adjacency_ids = set(ref_acnp.keys()).union(set(inf_acnp.keys()))
    for a_id in adjacency_ids:
        ref_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=ref_acnp)
        inf_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=inf_acnp)
        if ref_total_cnt == 0 and inf_total_cnt == 0:
            result.add(a_id)
    return result


def get_correctly_inferred_unphased_adjacencies(ref_acnp, inf_acnp):
    result = set()
    adjacency_ids = set(ref_acnp.keys()).union(set(inf_acnp.keys()))
    for a_id in adjacency_ids:
        ref_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=ref_acnp)
        inf_total_cnt = get_unphased_adjacency_cn(adjacency_id=a_id, acnp=inf_acnp)
        if int(ref_total_cnt) == int(inf_total_cnt):
            result.add(a_id)
    return result


def get_correctly_inferred_unphased_segments(ref_scnp, inf_scnp):
    result = set()
    segment_ids = set(ref_scnp.keys()).union(inf_scnp.keys())
    for s_id in segment_ids:
        ref_total_cnt = get_unphased_segment_cn(segment_id=s_id, scnp=ref_scnp)
        inf_total_cnt = get_unphased_segment_cn(segment_id=s_id, scnp=inf_scnp)
        if ref_total_cnt == inf_total_cnt:
            result.add(s_id)
    return result


def generate_fake_adjacencies(positions, number, real_adjacencies=None):
    if real_adjacencies is None:
        real_adjacencies = []
    existing_adjacencies_ids = {a.idx for a in real_adjacencies}
    result = {}
    for _ in range(number):
        while True:
            p1, p2 = np.random.choice(positions, size=2)
            adj = Adjacency(position1=p1, position2=p2, adjacency_type=AdjacencyType.NOVEL, extra={"fake": True})
            if adj.idx in existing_adjacencies_ids:
                continue
            result[adj.idx] = adj
            existing_adjacencies_ids.add(adj.idx)
            break
    return result
