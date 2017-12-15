from collections import Counter, defaultdict
from copy import deepcopy

from dassp.core.structures import Adjacency, AdjacencyType, Haplotype, Phasing, HAPLOTYPE_PAIRS_TO_PHASING
from dassp.simulation.parts import ChromosomeGenerator, MutationType, Mutation, reverse_segments, duplicate_segments, delete_segments, reverse_segment, HAPLOTYPE
import numpy as np

from simulation.parts import translocate_segments

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
MUTATION_PROBABILITIES = "mutation_probabilities"
DUPLICATION_TANDEM_PROBABILITY = "tandem_duplication_probability"
DUPLICATION_REVERSED_PROBABILITY = "reversed_duplication_probability"

CHROMOSOME_INDEX = "chromosome_index"
REVERSAL_START_INDEX = "reversal_start_index"
REVERSAL_END_INDEX = "reversal_end_index"
POSITIONS = "positions"

DUPLICATION_START_INDEX = "duplication_start_index"
DUPLICATION_END_INDEX = "duplication_end_index"
DUPLICATION_INSERT_SHIFT = "duplication_insert_shift"
DUPLICATION_REVERSED = "duplication_reversed"

MUTATION_CONFIG = {
    MUTATIONS_TYPES: [MutationType.REVERSAL,
                      MutationType.DUPLICATION,
                      MutationType.DELETION,
                      MutationType.TRANSLOCATION],
    MUTATION_PROBABILITIES: [0.1,
                             0.425,
                             0.425,
                             0.05],
    PRESERVE_TELOMERES: True,
    HIIS: True,
    HSIS: True,
    MUTATIONS_TYPES_SPECIFICATIONS: {
        MutationType.DELETION: {
            DELETION_TELOMERE: False
        },
        MutationType.REVERSAL: {
            REVERSAL_TELOMERE: False
        },
        MutationType.DUPLICATION: {
            DUPLICATION_TELOMERE: False,
            DUPLICATION_TANDEM_PROBABILITY: 1,
            DUPLICATION_REVERSED_PROBABILITY: 0.25
        },
        MutationType.TRANSLOCATION: {
            TRANSLOCATION_CHR1_TELOMERE: False,
            TRANSLOCATION_CHR2_TELOMERE: False,
            HOMOZYGOUS: False
        }
    }
}


class SimulationManager(object):
    def __init__(self, chrs_cnt, clones_cnt, chrs_size=CHROMOSOMES_SIZE, ab=AB):
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
        self.majority_clone_history = generate_mutated_genome(starting_genome=self.initial_genome, mutation_cnt=self.clonal_mutation_cnt)
        self.majority_clone = self.majority_clone_history["genomes"][-1]
        self.minority_clone_history = generate_mutated_genome(starting_genome=self.majority_clone, mutation_cnt=self.subclonal_mutation_cnt)
        self.minority_clone = self.minority_clone_history["genomes"][-1]


def generate_mutation(mutation_type, genome, mutations_config):
    if mutation_type not in MUTATION_TYPES_to_GENERATING_FUNCTIONS:
        raise ValueError("Unsupported mutation type {mt}".format(mt=repr(mutation_type)))
    generating_function = MUTATION_TYPES_to_GENERATING_FUNCTIONS[mutation_type]
    return generating_function(genome=genome, config=mutations_config)


def get_available_breakage_indexes(chromosome, config):
    result = []
    if len(chromosome) == 0:
        return result
    ###
    #
    # processing left telomere
    #
    ###
    start_position, haplotype = chromosome[0].start_position, chromosome[0].extra[HAPLOTYPE]
    if not config[PRESERVE_TELOMERES]:
        if config[HIIS]:
            if start_position not in config[MUTATED_EXTREMITIES]:
                result.append(0)
        elif config[HSIS]:
            if (start_position, haplotype) not in config[MUTATED_EXTREMITIES]:
                result.append(0)
        else:
            result.append(0)
    ###
    #
    # processing all but telomeres in a chromosome
    #
    ###
    for i in range(len(chromosome) - 1):
        p1 = chromosome[i].end_position
        p1_haplotype = chromosome[i].extra[HAPLOTYPE]
        p2 = chromosome[i + 1].start_position
        p2_haplotype = chromosome[i].extra[HAPLOTYPE]
        if config[HIIS]:
            if p1 not in config[MUTATED_EXTREMITIES] and p2 not in config[MUTATED_EXTREMITIES]:
                result.append(i + 1)
        elif config[HSIS]:
            if (p1, p1_haplotype) not in config[MUTATED_EXTREMITIES] and (p2, p2_haplotype) not in config[MUTATED_EXTREMITIES]:
                result.append(i + 1)
        else:
            result.append(i + 1)
    ###
    #
    # processing right telomere
    #
    ###
    end_position, haplotype = chromosome[-1].end_position, chromosome[-1].extra[HAPLOTYPE]
    if not config[PRESERVE_TELOMERES]:
        if config[HIIS]:
            if end_position not in config[MUTATED_EXTREMITIES]:
                result.append(len(chromosome))
        elif config[HSIS]:
            if (end_position, haplotype) not in config[MUTATED_EXTREMITIES]:
                result.append(len(chromosome))
        else:
            result.append(len(chromosome))
    return result


def generate_reversal(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = get_available_breakage_indexes(chromosome=chromosome, config=config)
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.REVERSAL][REVERSAL_TELOMERE]:
            if br_indexes[0] == 0:
                br_indexes = br_indexes[1:]
            if br_indexes[-1] == len(chromosome):
                br_indexes = br_indexes[:-1]
        if len(br_indexes) >= 2:
            break
    reversal_start_index_i = np.random.randint(low=0, high=len(br_indexes) - 1)
    reversal_start_index = br_indexes[reversal_start_index_i]
    if reversal_start_index == 0:
        sp1 = None
        sp2 = chromosome[reversal_start_index].start_position
    else:
        sp1 = chromosome[reversal_start_index - 1].end_position
        sp2 = chromosome[reversal_start_index].start_position
    if reversal_start_index_i + 1 == len(br_indexes) - 1:
        reversal_end_index_i = len(br_indexes) - 1
    else:
        reversal_end_index_i = np.random.randint(low=reversal_start_index_i + 1, high=len(br_indexes))
    reversal_end_index = br_indexes[reversal_end_index_i]
    if reversal_end_index == len(chromosome):
        ep1 = chromosome[-1].end_position
        ep2 = None
    else:
        ep1 = chromosome[reversal_end_index - 1].end_position
        ep2 = chromosome[reversal_end_index].start_position
    return Mutation(mutation_type=MutationType.REVERSAL,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        REVERSAL_START_INDEX: reversal_start_index,
                        REVERSAL_END_INDEX: reversal_end_index,
                        POSITIONS: [sp1, sp2, ep1, ep2]
                    })


def generate_duplication(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = get_available_breakage_indexes(chromosome=chromosome, config=config)
        if not config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_TELOMERE]:
            if br_indexes[0] == 0:
                br_indexes = br_indexes[1:]
            if br_indexes[-1] == len(chromosome):
                br_indexes = br_indexes[:-1]
        if len(br_indexes) >= 2:
            break
    duplication_start_index_i = np.random.randint(low=0, high=len(br_indexes) - 1)
    duplication_start_index = br_indexes[duplication_start_index_i]
    if duplication_start_index == 0:
        sp1 = None
        sp2 = chromosome[duplication_start_index].start_position
    else:
        sp1 = chromosome[duplication_start_index - 1].end_position
        sp2 = chromosome[duplication_start_index].start_position
    if duplication_start_index_i + 1 == len(br_indexes) - 1:
        duplication_end_index_i = len(br_indexes) - 1
    else:
        duplication_end_index_i = np.random.randint(low=duplication_start_index_i + 1, high=len(br_indexes))
    duplication_end_index = br_indexes[duplication_end_index_i]
    if duplication_end_index == len(chromosome):
        ep1 = chromosome[-1].end_position
        ep2 = None
    else:
        ep1 = chromosome[duplication_end_index - 1].end_position
        ep2 = chromosome[duplication_end_index].start_position
    tandem_dup_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_TANDEM_PROBABILITY]
    tandem_duplication = np.random.choice(a=[True, False], p=[tandem_dup_prob, 1 - tandem_dup_prob])
    if tandem_duplication:
        insert_shift = 0
    else:
        # TODO: CHANGE, PLACEHOLDER
        insert_shift = 0
    reversed_dup_prob = config[MUTATIONS_TYPES_SPECIFICATIONS][MutationType.DUPLICATION][DUPLICATION_REVERSED_PROBABILITY]
    reversed_duplication = np.random.choice(a=[True, False], p=[reversed_dup_prob, 1 - reversed_dup_prob])
    return Mutation(mutation_type=MutationType.DUPLICATION,
                    mutation_data={
                        CHROMOSOME_INDEX: chromosome_index,
                        DUPLICATION_START_INDEX: duplication_start_index,
                        DUPLICATION_END_INDEX: duplication_end_index,
                        DUPLICATION_INSERT_SHIFT: insert_shift,
                        DUPLICATION_REVERSED: reversed_duplication,
                        POSITIONS: [sp1, sp2, ep1, ep2]
                    })


def generate_deletion(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = get_available_breakage_indexes(chromosome=chromosome, config=config)
        if not config["mut_types_spec"][MutationType.DELETION]["arm_deletion"]:
            if br_indexes[0] == 0:
                br_indexes = br_indexes[1:]
            if br_indexes[-1] == len(chromosome):
                br_indexes = br_indexes[:-1]
        if len(br_indexes) >= 2:
            break
    deletion_start_index_i = np.random.randint(low=0, high=len(br_indexes) - 1)
    deletion_start_index = br_indexes[deletion_start_index_i]
    if deletion_start_index == 0:
        sp1 = None
        sp2 = chromosome[deletion_start_index].start_position
    else:
        sp1 = chromosome[deletion_start_index - 1].end_position
        sp2 = chromosome[deletion_start_index].start_position
    deletion_end_index_i = np.random.randint(low=deletion_start_index_i + 1, high=len(br_indexes))
    deletion_end_index = br_indexes[deletion_end_index_i]
    if deletion_end_index == len(chromosome):
        ep1 = chromosome[-1].end_position
        ep2 = None
    else:
        ep1 = chromosome[deletion_end_index - 1].end_position
        ep2 = chromosome[deletion_end_index].start_position
    return Mutation(mutation_type=MutationType.DELETION,
                    mutation_data={
                        "chromosome_index": chromosome_index,
                        "deletion_start_index": deletion_start_index,
                        "deletion_end_index": deletion_end_index,
                        "positions": [sp1, sp2, ep1, ep2]
                    })


def chromosomes_are_mates(chromosome1, chromosome2):
    chr1_chr = [s.chromosome for s in chromosome1]
    chr2_chr = [s.chromosome for s in chromosome2]
    chr1_counter = Counter(chr1_chr)
    chr2_counter = Counter(chr2_chr)
    chr1_most_common_chr = chr1_counter.most_common(n=1)[0]
    chr2_most_common_chr = chr2_counter.most_common(n=1)[0]
    if chr1_most_common_chr != chr2_most_common_chr:
        return False


def generate_translocation(genome, config):
    while True:
        chromosome_1_index, chromosome_2_index = np.random.randint(low=0, high=len(genome), size=2)
        if chromosome_1_index == chromosome_2_index:
            continue
        chromosome_1 = genome[chromosome_1_index]
        chromosome_2 = genome[chromosome_2_index]
        if not config["mut_types_spec"][MutationType.TRANSLOCATION]["homozygous"] and \
                chromosomes_are_mates(chromosome1=chromosome_1, chromosome2=chromosome_2):
            continue
        chr1_br_indexes = get_available_breakage_indexes(chromosome=chromosome_1, config=config)
        if not config["mut_types_spec"][MutationType.TRANSLOCATION]["chr1_empty_arm"]:
            if chr1_br_indexes[0] == 0:
                chr1_br_indexes = chr1_br_indexes[1:]
            if chr1_br_indexes[-1] == len(chromosome_1):
                chr1_br_indexes = chr1_br_indexes[:-1]
        chr2_br_indexes = get_available_breakage_indexes(chromosome=chromosome_2, config=config)
        if not config["mut_types_spec"][MutationType.TRANSLOCATION]["chr2_empty_arm"]:
            if chr2_br_indexes[0] == 0:
                chr2_br_indexes = chr2_br_indexes[1:]
            if chr2_br_indexes[-1] == len(chromosome_2):
                chr2_br_indexes = chr2_br_indexes[:-1]
        if len(chr1_br_indexes) >= 1 and len(chr2_br_indexes) >= 1:
            break
    chr1_index = np.random.randint(low=0, high=len(chr1_br_indexes) - 1)
    chr2_index = np.random.randint(low=0, high=len(chr2_br_indexes) - 1)
    cc = np.random.choice([True, False])
    if chr1_index == 0:
        chr1_p1 = None
        chr1_p2 = chromosome_1[chr1_index].start_position
    elif chr1_index == len(chromosome_1):
        chr1_p1 = chromosome_1[chr1_index - 1].end_position
        chr1_p2 = None
    else:
        chr1_p1 = chromosome_1[chr1_index - 1].end_position
        chr1_p2 = chromosome_1[chr1_index].start_position
    if chr2_index == 0:
        chr2_p1 = None
        chr2_p2 = chromosome_2[chr2_index].start_position
    elif chr2_index == len(chromosome_2):
        chr2_p1 = chromosome_2[chr2_index - 1].end_position
        chr2_p2 = None
    else:
        chr2_p1 = chromosome_2[chr2_index - 1].end_position
        chr2_p2 = chromosome_2[chr2_index].start_position
    return Mutation(mutation_type=MutationType.TRANSLOCATION,
                    mutation_data={
                        "chr1_index": chromosome_1_index,
                        "chr2_index": chromosome_2_index,
                        "chr1_transl_index": chr1_index,
                        "chr2_transl_index": chr2_index,
                        "cc": cc,
                        "positions": [chr1_p1, chr1_p2, chr2_p1, chr2_p2]
                    })


def apply_mutation(genome, mutation):
    if mutation.mutation_type not in MUTATION_TYPES_to_MUTATING_FUNCTIONS:
        raise ValueError("Unsupported mutation type {mt}".format(mt=repr(mutation.mutation_type)))
    mutating_function = MUTATION_TYPES_to_MUTATING_FUNCTIONS[mutation.mutation_type]
    return mutating_function(genome=genome, mutation=mutation)


def apply_reversal(genome, mutation):
    new_genome = deepcopy(genome)
    chr_index = mutation.mutation_data["chromosome_index"]
    chromosome = new_genome[chr_index]
    rev_start_index = mutation.mutation_data["reversal_start_index"]
    rev_end_index = mutation.mutation_data["reversal_end_index"]
    new_chromosome = reverse_segments(chromosome=chromosome,
                                      start_segment_index=rev_start_index,
                                      end_segment_index=rev_end_index)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_duplication(genome, mutation):
    new_genome = deepcopy(genome)
    chr_index = mutation.mutation_data["chromosome_index"]
    chromosome = new_genome[chr_index]
    dup_start_index = mutation.mutation_data["duplication_start_index"]
    dup_end_index = mutation.mutation_data["duplication_end_index"]
    new_chromosome = duplicate_segments(chromosome=chromosome,
                                        start_segment_index=dup_start_index,
                                        end_segment_index=dup_end_index)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_deletion(genome, mutation):
    new_genome = deepcopy(genome)
    chr_index = mutation.mutation_data["chromosome_index"]
    chromosome = new_genome[chr_index]
    del_start_index = mutation.mutation_data["deletion_start_index"]
    del_end_index = mutation.mutation_data["deletion_end_index"]
    new_chromosome = delete_segments(chromosome=chromosome,
                                     start_segment_index=del_start_index,
                                     end_segment_index=del_end_index)
    new_genome[chr_index] = new_chromosome
    return new_genome


def apply_translocation(genome, mutation):
    new_genome = deepcopy(genome)
    chr1_index = mutation.mutation_data["chr1_index"]
    chr2_index = mutation.mutation_data["chr2_index"]
    chromosome_1 = new_genome[chr1_index]
    chromosome_2 = new_genome[chr2_index]
    chr1_trans_index = mutation.mutation_data["chr1_transl_index"]
    chr2_trans_index = mutation.mutation_data["chr2_transl_index"]
    cc = mutation.mutation_data["cc"]
    new_chr1, new_chr2 = translocate_segments(chromosome1=chromosome_1, chromosome2=chromosome_2,
                                              chromosome1_segment_index=chr1_trans_index,
                                              chromosome2_segment_index=chr2_trans_index, cc=cc)
    new_genome[chr1_index] = new_chr1
    new_genome[chr2_index] = new_chr2
    return new_genome


MUTATION_TYPES_to_GENERATING_FUNCTIONS = {
    MutationType.REVERSAL: generate_reversal,
    MutationType.DELETION: generate_deletion,
    MutationType.DUPLICATION: generate_duplication,
    MutationType.TRANSLOCATION: generate_translocation
}

MUTATION_TYPES_to_MUTATING_FUNCTIONS = {
    MutationType.REVERSAL: apply_reversal,
    MutationType.DELETION: apply_deletion,
    MutationType.DUPLICATION: apply_duplication,
    MutationType.TRANSLOCATION: apply_translocation
}


def generate_mutated_genome(starting_genome, mutation_cnt, mutations_config=MUTATION_CONFIG):
    if MUTATED_EXTREMITIES not in mutations_config:
        mutations_config[MUTATED_EXTREMITIES] = set()
    current_genome = deepcopy(starting_genome)
    history = {
        "genomes": [current_genome],
        "mutations": []
    }
    for _ in range(mutation_cnt):
        mut_types = mutations_config["mut_types"]
        mut_prob = mutations_config["mut_probs"]
        mutation_type = np.random.choice(a=mut_types,
                                         p=mut_prob)
        mutation = generate_mutation(mutation_type, current_genome, mutations_config)
        history["mutations"].append(mutation)
        current_genome = apply_mutation(current_genome, mutation)
        for p in mutation.mutation_data["positions"]:
            if p is not None:
                mutations_config[MUTATED_EXTREMITIES].add(p)
        history["genomes"].append(current_genome)
    return history


def get_adjacencies_from_genome(genome, is_reference=False, ref_adjacencies=None):
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
