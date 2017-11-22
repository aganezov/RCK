from collections import Counter
from copy import deepcopy

from dassp.simulation.parts import ChromosomeGenerator, MutationType, Mutation, reverse_segments, duplicate_segments, delete_segments
import numpy as np

from simulation.parts import translocation_segments

CHROMOSOMES_SIZE = 100
AB = True
ME = "mutated_extremities"

MUTATION_CONFIG = {
    "mut_types": [MutationType.REVERSAL,
                  MutationType.DUPLICATION,
                  MutationType.DELETION,
                  MutationType.TRANSLOCATION],
    "mut_probs": [0.15,
                  0.4,
                  0.4,
                  0.05],
    "mut_types_spec": {
        MutationType.DELETION: {
            "arm_deletion": False
        },
        MutationType.REVERSAL: {
            "arm_reversal": False
        },
        MutationType.TRANSLOCATION: {
            "chr1_empty_arm": False,
            "chr2_empty_arm": False,
            "homozygous": False
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
    def __init__(self, chrs_cnt, clones_cnt,
                 clonal_mut_cnt, subclonal_mut_cnt,
                 chrs_size=CHROMOSOMES_SIZE, ab=AB):
        super(ClonalSubClonalSimulationManager, self).__init__(chrs_cnt=chrs_cnt,
                                                               chrs_size=chrs_size,
                                                               clones_cnt=clones_cnt,
                                                               ab=ab)
        self.clonal_mutation_cnt = clonal_mut_cnt
        self.subclonal_mutation_cnt = subclonal_mut_cnt


def generate_mutation(mutation_type, genome, mutations_config):
    if mutation_type not in MUTATION_TYPES_to_GENERATING_FUNCTIONS:
        raise ValueError("Unsupported mutation type {mt}".format(mt=repr(mutation_type)))
    generating_function = MUTATION_TYPES_to_GENERATING_FUNCTIONS[mutation_type]
    return generating_function(genome=genome, config=mutations_config)


def breakage_indexes(chromosome, config):
    result = []
    if len(chromosome) == 0:
        return result
    if chromosome[0].start_position not in config[ME]:
        result.append(0)
    for i in range(len(chromosome) - 1):
        p1 = chromosome[i].end_position
        p2 = chromosome[i + 1].start_position
        if p1 not in config[ME] and p2 not in config[ME]:
            result.append(i + 1)
    if chromosome[-1].end_position not in config[ME]:
        result.append(len(chromosome))
    return result


def generate_reversal(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = breakage_indexes(chromosome=chromosome, config=config)
        if not config["mut_types_spec"][MutationType.REVERSAL]["arm_reversal"]:
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
                        "chromosome_index": chromosome_index,
                        "reversal_start_index": reversal_start_index,
                        "reversal_end_index": reversal_end_index,
                        "positions": [sp1, sp2, ep1, ep2]
                    })


def generate_duplication(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = breakage_indexes(chromosome=chromosome, config=config)
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
    duplication_end_index_i = np.random.randint(low=duplication_start_index_i + 1, high=len(br_indexes))
    duplication_end_index = br_indexes[duplication_end_index_i]
    if duplication_end_index == len(chromosome):
        ep1 = chromosome[-1].end_position
        ep2 = None
    else:
        ep1 = chromosome[duplication_end_index - 1].end_position
        ep2 = chromosome[duplication_end_index].start_position
    return Mutation(mutation_type=MutationType.DUPLICATION,
                    mutation_data={
                        "chromosome_index": chromosome_index,
                        "duplication_start_index": duplication_start_index,
                        "duplication_end_index": duplication_end_index,
                        "positions": [sp1, sp2, ep1, ep2]
                    })


def generate_deletion(genome, config):
    while True:
        chromosome_index = np.random.randint(low=0, high=len(genome))
        chromosome = genome[chromosome_index]
        br_indexes = breakage_indexes(chromosome=chromosome, config=config)
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
        chr1_br_indexes = breakage_indexes(chromosome=chromosome_1, config=config)
        if not config["mut_types_spec"][MutationType.TRANSLOCATION]["chr1_empty_arm"]:
            if chr1_br_indexes[0] == 0:
                chr1_br_indexes = chr1_br_indexes[1:]
            if chr1_br_indexes[-1] == len(chromosome_1):
                chr1_br_indexes = chr1_br_indexes[:-1]
        chr2_br_indexes = breakage_indexes(chromosome=chromosome_2, config=config)
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
    new_chr1, new_chr2 = translocation_segments(chromosome1=chromosome_1, chromosome2=chromosome_2,
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
    if ME not in mutations_config:
        mutations_config[ME] = set()
    current_genome = deepcopy(starting_genome)
    history = {
        "genomes": [current_genome],
        "mutations": []
    }
    for _ in range(mutation_cnt):
        mutation_type = np.random.choice(a=mutations_config["mut_types"],
                                         p=mutations_config.get("mut_prob"))
        mutation = generate_mutation(mutation_type, current_genome, mutations_config)
        history["mutations"].append(mutation)
        current_genome = apply_mutation(current_genome, mutation)
        for p in mutation.mutation_data["positions"]:
            if p is not None:
                mutations_config[ME].add(p)
        history["genomes"].append(current_genome)
    return history
