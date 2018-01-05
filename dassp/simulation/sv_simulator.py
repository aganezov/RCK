import os
from copy import deepcopy
import networkx as nx

from dassp.core.graph import construct_iag, construct_hiag
from dassp.simulation.manager import generate_mutated_genome, MUTATION_CONFIG, HIIS, HSIS, Manager, assign_uniform_mutation_configs, GENOME
from dassp.simulation.parts import ChromosomeGenerator
from dassp.simulation.io import write_genome_to_file, write_mutation_history_to_file, IO_CONFIG
from dassp.simulation.io import read_genomes_from_file, write_phylogenetic_mutation_history_to_file
from dassp.simulation.manager import generate_random_phylogeny_tree, assign_reference_genome, tree_cnt
from simulation.io import read_phylogenetic_mutation_history_from_file


def genome_string_name_for_file(prefix, ab, mutation_cnt, config):
    result = [str(prefix), str(mutation_cnt), "AB" if ab else "A", "HIIS" if config[HIIS] else ("HSIS" if config[HSIS] else "nIS")]
    return "_".join(result) + ".txt"


def main(chrs_cnt=5, chrs_size=300, ab=True, mutations_cnt=60, sim_genomes_cnt=100, path_prefix="dassp/simulation/sim_5_300_ab_HIIS_60", mut_config=None, io_config=None):
    if mut_config is None:
        mut_config = MUTATION_CONFIG
    if io_config is None:
        io_config = IO_CONFIG
    initial_genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=chrs_size,
                                                              chromosomes_cnt=chrs_cnt,
                                                              ab=ab)
    write_genome_to_file(genome=initial_genome, file_name=os.path.join(path_prefix, "reference.txt"), config=io_config)
    for g_cnt in range(sim_genomes_cnt):
        mutation_history = generate_mutated_genome(starting_genome=deepcopy(initial_genome),
                                                   mutation_cnt=mutations_cnt,
                                                   config=deepcopy(mut_config))
        write_mutation_history_to_file(history=mutation_history,
                                       file_name=os.path.join(path_prefix, genome_string_name_for_file(prefix=g_cnt,
                                                                                                       ab=ab,
                                                                                                       mutation_cnt=mutations_cnt,
                                                                                                       config=mut_config)),
                                       config=io_config)


if __name__ == "__main__":
    # main()
    import sys
    sys.setrecursionlimit(10000)
    # manager = read_phylogenetic_mutation_history_from_file(file_name="dassp/simulation/instances/0_mut_history.txt")
    # a = 5
    for cnt in range(85, 100):
        print(cnt, end="", flush=True)
        ref_genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=300,
                                                              chromosomes_cnt=5,
                                                              ab=True)
        tree, root = generate_random_phylogeny_tree(n=3,
                                                    return_root=True)
        manager = Manager(phylo_tree=tree)
        assign_uniform_mutation_configs(manager=manager,
                                        mut_cnt=50,
                                        mut_config=MUTATION_CONFIG)
        assign_reference_genome(manager=manager,
                                ref_genome=ref_genome)
        manager.generate_mutated_genomes()
        write_phylogenetic_mutation_history_to_file(manager=manager, file_name="dassp/simulation/instances/{cnt}_mut_history.txt".format(cnt=cnt))
        manager.propagate_haplotypes_to_positions()
        ref_genome = manager.tree.nodes[manager.root][GENOME]
        mut_genomes = []
        for node in manager.tree:
            if node == manager.root:
                continue
            mut_genomes.append(manager.tree.nodes[node][GENOME])
        for mut_genome in mut_genomes:
            iag = construct_iag(ref_genome=ref_genome,
                                mut_genomes=[mut_genome],
                                build_graph=True)
            assert iag.complies_with_is
            hsiag = construct_hiag(ref_genome=ref_genome,
                                   mut_genomes=[mut_genome],
                                   build_graph=True)
            assert hsiag.complies_with_is
            assert hsiag.complies_with_hiis
            assert hsiag.complies_with_hsis
        iag = construct_iag(ref_genome=ref_genome,
                            mut_genomes=mut_genomes,
                            build_graph=True)
        assert iag.complies_with_is
        hsiag = construct_hiag(ref_genome=ref_genome,
                               mut_genomes=mut_genomes,
                               build_graph=True)
        assert hsiag.complies_with_is
        assert hsiag.complies_with_hiis
        assert hsiag.complies_with_hsis
        print("\r", end="", flush=True)

    # manager = read_phylogenetic_mutation_history_from_file(file_name="dassp/simulation/instances/83_mut_history.txt")
    # manager.propagate_haplotypes_to_positions()
    # ref_genome = manager.tree.nodes[manager.root][GENOME]
    # mut_genomes = []
    # for node in manager.tree:
    #     if node == manager.root:
    #         continue
    #     mut_genomes.append(manager.tree.nodes[node][GENOME])
    # for mut_genome in mut_genomes:
    #     iag = construct_iag(ref_genome=ref_genome,
    #                         mut_genomes=[mut_genome],
    #                         build_graph=True)
    #     assert iag.complies_with_is
    #     hsiag = construct_hiag(ref_genome=ref_genome,
    #                            mut_genomes=[mut_genome],
    #                            build_graph=True)
    #     assert hsiag.complies_with_is
    #     assert hsiag.complies_with_hiis
    #     assert hsiag.complies_with_hsis
    # iag = construct_iag(ref_genome=ref_genome,
    #                     mut_genomes=mut_genomes,
    #                     build_graph=True)
    # assert iag.complies_with_is
    # hsiag = construct_hiag(ref_genome=ref_genome,
    #                        mut_genomes=mut_genomes,
    #                        build_graph=True)
    # assert hsiag.complies_with_is
    # assert hsiag.complies_with_hiis
    # assert hsiag.complies_with_hsis