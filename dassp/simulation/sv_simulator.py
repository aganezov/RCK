import os
import pprint
from copy import deepcopy
import networkx as nx

from dassp.core.structures import get_telomeres_from_genome, strip_haplotype_from_positions
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
    initial_genome = ChromosomeGenerator.generate_genome(chromosome_size=chrs_size,
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
    manager = read_phylogenetic_mutation_history_from_file(file_name="dassp/simulation/instances/hiis/5_mut_history.txt")
    manager.propagate_haplotypes_to_positions()
    ref_genome = manager.tree.nodes[manager.root][GENOME]
    mut_genomes = []
    for node in manager.tree:
        if node == manager.root:
            continue
        mut_genomes.append(manager.tree.nodes[node][GENOME])
    for mut_genome in mut_genomes:
        iag = construct_iag(ref_genome=ref_genome, mut_genomes=[mut_genome])
        if not iag.complies_with_is:
            print("Violates HIIS")
        if not iag.compatible_with_hsis:
            print("Not compatible with HSIS")
        if not iag.topology_allows_for_genome(genome=mut_genome):
            print("Does not allow for the mutated_genome")
        iag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=False)
        if not iag.represents_given_genome(genome=mut_genome):
            print("Does not represent the mutated genome")
        genome_telomeres = [t.get_non_hap_copy() for t in ref_genome.iter_telomeres()]
        graph_telomeres = [t for t in iag.iter_telomeres(check_cn_awareness=False)]
        genome_telomeres = set(genome_telomeres)
        graph_telomeres = set(graph_telomeres)
        sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
        if len(sym_dif) > 0:
            print("Not all telomeres are preserved")
            print(sorted(sym_dif))
        hsiag = construct_hiag(ref_genome=ref_genome, mut_genomes=[mut_genome])
        if not hsiag.complies_with_is:
            print("Violates HSIS")
        if not hsiag.complies_with_hsis:
            print("Violates HSIS")
        if not hsiag.complies_with_hiis:
            print("Violates HIIS")
        if not hsiag.topology_allows_for_genome(genome=mut_genome):
            print("Does not allow for the mutated genome")
        hsiag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=False)
        if not hsiag.represents_given_genome(genome=mut_genome):
            print("Does not represent the mutated genome")
        genome_telomeres = [deepcopy(t) for t in ref_genome.iter_telomeres()]
        graph_telomeres = [t for t in hsiag.iter_telomeres(check_cn_awareness=False)]
        genome_telomeres = set(genome_telomeres)
        graph_telomeres = set(graph_telomeres)
        sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
        if len(sym_dif) > 0:
            print("Not all telomeres are preserved")
            print(sorted(sym_dif))
    iag = construct_iag(ref_genome=ref_genome, mut_genomes=mut_genomes)
    if not iag.complies_with_is:
        print("Violates HIIS")
    if not iag.compatible_with_hsis:
        print("Not compatible with HSIS")
    for mut_genome in mut_genomes:
        if not iag.topology_allows_for_genome(genome=mut_genome):
            print("Does not allow for the mutated_genome")
    hsiag = construct_hiag(ref_genome=ref_genome, mut_genomes=mut_genomes)
    if not hsiag.complies_with_is:
        print("Violates HSIS")
    if not hsiag.complies_with_hsis:
        print("Violates HSIS")
    if not hsiag.complies_with_hiis:
        print("Violates HIIS")
    for mut_genome in mut_genomes:
        if not hsiag.topology_allows_for_genome(genome=mut_genome):
            print("Does not allow for the mutated genome")

    # ref_genome = ChromosomeGenerator.generate_genome(chromosome_size=300,
    #                                                  chromosomes_cnt=5,
    #                                                  ab=True)
    # b = ref_genome[2:5]
    # a = 5
    # for cnt in range(750, 1000):
    #     print(cnt, end=" ", flush=True)
    #     ref_genome = ChromosomeGenerator.generate_chromosomes(chromosome_size=300,
    #                                                           chromosomes_cnt=5,
    #                                                           ab=True)
    #     tree, root = generate_random_phylogeny_tree(n=3,
    #                                                 return_root=True)
    #     manager = Manager(phylo_tree=tree)
    #     assign_uniform_mutation_configs(manager=manager,
    #                                     mut_cnt=50,
    #                                     mut_config=MUTATION_CONFIG)
    #     assign_reference_genome(manager=manager,
    #                             ref_genome=ref_genome)
    #     manager.generate_mutated_genomes()
    #     write_phylogenetic_mutation_history_to_file(manager=manager, file_name="dassp/simulation/instances/{cnt}_mut_history.txt".format(cnt=cnt))
    #     manager.propagate_haplotypes_to_positions()
    #     ref_genome = manager.tree.nodes[manager.root][GENOME]
    #     mut_genomes = []
    #     for node in manager.tree:
    #         if node == manager.root:
    #             continue
    #         mut_genomes.append(manager.tree.nodes[node][GENOME])
    #     for mut_genome in mut_genomes:
    #         iag = construct_iag(ref_genome=ref_genome,
    #                             mut_genomes=[mut_genome],
    #                             build_graph=True)
    #         if not iag.complies_with_is:
    #             print("IAG IS violation", end=" ", flush=True)
    #         assert iag.topology_allows_for_genome(genome=mut_genome)
    #         iag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=False, inherit_segment_topology=False, inherit_adjacency_topology=False)
    #         assert iag.represents_given_genome(genome=mut_genome)
    #         genome_telomeres = get_telomeres_from_genome(genome=mut_genome, copy=True, inherit_haplotypes=True)
    #         strip_haplotype_from_positions(positions=genome_telomeres, inplace=True)
    #         genome_telomeres = set(genome_telomeres)
    #         graph_telomeres = iag.get_telomeres(check_cn_awareness=False, sort=True, copy=False)
    #         graph_telomeres = set(graph_telomeres)
    #         sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
    #         # print(sorted(sym_dif))
    #         assert len(sym_dif) == 0
    #         hsiag = construct_hiag(ref_genome=ref_genome,
    #                                mut_genomes=[mut_genome],
    #                                build_graph=True)
    #         if not hsiag.complies_with_is:
    #             print("HSIAG IS violation", end=" ", flush=True)
    #         if not hsiag.complies_with_hiis:
    #             print("HSIAG HIIS violation", end=" ", flush=True)
    #         assert hsiag.complies_with_hsis
    #         assert hsiag.topology_allows_for_genome(genome=mut_genome)
    #         hsiag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=True, inherit_segment_topology=False, inherit_adjacency_topology=False)
    #         assert hsiag.represents_given_genome(genome=mut_genome)
    #         genome_telomeres = get_telomeres_from_genome(genome=mut_genome, copy=True, inherit_haplotypes=True)
    #         genome_telomeres = set(genome_telomeres)
    #         graph_telomeres = hsiag.get_telomeres(check_cn_awareness=False, sort=True, copy=False)
    #         graph_telomeres = set(graph_telomeres)
    #         sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
    #         # print(sorted(sym_dif))
    #         assert len(sym_dif) == 0
    #     assert not iag.represents_given_genome(genome=mut_genomes[0])
    #     assert not hsiag.represents_given_genome(genome=mut_genomes[0])
    #     iag = construct_iag(ref_genome=ref_genome,
    #                         mut_genomes=mut_genomes,
    #                         build_graph=True)
    #     if not iag.complies_with_is:
    #         print("IAG IS violation (all genomes NAs)", end=" ", flush=True)
    #     hsiag = construct_hiag(ref_genome=ref_genome,
    #                            mut_genomes=mut_genomes,
    #                            build_graph=True)
    #     if not hsiag.complies_with_is:
    #         print("HSIAG IS violation (all genomes NAs)", end=" ", flush=True)
    #     if not hsiag.complies_with_hiis:
    #         print("HSIAG HIIS violation (all genomes NAs)", end=" ", flush=True)
    #     assert hsiag.complies_with_hsis
    #     print("\r", end="", flush=True)
    # print(pprint.pformat(MUTATION_CONFIG))
    # manager = read_phylogenetic_mutation_history_from_file(file_name="dassp/simulation/instances/609_mut_history.txt")
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
    #     # assert iag.complies_with_is
    #     iag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=True, inherit_segment_topology=False, inherit_adjacency_topology=False)
    #     genome_telomeres = get_telomeres_from_genome(genome=mut_genome, copy=True, inherit_haplotypes=True)
    #     strip_haplotype_from_positions(positions=genome_telomeres, inplace=True)
    #     genome_telomeres = set(genome_telomeres)
    #     graph_telomeres = iag.get_telomeres(check_cn_awareness=False, sort=True, copy=False)
    #     graph_telomeres = set(graph_telomeres)
    #     sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
    #     # print(sorted(sym_dif))
    #     assert len(sym_dif) == 0
    #     hsiag = construct_hiag(ref_genome=ref_genome,
    #                            mut_genomes=[mut_genome],
    #                            build_graph=True)
    #     # assert hsiag.complies_with_is
    #     # assert hsiag.complies_with_hiis
    #     assert hsiag.complies_with_hsis
    #     assert hsiag.topology_allows_for_genome(genome=mut_genome)
    #
    # iag = construct_iag(ref_genome=ref_genome,
    #                     mut_genomes=mut_genomes,
    #                     build_graph=True)
    # # assert iag.complies_with_is
    # hsiag = construct_hiag(ref_genome=ref_genome,
    #                        mut_genomes=mut_genomes,
    #                        build_graph=True)
    # # assert hsiag.complies_with_is
    # # assert hsiag.complies_with_hiis
    # assert hsiag.complies_with_hsis
