import os
import sys
import pprint
import numpy as np
from copy import deepcopy
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


from dassp.core.graph import construct_hiag_inflate_from_haploid_data
from dassp.core.structures import get_unique_haploid_telomeres, Phasing, get_novel_adjacencies_from_ref_and_mut_genomes, get_shuffled_clone_specific_scnp, \
    segment_copy_number_tensors_are_compatible
from dassp.algo.ilp import SatModelSingleClone, SatModelMultiClone
from dassp.core.structures import SegmentCopyNumberProfile, get_unique_forward_haploid_segments, get_shuffled_scnp, get_unique_haploid_ref_and_novel_sorted_adjacencies
from dassp.core.structures import StructureProfile
from dassp.core.structures import get_telomeres_from_genome, strip_haplotype_from_positions
from dassp.core.graph import construct_iag, construct_hiag
from dassp.simulation.manager import generate_mutated_genome, MUTATION_CONFIG, HIIS, HSIS, Manager, assign_uniform_mutation_configs, GENOME
from dassp.simulation.parts import ChromosomeGenerator
from dassp.simulation.io import write_genome_to_file, write_mutation_history_to_file, IO_CONFIG
from dassp.simulation.io import read_genomes_from_file, write_phylogenetic_mutation_history_to_file
from dassp.simulation.manager import generate_random_phylogeny_tree, assign_reference_genome, tree_cnt
from dassp.simulation.io import read_phylogenetic_mutation_history_from_file


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
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description="Multi-clone ILP feasibility test on simulated datasets")
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--end', type=int, default=2000)
    parser.add_argument('--dataset-dir', default=os.path.join(current_dir, "instances", "hiis"))
    parser.add_argument('--dataset-file-template', type=str, default="{cnt}_mut_history.txt")
    parser.add_argument('--no-gurobi-silent', dest="gurobi_silent", action="store_false")
    parser.add_argument('--clone-cnt', type=int, choices=[2, 3], default=3)
    args = parser.parse_args()

    for cnt in range(args.start, args.end):
        print("\r{cnt}".format(cnt=cnt), flush=True, end="")
        file_name = args.dataset_file_template.format(cnt=cnt)
        full_file_name = os.path.join(args.dataset_dir, file_name)
        manager = read_phylogenetic_mutation_history_from_file(file_name=full_file_name)
        manager.propagate_haplotypes_to_positions()
        ref_genome = manager.tree.nodes[manager.root][GENOME]

        hapl_segments = get_unique_forward_haploid_segments(genome=ref_genome)

        mut_genomes = {}
        mgc = 0
        for node in manager.tree:
            if node == manager.root:
                continue
            mut_genomes[node] = manager.tree.nodes[node][GENOME]
        mut_genomes_keys = np.random.choice(a=list(mut_genomes.keys()), size=args.clone_cnt, replace=False)
        mut_genomes = {key: mut_genomes[key] for key in mut_genomes_keys}

        mut_clone_specific_scnp = {}
        for mut_genome_name, mut_genome in mut_genomes.items():
            mut_clone_specific_scnp[mut_genome_name] = SegmentCopyNumberProfile.from_genome(genome=mut_genome)
        shuffled_mut_clone_specific_scnp = get_shuffled_clone_specific_scnp(clone_specific_scnp=mut_clone_specific_scnp, segments=hapl_segments)

        dipl_ref_telomeres = list(ref_genome.iter_telomeres())
        dipl_mut_telomeres = []
        for mut_genome in mut_genomes.values():
            dipl_mut_telomeres.extend(list(mut_genome.iter_telomeres()))
        hapl_telomeres = get_unique_haploid_telomeres(genome=ref_genome)
        assert len(set(dipl_ref_telomeres)) == len(set(dipl_mut_telomeres))
        for telomere in dipl_mut_telomeres:
            assert telomere in dipl_ref_telomeres

        ref_hapl_adjacencies, nov_hapl_adjacencies = get_unique_haploid_ref_and_novel_sorted_adjacencies(ref_genome=ref_genome, mut_genomes=list(mut_genomes.values()))
        hapl_adjacencies = ref_hapl_adjacencies + nov_hapl_adjacencies

        model = SatModelMultiClone(hapl_segments=hapl_segments, hapl_adjacencies=hapl_adjacencies, clone_specific_scnp=shuffled_mut_clone_specific_scnp, hapl_telomeres=hapl_telomeres)
        model.build_gurobi_model()
        if args.gurobi_silent:
            model.gm.setParam("OutputFlag", False)
        model.solve_model()
        inferred_clone_specific_scnp = model.get_clone_specific_scnp_from_model()
        inferred_clone_specific_acnp = model.get_clone_specific_acnp_from_models()

        # checking that reference adjacencies don't have inconsistent (i.e., AB/BA phasing values)
        for clone_id in mut_genomes_keys:
            for adj in ref_hapl_adjacencies:
                for ph in [Phasing.AB, Phasing.BA]:
                    assert inferred_clone_specific_acnp[clone_id].get_phase_aware_cn_by_adj_and_phasing(adjacency=adj, phasing=ph) == 0

        # checking that for every haploid novel adjacency exactly one underlying diploid novel adjacency has a positive values (across all clones)
        for adj in nov_hapl_adjacencies:
            cns = []
            present_ph = 0
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                inferred_cns = [inferred_clone_specific_acnp[clone_id].get_phase_aware_cn_by_adj_and_phasing(adjacency=adj, phasing=ph) for clone_id in mut_genomes_keys]
                if any(map(lambda entry: entry >= 1, inferred_cns)):
                    present_ph += 1
                    cns = inferred_cns
            assert present_ph == 1
            assert sum(cns) > 0

        # checking that self-loop adjacencies don't have AB/BA phasing
        for adj in filter(lambda a: a.is_self_loop_hapl, nov_hapl_adjacencies):
            for clone_id in mut_genomes_keys:
                for ph in [Phasing.AB, Phasing.BA]:
                    assert inferred_clone_specific_acnp[clone_id].get_phase_aware_cn_by_adj_and_phasing(adjacency=adj, phasing=ph) == 0

        # checking that inferred clone specific segment copy number tensor is compatible with the real diploid segment copy number tensor
        assert segment_copy_number_tensors_are_compatible(tensor1=mut_clone_specific_scnp, tensor2=inferred_clone_specific_scnp, segments=hapl_segments)

        # checking if every non-telomere vertex is balanced and every telomere vertex is unbalanced.
        for clone_id in mut_genomes_keys:
            inflated_hiag = construct_hiag_inflate_from_haploid_data(hapl_segments=hapl_segments, hapl_adjacencies=hapl_adjacencies)
            inflated_hiag.assign_copy_numbers_from_scn_profile(scn_profile=inferred_clone_specific_scnp[clone_id])
            inflated_hiag.assign_copy_numbers_from_acn_profile(acn_profile=inferred_clone_specific_acnp[clone_id])
            assert inflated_hiag.represents_a_genome
            inferred_telomeres = inflated_hiag.get_telomeres(check_cn_awareness=False)
            assert len(set(inferred_telomeres)) == len(set(dipl_ref_telomeres))
            for telomere in inferred_telomeres:
                assert telomere in dipl_ref_telomeres



    # sys.setrecursionlimit(10000)
    # for cnt in range(1800, 2000):
    #     print(cnt, end=" ", flush=True)
    #     ref_genome = ChromosomeGenerator.generate_genome(chromosome_size=300,
    #                                                      chromosomes_cnt=5,
    #                                                      ab=True)
    #     tree, root = generate_random_phylogeny_tree(n=3,
    #                                                 return_root=True)
    #     manager = Manager(phylo_tree=tree)
    #     assign_uniform_mutation_configs(manager=manager,
    #                                     mut_cnt=50,
    #                                     mut_config=MUTATION_CONFIG)
    #     assign_reference_genome(manager=manager,
    #                             ref_genome=ref_genome)
    #     manager.generate_mutated_genomes()
    #     write_phylogenetic_mutation_history_to_file(manager=manager, file_name="dassp/simulation/instances/nr/{cnt}_mut_history.txt".format(cnt=cnt))
    #     # manager = read_phylogenetic_mutation_history_from_file(file_name="dassp/simulation/instances/hiis/1000_mut_history.txt")
    #     manager.propagate_haplotypes_to_positions()
    #     ref_genome = manager.tree.nodes[manager.root][GENOME]
    #     mut_genomes = []
    #     for node in manager.tree:
    #         if node == manager.root:
    #             continue
    #         mut_genomes.append(manager.tree.nodes[node][GENOME])
    #     keep = False
    #     for node in manager.tree:
    #         if node == manager.root:
    #             continue
    #         mut_genome = manager.tree.nodes[node][GENOME]
    #         mut_genome_name = str(node)
    #         structure_profile = StructureProfile.from_genome(genome=mut_genome)
    #         iag = construct_iag(ref_genome=ref_genome, mut_genomes=[mut_genome])
    #         if not iag.complies_with_is:
    #             keep = False
    #             print("IAG for {gn} violates HIIS".format(gn=mut_genome_name), end="; ", flush=True)
    #         if not iag.compatible_with_hsis:
    #             keep = False
    #             print("IAG for {gn} not compatible with HSIS".format(gn=mut_genome_name), end="; ", flush=True)
    #         if not iag.topology_allows_for_genome(genome=mut_genome):
    #             keep = True
    #             print("IAG for {gn} does not allow for the mutated_genome".format(gn=mut_genome_name), end="; ", flush=True)
    #         iag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=False)
    #         if not iag.represents_given_genome(genome=mut_genome):
    #             keep = True
    #             print("IAG for {gn} does not represent the mutated genome".format(gn=mut_genome_name), end="; ", flush=True)
    #         genome_telomeres = [t.get_non_hap_copy() for t in ref_genome.iter_telomeres()]
    #         graph_telomeres = [t for t in iag.iter_telomeres(check_cn_awareness=False)]
    #         genome_telomeres = set(genome_telomeres)
    #         graph_telomeres = set(graph_telomeres)
    #         sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
    #         if len(sym_dif) > 0:
    #             keep = True
    #             print("Not all telomeres are preserved in IAG for {gn}".format(gn=mut_genome_name), end="; ", flush=True)
    #             # print(sorted(sym_dif))
    #         hsiag = construct_hiag(ref_genome=ref_genome, mut_genomes=[mut_genome])
    #         if not hsiag.complies_with_is:
    #             keep = False
    #             print("HSIAG for {gn} violates HSIS".format(gn=mut_genome_name), end="; ", flush=True)
    #         if not hsiag.complies_with_hsis:
    #             keep = False
    #             print("HSIAG for {gn} violates HSIS".format(gn=mut_genome_name), end="; ", flush=True)
    #         if not hsiag.complies_with_hiis:
    #             keep = False
    #             print("HSIAG for {gn} violates HIIS".format(gn=mut_genome_name), end="; ", flush=True)
    #         if not hsiag.topology_allows_for_genome(genome=mut_genome):
    #             keep = True
    #             print("HSIAG for {gn} does not allow for the mutated genome".format(gn=mut_genome_name), end="; ", flush=True)
    #         hsiag.assign_copy_numbers_from_genome(genome=mut_genome, ensure_topology=False)
    #         if not hsiag.represents_given_genome(genome=mut_genome):
    #             keep = True
    #             print("HSIAG for {gn} does not represent the mutated genome".format(gn=mut_genome_name), end="; ", flush=True)
    #         genome_telomeres = [deepcopy(t) for t in ref_genome.iter_telomeres()]
    #         graph_telomeres = [t for t in hsiag.iter_telomeres(check_cn_awareness=False)]
    #         genome_telomeres = set(genome_telomeres)
    #         graph_telomeres = set(graph_telomeres)
    #         sym_dif = genome_telomeres.symmetric_difference(graph_telomeres)
    #         if len(sym_dif) > 0:
    #             keep = True
    #             print("Not all telomeres are preserved in HSIAG for {gn}".format(gn=mut_genome_name), end="; ", flush=True)
    #             # print(sorted(sym_dif))
    #     iag = construct_iag(ref_genome=ref_genome, mut_genomes=mut_genomes)
    #     if not iag.complies_with_is:
    #         keep = False
    #         print("Violates HIIS", end="; ", flush=True)
    #     if not iag.compatible_with_hsis:
    #         keep = False
    #         print("Not compatible with HSIS", end="; ", flush=True)
    #     for mut_genome in mut_genomes:
    #         if not iag.topology_allows_for_genome(genome=mut_genome):
    #             keep = True
    #             print("Does not allow for the mutated_genome", end="; ", flush=True)
    #     hsiag = construct_hiag(ref_genome=ref_genome, mut_genomes=mut_genomes)
    #     if not hsiag.complies_with_is:
    #         keep = False
    #         print("Violates HSIS", end="; ", flush=True)
    #     if not hsiag.complies_with_hsis:
    #         keep = False
    #         print("Violates HSIS", end="; ", flush=True)
    #     if not hsiag.complies_with_hiis:
    #         keep = False
    #         print("Violates HIIS", end="; ", flush=True)
    #     for mut_genome in mut_genomes:
    #         if not hsiag.topology_allows_for_genome(genome=mut_genome):
    #             keep = True
    #             print("Does not allow for the mutated genome", end="; ", flush=True)
    #     if keep:
    #         print(flush=True)
    #     print("\r", end="", flush=True)

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
