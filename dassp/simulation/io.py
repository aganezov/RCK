import pprint
import itertools

import networkx as nx

from dassp.core.structures import Segment, Genome, Chromosome
from dassp.simulation.manager import GENOMES, MUTATIONS, MUT_CONFIG, GENOME, HISTORY, Manager

SEGMENTS_SEPARATOR = "segment_separator"
CHROMOSOME_PREFIX = "chromosome_prefix"
CHROMOSOME_SUFFIX = "chromosome_suffix"
GENOME_HEADER = "genome_header"
GENOME_DATA_HEADER = "genome_data_header"
GENOME_FOOTER = "genome_footer"
GENOME_DATA_FOOTER = "genome_data_footer"
CHROMOSOMES_SEPARATOR = "chromosome_separator"
COMMENT_INDICATOR = "comment_indicator"
MUTATION_HEADER = "mutation_header"
MUTATION_FOOTER = "mutation_footer"
MUTATIONS_SEPARATOR = "mutation_separator"
TREE_PREFIX = "tree_prefix"
TREE_DATA_PREFIX = "tree_data_prefix"
TREE_SUFFIX = "tree_suffix"
TREE_DATA_SUFFIX = "tree_data_suffix"
TREE_LINE_SEPARATOR = "tree_line_separator"
TREE_ADJ_LIST_SEPARATOR = "tree_adj_list_separator"
EDGE_HEADER = "edge_header"
EDGE_DATA_HEADER = "edge_data_header"
EDGE_FOOTER = "edge_footer"
EDGE_DATA_FOOTER = "edge_data_footer"
MUT_CONFIG_HEADER = "mut_config_header"
MUT_CONFIG_FOOTER = "mut_config_footer"

IO_CONFIG = {
    SEGMENTS_SEPARATOR: "\t",
    CHROMOSOME_PREFIX: "",
    CHROMOSOME_SUFFIX: "",
    CHROMOSOMES_SEPARATOR: "\n",
    GENOME_HEADER: ">",
    GENOME_DATA_HEADER: ">>",
    GENOME_FOOTER: "<",
    GENOME_DATA_FOOTER: "<<",
    COMMENT_INDICATOR: "#",
    EDGE_HEADER: "\\",
    EDGE_DATA_HEADER: "\\\\",
    EDGE_FOOTER: "/",
    EDGE_DATA_FOOTER: "//",
    MUTATION_HEADER: "",
    MUTATION_FOOTER: "",
    MUTATIONS_SEPARATOR: "\n",
    TREE_PREFIX: "[",
    TREE_DATA_PREFIX: "[[",
    TREE_DATA_SUFFIX: "]]",
    TREE_SUFFIX: "]",
    TREE_LINE_SEPARATOR: "\n",
    TREE_ADJ_LIST_SEPARATOR: " ",
    MUT_CONFIG_HEADER: "\n",
    MUT_CONFIG_FOOTER: "\n",
}


class ParseError(Exception):
    pass


def comment_aware_string(string, commented, config, multistring_aware=True):
    if commented:
        if multistring_aware:
            return "\n".join(config[COMMENT_INDICATOR] + " " + entry for entry in string.split("\n"))
        else:
            return config[COMMENT_INDICATOR] + " " + string
    return string


def chromosome_to_string(chromosome, config, commented=False):
    result = [comment_aware_string(string=config[CHROMOSOME_PREFIX], commented=commented, config=config)]
    for segment in chromosome:
        result.append(comment_aware_string(string=str(segment), commented=commented, config=config))
    return config[SEGMENTS_SEPARATOR].join(result)


def genome_to_string(genome, config, commented=False):
    result = []
    for chromosome in genome:
        result.append(comment_aware_string(string=chromosome_to_string(chromosome=chromosome, config=config), commented=commented, config=config))
    return config[CHROMOSOMES_SEPARATOR].join(result)


def mutation_to_string(mutation, config, commented=False):
    result = [comment_aware_string(string=config[MUTATION_HEADER], commented=commented, config=config, multistring_aware=True),
              comment_aware_string(string=str(mutation.mutation_type), commented=commented, config=config),
              comment_aware_string(string=pprint.pformat(object=mutation.mutation_data), commented=commented, config=config, multistring_aware=True),
              comment_aware_string(string=config[MUTATION_FOOTER], commented=commented, config=config)]
    return config[MUTATIONS_SEPARATOR].join(result)


def tree_to_string(tree, config, commented=False):
    result = [comment_aware_string(string=config[TREE_PREFIX], config=config, commented=commented),
              comment_aware_string(string=config[TREE_DATA_PREFIX], commented=commented, config=config)]
    for line in nx.generate_adjlist(G=tree, delimiter=config[TREE_ADJ_LIST_SEPARATOR]):
        result.append(comment_aware_string(string=line, commented=commented, config=config))
    result.append(comment_aware_string(string=config[TREE_DATA_SUFFIX], commented=commented, config=config))
    result.append(comment_aware_string(string=config[TREE_SUFFIX], commented=commented, config=config))
    return "\n".join(result)


def mutation_config_to_string(mut_config, config, commented=True):
    result = [
        comment_aware_string(string=config[MUT_CONFIG_HEADER], config=config, commented=commented),
        comment_aware_string(string=pprint.pformat(object=mut_config), commented=commented, config=config),
        comment_aware_string(string=config[MUT_CONFIG_FOOTER], config=config, commented=commented),
    ]
    return "\n".join(result)


def edge_to_string(manager, edge, config, commented=False, comment_data=True, output_intermediate_genomes=False, output_mut_config=False):
    result = [comment_aware_string(string=config[EDGE_HEADER], commented=commented, config=config),
              comment_aware_string(string=config[TREE_ADJ_LIST_SEPARATOR].join(map(str, edge)), commented=commented, config=config),
              comment_aware_string(string=config[EDGE_DATA_HEADER], commented=commented, config=config)]
    comment_data = commented or comment_data
    u, v = edge
    mut_config = manager.tree.edges[u, v][MUT_CONFIG]
    if output_mut_config:
        result.append(mutation_config_to_string(mut_config=mut_config, config=config, commented=True))
    history = manager.tree.edges[u, v][HISTORY]
    intermediate_genomes = [] if len(history[GENOMES]) == 1 else history[GENOMES][:-1]
    mutations = history[MUTATIONS]
    for mutation, genome in itertools.zip_longest(mutations, intermediate_genomes):
        result.append(mutation_to_string(mutation=mutation, config=config, commented=comment_data))
        if genome is None or not output_intermediate_genomes:
            continue
        result.append(genome_to_string(genome=genome, config=config, commented=comment_data))
    result.append(comment_aware_string(string=config[EDGE_DATA_FOOTER], commented=commented, config=config))
    result.append(comment_aware_string(string=config[EDGE_FOOTER], commented=commented, config=config))
    return "\n".join(result)


def node_to_string(manager, node, config, commented=False):
    genome = manager.tree.nodes[node][GENOME]
    result = [
        comment_aware_string(string=config[GENOME_HEADER], config=config, commented=commented),
        comment_aware_string(string=str(node), config=config, commented=commented),
        comment_aware_string(string=config[GENOME_DATA_HEADER], config=config, commented=commented),
        genome_to_string(genome=genome, config=config, commented=commented),
        comment_aware_string(string=config[GENOME_DATA_FOOTER], config=config, commented=commented),
        comment_aware_string(string=config[GENOME_FOOTER], config=config, commented=commented),
    ]
    return "\n".join(result)


def write_genome_to_file(genome, file_name, config=None):
    if config is None:
        config = IO_CONFIG
    with open(file_name, "wt") as destination:
        print(genome_to_string(genome, config=config), file=destination)


def write_phylogenetic_mutation_history_to_file(manager, file_name, config=None):
    if config is None:
        config = IO_CONFIG
    processed_nodes = set()
    with open(file_name, "wt") as dest:
        print(tree_to_string(tree=manager.tree, config=config, commented=False), file=dest)
        for edge in nx.bfs_edges(G=manager.tree, source=manager.root):
            u, v = edge
            if u not in processed_nodes:
                print(node_to_string(manager=manager, node=u, config=config, commented=False), file=dest)
                processed_nodes.add(u)
            print(edge_to_string(manager=manager, edge=edge, config=config, commented=False), file=dest)
            if v not in processed_nodes:
                print(node_to_string(manager=manager, node=v, config=config, commented=False), file=dest)
                processed_nodes.add(v)


def read_phylogenetic_mutation_history_from_file(file_name, config=None):
    if config is None:
        config = IO_CONFIG
    genomes_by_ids = {}
    tree = None
    in_tree = False
    in_tree_data = False
    in_edge = False
    in_edge_data = False
    in_node = False
    in_node_data = False
    with open(file_name, "rt") as source:
        current_genome_name = None
        current_genome = Genome()
        adj_list = []
        for line in source:
            line = line.strip()
            if len(line) == 0 or line.startswith(config[COMMENT_INDICATOR]):
                continue
            if line == config[TREE_PREFIX]:
                in_tree = True
                continue
            if line == config[TREE_SUFFIX]:
                in_tree = False
                continue
            if line == config[TREE_DATA_PREFIX]:
                in_tree_data = True
                continue
            if line == config[TREE_DATA_SUFFIX]:
                tree = nx.parse_adjlist(lines=adj_list, create_using=nx.DiGraph(), nodetype=str)
                in_tree_data = False
                continue
            if line == config[EDGE_HEADER]:
                in_edge = True
                continue
            if line == config[EDGE_FOOTER]:
                in_edge = False
                continue
            if line == config[EDGE_DATA_HEADER]:
                in_edge_data = True
                continue
            if line == config[EDGE_DATA_FOOTER]:
                in_edge_data = False
                continue
            if line == config[GENOME_HEADER]:
                current_genome_name = None
                current_genome = Genome()
                in_node = True
                continue
            if line == config[GENOME_FOOTER]:
                if current_genome_name is None:
                    raise ParseError()
                genomes_by_ids[current_genome_name] = current_genome
                in_node = False
                continue
            if line == config[GENOME_DATA_HEADER]:
                in_node_data = True
                continue
            if line == config[GENOME_DATA_FOOTER]:
                in_node_data = False
                continue
            ###
            #
            # parsing data for an entry is only possible if we are in a respectful entry
            #
            ###
            if in_tree_data and not in_tree:
                raise ParseError()
            if in_edge_data and not in_edge:
                raise ParseError()
            if in_node_data and not in_node:
                raise ParseError()
            entries = sum([in_tree, in_edge, in_node])
            if entries == 0:
                # not in an entry, skipping line
                continue
            if entries > 1:
                # can only be in at most one entry at a time
                raise ParseError()
            ###
            #
            # we are in one of the three cases (tree | edge | node)
            #
            ###
            if in_tree:
                if in_tree_data:
                    adj_list.append(line)
                    continue
                else:
                    continue
            if in_edge:
                if in_edge_data:
                    continue
                else:
                    continue
            if in_node:
                if in_node_data:
                    chromosome = parse_chromosome_string(string=line, config=config)
                    current_genome.append(chromosome)
                    continue
                else:
                    # first line in "node" meta-data section is set as current_genome_name
                    current_genome_name = line if current_genome_name is None else current_genome_name
                    continue
        if tree is None:
            raise ParseError()
        for genome_name in genomes_by_ids:
            if genome_name not in tree:
                raise ParseError()
            tree.nodes[genome_name][GENOME] = genomes_by_ids[genome_name]
        result = Manager(phylo_tree=tree, assign_artificial_root=False)
        return result


def write_mutation_history_to_file(history, file_name, config=None):
    if config is None:
        config = IO_CONFIG
    with open(file_name, "wt") as destination:
        for genome, mutation in itertools.zip_longest(history[GENOMES], history[MUTATIONS]):
            comment_genome = mutation is not None  # last genome must be output as a single non-commented entry
            genome_string = genome_to_string(genome=genome, commented=comment_genome, config=config)
            mutation_string = "" if mutation is None else mutation_to_string(mutation=mutation, config=config, commented=True)
            print(genome_string, file=destination)
            print(mutation_string, file=destination)


def parse_chromosome_string(string, config):
    result = Chromosome()
    if string.startswith(config[CHROMOSOME_PREFIX]):
        string = string[len(config[CHROMOSOME_PREFIX]):]
    if string.endswith(config[CHROMOSOME_SUFFIX]) and len(config[CHROMOSOME_SUFFIX]) > 0:
        string = string[:-len(config[CHROMOSOME_SUFFIX])]
    string = string.strip()
    segments_entries = string.split(config[SEGMENTS_SEPARATOR])
    for string_entry in segments_entries:
        segment = Segment.from_string(string=string_entry)
        result.append(segment)
    return result


def read_genomes_from_file(file_name, config=None):
    if config is None:
        config = IO_CONFIG
    current_genome = []
    genomes = []
    in_genome = False
    with open(file_name, "rt") as source:
        for line in source:
            line = line.strip()
            if len(line) == 0 or line.startswith(config[COMMENT_INDICATOR]):
                continue
            if line == config[GENOME_HEADER]:
                in_genome = True
                current_genome = []
            elif line == config[GENOME_FOOTER]:
                genomes.append(current_genome)
                current_genome = []
                in_genome = False
            elif in_genome:
                chromosome = parse_chromosome_string(string=line, config=config)
                current_genome.append(chromosome)
    return genomes
