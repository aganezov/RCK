from copy import deepcopy

from dassp.core.graph import IntervalAdjacencyGraph
from dassp.core.structures import strip_phasing_from_adjacencies, get_segments_list_from_genome, get_adjacencies_from_genome, AdjacencyType, assign_adjacency_status, \
    assign_ref_adjacency_status, get_unique_adjacencies, strip_haplotype_from_segments, get_telomeres_from_genome, strip_haplotype_from_positions
from dassp.core.graph import construct_iag
from dassp.core.structures import StructureProfile
import networkx as nx


def unique_ref_segments_presence(iag):
    for node in iag.nodes(data=False):
        adjs = list(iag.ref_adjacency_edges(data=True, nbunhc=node, sort=True))
        if len(adjs) > 1:
            return False
    return True


def linear_chromosomes_only(iag):
    for iag_cc in nx.connected_component_subgraphs(G=iag.graph, copy=False):
        telomere_cnt = 0
        for node in iag_cc:
            adjs = list(iag.adjacency_edges(data=True, nbunch=node, sort=True))
            if len(adjs) == 0:
                telomere_cnt += 1
        if telomere_cnt != 2:
            return False
    return True


class Window(object):
    def __init__(self, left_segment, right_segment, adjacency=None, sort=True, copy=False):
        if copy:
            adjacency = deepcopy(adjacency)
            left_segment = deepcopy(left_segment)
            right_segment = deepcopy(right_segment)
        self.adjacency = adjacency
        self.left_segment = left_segment
        self.right_segment = right_segment
        if sort and None not in [self.left_segment, self.right_segment] and self.left_segment > self.right_segment:
            self.left_segment, self.right_segment = self.right_segment, self.left_segment


class WindowsGroup(object):
    def __init__(self, windows, graph, sort=True, copy=False):
        if copy:
            windows = deepcopy(windows)
            graph = deepcopy(graph)
        self.windows = windows
        self.graph = graph
        if sort:
            self.windows = sorted(self.windows, key=lambda w: (w.left_segment, w.right_segment))


def windows_group_from_adjacency_iag(iag, check_connectedness=False, check_hiis_topology=False):
    if check_connectedness:
        if not nx.is_connected(G=iag.graph):
            raise ValueError()
    if check_hiis_topology and not iag.complies_with_is:
        raise ValueError()


def solve_single_hiis_decision(ref_genome,
                               novel_adjacencies,
                               mut_scn_profile,
                               telomere_vertices=None,
                               check_ref_topology=True,
                               check_nov_topology=True,
                               allow_hsis_compatible_topology=True):
    if check_ref_topology:
        iag = construct_iag(ref_genome=ref_genome, mut_genomes=[])
        if not unique_ref_segments_presence(iag=iag):
            raise ValueError()
        if not linear_chromosomes_only(iag=iag):
            raise ValueError()
    segments = get_segments_list_from_genome(genome=ref_genome, copy=True, make_all_non_reversed=True)
    nh_segments = strip_haplotype_from_segments(segments=segments, inplace=False, strip_positions_haplotypes=True)
    ref_adjacencies = get_adjacencies_from_genome(genome=ref_genome, copy=True, inherit_haplotypes=True, default_adjacency_type=AdjacencyType.REFERENCE)
    nh_ref_adjacencies = strip_phasing_from_adjacencies(adjacencies=ref_adjacencies, inplace=False, strip_positions_haplotypes=True, sort=True)
    assign_ref_adjacency_status(adjacencies=nh_ref_adjacencies, inplace=True)
    nh_ref_adjacencies_set = set(nh_ref_adjacencies)
    nh_nov_adjacencies = strip_phasing_from_adjacencies(adjacencies=novel_adjacencies, inplace=False, strip_positions_haplotypes=True, sort=True)
    assign_adjacency_status(adjacencies=nh_nov_adjacencies, ref_adjacencies=nh_ref_adjacencies_set, inplace=True)
    nh_nov_adjacencies = [a for a in nh_nov_adjacencies if a.adjacency_type == AdjacencyType.NOVEL]
    unique_nh_ref_adjacencies = get_unique_adjacencies(adjacencies=nh_ref_adjacencies, copy=False)
    unique_nh_nov_adjacencies = get_unique_adjacencies(adjacencies=nh_nov_adjacencies, copy=False)
    iag = IntervalAdjacencyGraph(segments=nh_segments, adjacencies=unique_nh_ref_adjacencies + unique_nh_nov_adjacencies)
    iag.build_graph()
    if check_nov_topology:
        if not allow_hsis_compatible_topology and not iag.complies_with_is:
            raise ValueError()
        elif allow_hsis_compatible_topology and not iag.compatible_with_hsis:
            raise ValueError()
    ref_telomeres = get_telomeres_from_genome(genome=ref_genome, copy=True, inherit_haplotypes=True)
    strip_haplotype_from_positions(positions=ref_telomeres, inplace=True)
    if telomere_vertices is not None:
        telomere_vertices = []
    telomere_vertices += ref_telomeres
    for adj_subgraph_iag in iag.adjacency_edges_connected_components_subgraphs(ref=True, nov=True, copy=False):
        pass
    result = StructureProfile()


def solve_single_hsis_decision(ref_genome,
                               novel_adjacencies,
                               mut_scn_profile,
                               telomere_vertices=None,
                               check_ref_topology=True,
                               check_nov_adjacency_topology=True):
    pass


def solve_single_decision():
    pass


