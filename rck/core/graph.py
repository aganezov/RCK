from copy import deepcopy

import networkx as nx

from rck.core.structures import propagate_haplotype_segment_to_positions, propagate_phasing_adjacency_to_positions
from rck.core.structures import SegmentCopyNumberProfile, AdjacencyCopyNumberProfile, StructureProfile, Segment, Adjacency, AdjacencyType, HAPLOTYPE, Haplotype, Phasing, PHASING

COPY_NUMBER = "copy_number"


def edge_tuple_based_on_flag(u, v, attr, data):
    if data:
        return u, v, attr
    else:
        return u, v


def node_tuple_based_on_flag(n, attr, data):
    if data:
        return n, attr
    return n


class IntervalAdjacencyGraph(object):
    def __init__(self, segments=None, adjacencies=None):
        segments = segments if segments is not None else []
        adjacencies = adjacencies if adjacencies is not None else []
        self.segments = segments
        self.adjacencies = adjacencies
        self.graph = nx.MultiGraph()

    @property
    def has_proper_topology(self):
        for node, n_data in self.nodes(data=True):
            s_edges = list(self.segment_edges(nbunch=node, data=True))
            if len(s_edges) != 1:    # every vertex (i.e., segment's extremity, can have at most 1 segment edge incident to it)
                return False
            if s_edges[0][0] == s_edges[0][1]:  # segments can not be represented in a form of self-loops
                return False
        return True

    @classmethod
    def get_segment_object_from_segment(cls, segment, copy=True):
        return segment.get_non_hap_copy()

    @classmethod
    def get_adjacency_object_from_adjacency(cls, adjacency, copy=True):
        return adjacency.get_non_phased_copy()

    def add_segment_edge(self, segment, sort=True, copy_segment=True):
        u, v = self.get_edge_vertices_pair_from_segment(segment=segment, sort=sort)
        if self.has_segment_edge_by_edge(edge=(u, v), sort=True):
            return
        obj = self.get_segment_object_from_segment(segment=segment, copy=copy_segment)
        self.graph.add_edge(u, v, object=obj)

    def add_adjacency_edge(self, adjacency, sort=True, copy_adjacency=True):
        u, v = self.get_edge_vertices_pair_from_adjacency(adjacency=adjacency, sort=sort)
        if u not in self.graph:
            raise ValueError()
        if v not in self.graph:
            raise ValueError()
        obj = self.get_adjacency_object_from_adjacency(adjacency=adjacency, copy=copy_adjacency)
        for adj_edge_w_data in self.adjacency_edges(nbunch=u):
            adj_u, adj_v, adj_data = adj_edge_w_data
            if {adj_u, adj_v} == {u, v} and adj_data["object"].adjacency_type == obj.adjacency_type:
                return
        # if self.has_adjacency_edge(edge=(u, v), sort=True):   # can not have parallel adjacency edges. Only possible parallel edges are pairs of segment/adjacency ones
        #     return
        self.graph.add_edge(u, v, object=obj)

    @staticmethod
    def get_edge_vertices_pair_from_adjacency(adjacency, sort=True):
        u, v = adjacency.position1.get_non_hap_copy(), adjacency.position2.get_non_hap_copy()
        if sort:
            u, v = tuple(sorted([u, v]))
        return u, v

    @staticmethod
    def get_edge_vertices_pair_from_segment(segment, sort=True):
        u, v = segment.start_position.get_non_hap_copy(), segment.end_position.get_non_hap_copy()
        if sort and segment.is_reversed:
            u, v = v, u
        return u, v

    def build_graph(self):
        for s in self.segments:
            self.add_segment_edge(segment=s)
        for a in self.adjacencies:
            self.add_adjacency_edge(adjacency=a)

    def nodes(self, data=True):
        for n, attr in self.graph.nodes(data=True):
            yield node_tuple_based_on_flag(n=n, attr=attr, data=data)

    def edges(self, data=True, nbunch=None):
        for value in self.graph.edges(nbunch=nbunch, data=data):
            yield value

    def segment_edges(self, data=True, nbunch=None, sort=True):
        for u, v, attr in self.edges(nbunch=nbunch, data=True):
            if isinstance(attr["object"], Segment):
                if sort:
                    u, v = tuple(sorted([u, v]))
                yield edge_tuple_based_on_flag(u, v, attr, data)

    def adjacency_edges(self, data=True, nbunch=None, sort=True):
        for u, v, attr in self.edges(nbunch=nbunch, data=True):
            if isinstance(attr["object"], Adjacency):
                if sort:
                    u, v = tuple(sorted([u, v]))
                yield edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)

    def ref_adjacency_edges(self, data=True, nbunch=None, sort=True):
        for u, v, attr in self.adjacency_edges(data=True, nbunch=nbunch, sort=sort):
            if attr["object"].adjacency_type == AdjacencyType.REFERENCE:
                yield edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)

    def nov_adjacency_edges(self, data=True, nbunch=None, sort=True):
        for u, v, attr in self.adjacency_edges(data=True, nbunch=nbunch, sort=sort):
            if attr["object"].adjacency_type == AdjacencyType.NOVEL:
                yield edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)

    def get_segment_edge(self, node, data=True, sort=True):
        segment_edges = list(self.segment_edges(data=True, nbunch=node, sort=sort))
        if len(segment_edges) != 1:
            raise ValueError("Vertex {vertex} does no either 0 or more than 1 segment edge incident to it".format(vertex=str(node)))
        u = segment_edges[0][0]
        v = segment_edges[0][1]
        attr = segment_edges[0][2]
        return edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)

    def has_segment_edge_by_edge(self, edge, sort=True):
        u, v = edge
        if sort:
            u, v = tuple(sorted([u, v]))
        if u not in self.graph or v not in self.graph:
            return False
        try:
            _, __, u_s_edge_data = self.get_segment_edge(node=u, data=True, sort=True)
            _, __, v_s_edge_data = self.get_segment_edge(node=v, data=True, sort=True)
            if u_s_edge_data["object"] != v_s_edge_data["object"]:
                return False
            return True
        except ValueError:
            return False

    def has_adjacency_edge_by_edge(self, edge, sort=True):
        u, v, = edge
        if sort:
            u, v = tuple(sorted([u, v]))
        if u not in self.graph or v not in self.graph:
            return False
        try:
            exists = False
            for u_, u__, u_a_edge_data in self.adjacency_edges(nbunch=u, data=True, sort=True):
                for v_, v__, v_a_edge_data in self.adjacency_edges(nbunch=v, data=True, sort=True):
                    if u_a_edge_data["object"] == v_a_edge_data["object"]:
                        exists = True
            return exists
        except ValueError:
            return False

    @property
    def complies_with_is(self):
        for node in self.nodes(data=False):
            n_as = list(self.nov_adjacency_edges(nbunch=node))
            if len(n_as) > 1:
                return False
        return True

    @property
    def compatible_with_hsis(self):
        for node in self.nodes(data=False):
            n_as = list(self.nov_adjacency_edges(nbunch=node))
            if len(n_as) > 2:
                return False
        return True

    @property
    def complies_with_extremity_exclusivity(self):
        return len(self.violations_of_extremity_exclusivity()) == 0

    def violations_of_extremity_exclusivity(self):
        violations = []
        for node in self.nodes(data=False):
            nas = list(self.nov_adjacency_edges(nbunch=node))
            if len(nas) > 1:
                violations.append(node)
        return violations

    def adjacency_edges_connected_components_subgraphs(self, ref=True, nov=True, copy=True):
        adjacency_edge_only_iag = self.__class__()
        if ref:
            for u, v, data in self.ref_adjacency_edges(data=True, sort=True):
                adjacency_edge_only_iag.graph.add_edge(u, v, **data)
        if nov:
            for u, v, data in self.nov_adjacency_edges(data=True, sort=True):
                adjacency_edge_only_iag.graph.add_edge(u, v, **data)
        for adj_iag_cc in nx.connected_component_subgraphs(G=adjacency_edge_only_iag.graph, copy=copy):
            iag = self.__class__()
            iag.graph = adj_iag_cc
            yield iag

    def ref_adjacency_edges_connected_components_subgraphs(self, copy=True):
        for iag in self.adjacency_edges_connected_components_subgraphs(ref=True, nov=False, copy=copy):
            yield iag

    def nov_adjacency_edges_connected_components_subgraphs(self, copy=True):
        for iag in self.adjacency_edges_connected_components_subgraphs(ref=False, nov=True, copy=copy):
            yield iag

    def get_set_self_segment_edges(self, sort=True):
        result = set()
        for u, v in self.segment_edges(data=False, sort=sort):
            result.add((u, v))
        return result

    def get_set_self_adjacency_edges(self, sort=True):
        result = set()
        for u, v in self.adjacency_edges(data=False, sort=sort):
            result.add((u, v))
        return result

    def intersection_on_segment_edges(self, segments):
        """iterates over `segments` once"""
        other_edges = set()
        for s in segments:
            u, v = self.get_edge_vertices_pair_from_segment(segment=s, sort=True)
            other_edges.add((u, v))
        return self.get_set_self_segment_edges(sort=True).intersection(other_edges)

    def intersection_on_adjacency_edges(self, adjacencies):
        """iterates over `adjacencies` once"""
        other_edges = set()
        for a in adjacencies:
            u, v = self.get_edge_vertices_pair_from_adjacency(adjacency=a, sort=True)
            other_edges.add((u, v))
        return self.get_set_self_adjacency_edges(sort=True).intersection(other_edges)

    def relative_complements_on_segment_edges(self, segments):
        """iterates over `segments` once"""
        self_edges = self.get_set_self_segment_edges(sort=True)
        other_edges = set()
        for s in segments:
            u, v = self.get_edge_vertices_pair_from_segment(segment=s, sort=True)
            other_edges.add((u, v))
        intersection = self_edges.intersection(other_edges)
        self_rel_compl = self_edges - intersection
        other_rel_compl = other_edges - intersection
        return self_rel_compl, other_rel_compl

    def relative_complements_on_adjacency_edges(self, adjacencies):
        """iterates over `adjacencies` once"""
        self_edges = self.get_set_self_adjacency_edges(sort=True)
        other_edges = set()
        for a in adjacencies:
            u, v = self.get_edge_vertices_pair_from_adjacency(adjacency=a, sort=True)
            other_edges.add((u, v))
        intersection = self_edges.intersection(other_edges)
        self_rel_compl = self_edges - intersection
        other_rel_compl = other_edges - intersection
        return self_rel_compl, other_rel_compl

    def topology_matches_for_genome(self, genome):
        self_rel_compl_segments, genome_rel_compl_segments = self.relative_complements_on_segment_edges(segments=genome.iter_segments())
        if len(self_rel_compl_segments) != 0 or len(genome_rel_compl_segments) != 0:
            return False
        self_rel_compl_adjs, genome_rel_compl_adjs = self.relative_complements_on_adjacency_edges(adjacencies=genome.iter_adjacencies())
        if len(self_rel_compl_adjs) != 0 or len(genome_rel_compl_adjs) != 0:
            return False
        return True

    def topology_allows_for_genome(self, genome):
        self_rel_compl_segments, genome_rel_compl_segments = self.relative_complements_on_segment_edges(segments=genome.iter_segments())
        if len(genome_rel_compl_segments) != 0:
            return False
        self_rel_compl_adjs, genome_rel_compl_adjs = self.relative_complements_on_adjacency_edges(adjacencies=genome.iter_adjacencies())
        if len(genome_rel_compl_adjs) != 0:
            return False
        return True

    def represents_given_genome(self, genome):
        if not self.is_copy_number_aware:
            raise ValueError()
        scn_profile = SegmentCopyNumberProfile.from_genome(genome=genome)
        if not self.matches_segment_copy_number_profile(scn_profile=scn_profile):
            return False
        acn_profile = AdjacencyCopyNumberProfile.from_genome(genome=genome)
        if not self.matches_adjacency_copy_number_profile(acn_profile=acn_profile):
            return False
        return True

    @property
    def is_copy_number_aware(self):
        for u, v, data in self.edges(data=True):
            if COPY_NUMBER not in data:
                return False
        return True

    @property
    def represents_a_genome(self):
        if not self.is_copy_number_aware:
            return False
        for node in self.graph:
            if not self.has_non_negative_imbalance(node=node):
                return False
        return True

    def node_imbalance(self, node):
        u, v, data = self.get_segment_edge(node=node, data=True, sort=True)
        if COPY_NUMBER not in data:
            raise ValueError("Segment edge {p1}-{p2} does not have a copy numbe associated with it".format(p1=str(u), p2=str(u)))
        segment_copy_number = data[COPY_NUMBER]
        imbalance = segment_copy_number
        for u, v, data in self.adjacency_edges(data=True, nbunch=node, sort=True):
            if COPY_NUMBER not in data:
                raise ValueError()
            adj_copy_number = data[COPY_NUMBER]
            adj_copy_number_multiplier = 2 if u == v else 1
            imbalance -= (adj_copy_number_multiplier * adj_copy_number)
        return imbalance

    def has_non_negative_imbalance(self, node):
        return self.node_imbalance(node=node) >= 0

    def has_positive_imbalance(self, node):
        return self.node_imbalance(node=node) > 0

    def has_zero_imbalance(self, node):
        return self.node_imbalance(node=node) == 0

    def is_telomere(self, node):
        return self.has_positive_imbalance(node=node)

    def is_non_telomere(self, node):
        return self.has_zero_imbalance(node=node)

    def matches_segment_copy_number_profile(self, scn_profile):
        sid_keys_set_from_profile = set(scn_profile.sid_keys())
        observed_sid_keys_set_from_graph = set()
        for u, v, data in self.segment_edges(data=True, sort=True):
            segment_obj = data["object"]
            cn = data[COPY_NUMBER]
            sid = segment_obj.stable_id_non_hap
            if scn_profile.get_combined_cn(sid=sid, default=0) != cn:
                return False
            observed_sid_keys_set_from_graph.add(sid)
        non_observed_sid_keys_set = sid_keys_set_from_profile - observed_sid_keys_set_from_graph
        for non_observed_sid in non_observed_sid_keys_set:
            if scn_profile.get_combined_cn(sid=non_observed_sid, default=0) != 0:
                return False
        return True

    def remove_segment_edge(self, u, v, silent=True):
        has_edge = self.has_segment_edge_by_edge(edge=(u, v), sort=True)
        if not has_edge:
            if not silent:
                raise ValueError("Attempted to remove segment edge {p1}-{p2}, but it is not present in the graph".format(p1=str(u), p2=str(v)))
            return
        edges_to_remove = []
        for u, v, key, data in self.graph.edges(nbunch=u, data=True, keys=True):
            obj = data["object"]
            if isinstance(obj, Segment):
                edges_to_remove.append((u, v, key))
        for u, v, key in edges_to_remove:
            self.graph.remove_edge(u=u, v=v, key=key)

    def remove_adjacency_edge(self, u, v, silent=True):
        has_edge = self.has_adjacency_edge_by_edge(edge=(u, v), sort=True)
        if not has_edge:
            if not silent:
                raise ValueError("Attempted to remove adjacency edge {p1}-{p2}, but it is not present in the graph".format(p1=str(u), p2=str(v)))
            return
        edges_to_remove = []
        for u_, u__, u_key, u_data in self.graph.edges(nbunch=u, data=True, keys=True):
            if not isinstance(u_data["object"], Adjacency):
                continue
            for v_, v__, v_key, v_data in self.graph.edges(nbunch=v, data=True, keys=True):
                if not isinstance(v_data["object"], Adjacency):
                    continue
                if u_data["object"] == v_data["object"] and u_data["object"].adjacency_type == v_data["object"].adjacency_type:
                    assert v_key == u_key
                    edges_to_remove.append((u, v, v_key))
        for u, v, key in edges_to_remove:
            self.graph.remove_edge(u=u, v=v, key=key)

    def remove_edges_with_zero_cn(self, check_cn_awareness=False, clear_nodes_after=True):
        self.remove_adjacency_edges_with_zero_cn(check_cn_awareness=check_cn_awareness, clear_nodes_after=clear_nodes_after)
        self.remove_segment_edges_with_zero_cn(check_cn_awareness=check_cn_awareness, clear_nodes_after=clear_nodes_after)

    def remove_segment_edges_with_zero_cn(self, check_cn_awareness=False, clear_nodes_after=True):
        if check_cn_awareness and not self.is_copy_number_aware:
            raise ValueError("Graph does not have copy number information associated with (at least some) its edges.")
        segment_edges_to_remove = []
        for u, v, data in self.segment_edges(data=True, sort=True):
            cn = data[COPY_NUMBER]
            if cn == 0:
                segment_edges_to_remove.append((u, v))
        for u, v in segment_edges_to_remove:
            self.remove_segment_edge(u=u, v=v)
            if clear_nodes_after:
                if self.graph.degree[u] == 0:
                    self.graph.remove_node(n=u)
                if self.graph.degree[v] == 0:
                    self.graph.remove_node(n=v)

    def remove_adjacency_edges_with_zero_cn(self, check_cn_awareness=False, clear_nodes_after=True):
        self.remove_ref_adjacency_edges_with_zero_cn(check_cn_awareness=check_cn_awareness, clear_nodes_after=clear_nodes_after)
        self.remove_nov_adjacency_edges_with_zero_cn(check_cn_awareness=check_cn_awareness, clear_nodes_after=clear_nodes_after)

    def remove_ref_adjacency_edges_with_zero_cn(self, check_cn_awareness=False, clear_nodes_after=True):
        if check_cn_awareness and not self.is_copy_number_aware:
            raise ValueError("Graph does not have copy number data associated with its edges")
        adjacency_edges_to_remove = []
        for u, v, data in self.ref_adjacency_edges(data=True, sort=True):
            cn = data[COPY_NUMBER]
            if cn == 0:
                adjacency_edges_to_remove.append((u, v))
                # self.remove_adjacency_edge(u=u, v=v)
        for u, v in adjacency_edges_to_remove:
            self.remove_adjacency_edge(u=u, v=v)
            if clear_nodes_after:
                if self.graph.degree[u] == 0:
                    self.graph.remove_node(n=u)
                if self.graph.degree[v] == 0:
                    self.graph.remove_node(n=v)

    def remove_nov_adjacency_edges_with_zero_cn(self, check_cn_awareness=False, clear_nodes_after=True):
        if check_cn_awareness and not self.is_copy_number_aware:
            raise ValueError()
        adjacency_edges_to_remove = []
        for u, v, data in self.nov_adjacency_edges(data=True, sort=True):
            cn = data[COPY_NUMBER]
            if cn == 0:
                adjacency_edges_to_remove.append((u, v))
                # self.remove_adjacency_edge(u=u, v=v)
        for u, v, in adjacency_edges_to_remove:
            self.remove_adjacency_edge(u=u, v=v)
            if clear_nodes_after:
                if self.graph.degree[u] == 0:
                    self.graph.remove_node(n=u)
                if self.graph.degree[v] == 0:
                    self.graph.remove_node(n=v)

    def matches_adjacency_copy_number_profile(self, acn_profile):
        aid_keys_set_from_profile = set(acn_profile.aid_keys())
        observed_aid_keys_from_graph = set()
        for u, v, data in self.adjacency_edges(data=True, sort=True):
            adjacency_obj = data["object"]
            cn = data[COPY_NUMBER]
            aid = adjacency_obj.stable_id_non_phased
            if acn_profile.get_combined_cn(aid=aid, default=0) != cn:
                return False
            observed_aid_keys_from_graph.add(aid)
        non_observed_aid_keys_set = aid_keys_set_from_profile - observed_aid_keys_from_graph
        for non_observed_aid in non_observed_aid_keys_set:
            if acn_profile.get_combined_cn(aid=non_observed_aid, default=0) != 0:
                return False
        return True

    def assign_copy_numbers_from_genome(self, genome, ensure_topology=True):
        if ensure_topology and not self.topology_allows_for_genome(genome=genome):
            raise ValueError()
        self.assign_copy_numbers_from_structure_profile(structure_profile=StructureProfile.from_genome(genome=genome))

    def assign_copy_numbers_from_structure_profile(self, structure_profile):
        self.assign_copy_numbers_from_scn_profile(scn_profile=structure_profile.scn_profile)
        self.assign_copy_numbers_from_acn_profile(acn_profile=structure_profile.acn_profile)

    def assign_copy_numbers_from_scn_profile(self, scn_profile):
        for u, v, data in self.segment_edges(data=True, sort=True):
            data[COPY_NUMBER] = scn_profile.get_combined_cn(sid=data["object"].stable_id_non_hap, default=0)

    def assign_copy_numbers_from_acn_profile(self, acn_profile):
        for u, v, data in self.adjacency_edges(data=True, sort=True):
            data[COPY_NUMBER] = acn_profile.get_combined_cn(aid=data["object"].stable_id_non_phased, default=0)

    def has_adjacency_edge(self, edge, sort=True):
        try:
            self.get_adjacency_edge(edge=edge, data=False, sort=sort)
            return True
        except ValueError:
            return False

    def get_adjacency_edge(self, edge, data=True, sort=True):
        u, v = edge
        if sort:
            u, v = tuple(sorted([u, v]))
        u_adjacency_edges = list(self.adjacency_edges(data=True, nbunch=u, sort=True))
        if len(u_adjacency_edges) == 0:
            raise ValueError()
        v_adjacency_edges = list(self.adjacency_edges(data=True, nbunch=v, sort=True))
        if len(v_adjacency_edges) == 0:
            raise ValueError()
        shared_edge = None
        for uu, uv, uattr in u_adjacency_edges:
            if (uu, uv) != (u, v):
                continue
            for vu, vv, vattr in v_adjacency_edges:
                if (uu, uv) == (vu, vv) and uattr["object"] == vattr["object"]:
                    shared_edge = u, v, uattr
                    break
            if shared_edge is not None:
                break
        if shared_edge is None:
            raise ValueError()
        u, v, attr = shared_edge
        return edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)

    def get_telomeres(self, check_cn_awareness=True, sort=True, copy=True):
        result = [t for t in self.iter_telomeres(check_cn_awareness=check_cn_awareness)]
        if copy:
            result = deepcopy(result)
        if sort:
            result = sorted(result)
        return result

    def iter_telomeres(self, check_cn_awareness=True):
        if check_cn_awareness and not self.is_copy_number_aware:
            raise ValueError()
        for node in self.nodes(data=False):
            if self.is_telomere(node=node):
                yield node


IAG = IntervalAdjacencyGraph


class HaplotypeSpecificIntervalAdjacencyGraph(IntervalAdjacencyGraph):
    def __init__(self, segments=None, adjacencies=None):
        super(HaplotypeSpecificIntervalAdjacencyGraph, self).__init__(segments=segments, adjacencies=adjacencies)

    @classmethod
    def get_segment_object_from_segment(cls, segment, copy=True):
        return deepcopy(segment) if copy else segment

    @classmethod
    def get_adjacency_object_from_adjacency(cls, adjacency, copy=True):
        return deepcopy(adjacency) if copy else adjacency

    def add_segment_edge(self, segment, sort=True, copy_segment=True):
        if not segment.is_haplotype_specific:
            raise ValueError()
        super(HaplotypeSpecificIntervalAdjacencyGraph, self).add_segment_edge(segment=segment, sort=sort)

    @classmethod
    def get_edge_vertices_pair_from_segment(cls, segment, sort=True):
        u, v = deepcopy(segment.start_position), deepcopy(segment.end_position)
        if sort and segment.is_reversed:
            u, v = v, u
        return u, v

    @classmethod
    def get_edge_vertices_pair_from_adjacency(cls, adjacency, sort=True):
        u, v = deepcopy(adjacency.position1), deepcopy(adjacency.position2)
        if sort:
            u, v = tuple(sorted([u, v]))
        return u, v

    def add_adjacency_edge(self, adjacency, sort=True, copy_adjacency=True):
        if not adjacency.is_phased:
            raise ValueError()
        super(HaplotypeSpecificIntervalAdjacencyGraph, self).add_adjacency_edge(adjacency=adjacency, sort=sort)

    def build_graph(self):
        for s in self.segments:
            self.add_segment_edge(segment=s)
        for a in self.adjacencies:
            self.add_adjacency_edge(adjacency=a)

    @property
    def complies_with_hiis(self):
        processed_nodes = set()
        for node, attr in self.nodes(data=True):
            n_as = list(self.nov_adjacency_edges(nbunch=node))
            if node in processed_nodes:
                continue
            position = node
            processed_nodes.add(position)
            p_haplotype = position.extra[HAPLOTYPE]
            hap_mates = [hap for hap in Haplotype if hap != p_haplotype]
            mate_positions = [deepcopy(position) for _ in hap_mates]
            for p, mh in zip(mate_positions, hap_mates):
                p.extra[HAPLOTYPE] = mh
                processed_nodes.add(p)
            for mp in mate_positions:
                if mp in self.graph:
                    n_as.extend(list(self.nov_adjacency_edges(nbunch=mp)))
            if len(n_as) > 1:
                return False
        return True

    @property
    def complies_with_hsis(self):
        return self.complies_with_is

    @property
    def compatible_with_hsis(self):
        return self.complies_with_hsis

    @property
    def complies_with_homologous_extremity_exclusivity(self):
        return len(self.violations_of_homologous_extremity_exclusivity()) == 0

    def violations_of_homologous_extremity_exclusivity(self):
        violations = []
        processed_nodes = set()
        for node, attr in self.nodes(data=True):
            if node in processed_nodes:
                continue
            n_as = list(self.nov_adjacency_edges(nbunch=node))
            if len(n_as) == 0:
                continue
            position = node
            processed_nodes.add(position)
            p_haplotype = position.extra[HAPLOTYPE]
            hap_mates = [hap for hap in (Haplotype.A, Haplotype.B) if hap != p_haplotype]
            mate_positions = [position.get_non_hap_copy() for _ in hap_mates]
            for p, mh in zip(mate_positions, hap_mates):
                p.extra[HAPLOTYPE] = mh
                processed_nodes.add(p)
            for mp in mate_positions:
                if mp in self.graph:
                    mate_nas = list(self.nov_adjacency_edges(nbunch=mp))
                    if len(mate_nas) > 0:
                        violations.append((node, mate_positions))
                        break
        return violations

    @property
    def complies_with_homologous_reciprocal_extremity_exclusivity(self):
        return len(self.violations_of_homologous_reciprocal_extremity_exclusivity()) == 0

    def violations_of_homologous_reciprocal_extremity_exclusivity(self):
        violators = []
        processed_nodes = set()
        for node, attr in self.nodes(data=True):
            if node in processed_nodes:
                continue
            ref_adjacencies = list(self.ref_adjacency_edges(nbunch=node))
            if len(ref_adjacencies) == 0:
                continue
            nas = list(self.nov_adjacency_edges(nbunch=node))
            if len(nas) == 0:
                continue
            for ref_adjacency in ref_adjacencies:
                u, v, data = ref_adjacency
                reciprocal_mate = u if v == node else v
                p_haplotype = reciprocal_mate.extra[HAPLOTYPE]
                hap_mates = [hap for hap in (Haplotype.A, Haplotype.B) if hap != p_haplotype]
                mate_positions = [reciprocal_mate.get_non_hap_copy() for _ in hap_mates]
                for p, mh in zip(mate_positions, hap_mates):
                    p.extra[HAPLOTYPE] = mh
                    processed_nodes.add(p)
                for mp in mate_positions:
                    if mp in self.graph:
                        if len(list(self.nov_adjacency_edges(nbunch=mp))) > 0:
                            violators.append((node, mp))
        return violators

    def matches_segment_copy_number_profile(self, scn_profile):
        sid_hap_pairs_set_from_profile = set(scn_profile.sid_hap_pairs())
        observed_sid_hap_pairs_set_from_graph = set()
        for u, v, data in self.segment_edges(data=True, sort=True):
            segment_object = data["object"]
            cn = data[COPY_NUMBER]
            hap = segment_object.haplotype
            sid = segment_object.stable_id_non_hap
            if scn_profile.get_cn(sid=sid, haplotype=hap, default=0) != cn:
                return False
            observed_sid_hap_pairs_set_from_graph.add((sid, hap))
        non_observed_sid_hap_pairs_set = sid_hap_pairs_set_from_profile - observed_sid_hap_pairs_set_from_graph
        for sid, hap in non_observed_sid_hap_pairs_set:
            if scn_profile.get_cn(sid=sid, hap=hap, default=0) != 0:
                return False
        return True

    def matches_adjacency_copy_number_profile(self, acn_profile):
        aid_phase_pairs_set_from_profile = set(acn_profile.aid_phase_pairs())
        observed_aid_phase_pairs_set_from_graph = set()
        for u, v, data in self.adjacency_edges(data=True, sort=True):
            adj_object = data["object"]
            cn = data[COPY_NUMBER]
            aid = adj_object.stable_id_non_phased
            phasing = adj_object.stable_phasing
            if acn_profile.get_cn(aid=aid, phasing=phasing, default=0) != cn:
                return False
            observed_aid_phase_pairs_set_from_graph.add((aid, phasing))
        non_observed_aid_phase_pairs_set = aid_phase_pairs_set_from_profile - observed_aid_phase_pairs_set_from_graph
        for aid, phasing in non_observed_aid_phase_pairs_set:
            if acn_profile.get_cn(aid=aid, phasing=phasing, default=0) != 0:
                return False
        return True

    def assign_copy_numbers_from_scn_profile(self, scn_profile):
        for u, v, data in self.segment_edges(data=True, sort=True):
            sid = data["object"].stable_id_non_hap
            hap = data["object"].haplotype
            data[COPY_NUMBER] = scn_profile.get_cn(sid=sid, haplotype=hap, default=0)

    def assign_copy_numbers_from_acn_profile(self, acn_profile):
        for u, v, data in self.adjacency_edges(data=True, sort=True):
            adj = data["object"]
            aid = adj.stable_id_non_phased
            phasing = adj.stable_phasing
            data[COPY_NUMBER] = acn_profile.get_cn(aid=aid, phasing=phasing, default=0)


def construct_iag(ref_genome, mut_genomes):
    ref_adjacencies = {}
    for adj in ref_genome.iter_adjacencies(adjacency_type=AdjacencyType.REFERENCE):
        aid = adj.stable_id_non_phased
        if aid not in ref_adjacencies:
            ref_adjacencies[aid] = adj
    nov_adjacencies = {}
    for mut_genome in mut_genomes:
        for adj in mut_genome.iter_adjacencies(adjacency_type=AdjacencyType.NOVEL):
            aid = adj.stable_id_non_phased
            if aid not in ref_adjacencies and aid not in nov_adjacencies:
                nov_adjacencies[aid] = adj
    segments = {}
    for s in ref_genome.iter_segments():
        sid = s.stable_id_non_hap
        if sid not in segments:
            segments[sid] = s
    result = IntervalAdjacencyGraph()
    for segment in segments.values():
        result.add_segment_edge(segment=segment, sort=True, copy_segment=True)
    for adj in ref_adjacencies.values():
        result.add_adjacency_edge(adjacency=adj, sort=True, copy_adjacency=True)
    for adj in nov_adjacencies.values():
        result.add_adjacency_edge(adjacency=adj, sort=True, copy_adjacency=True)
    return result


def construct_hiag(ref_genome, mut_genomes):
    ref_adjacencies = {}
    for adj in ref_genome.iter_adjacencies(adjacency_type=AdjacencyType.REFERENCE):
        aid = adj.stable_id_phased
        if aid not in ref_adjacencies:
            ref_adjacencies[aid] = adj
    nov_adjacencies = {}
    for mut_genome in mut_genomes:
        for adj in mut_genome.iter_adjacencies(adjacency_type=AdjacencyType.NOVEL):
            aid = adj.stable_id_phased
            if aid not in ref_adjacencies and aid not in nov_adjacencies:
                nov_adjacencies[aid] = adj
    segments = {}
    for s in ref_genome.iter_segments():
        sid = s.stable_id_hap
        if sid not in segments:
            segments[sid] = s
    result = HaplotypeSpecificIntervalAdjacencyGraph()
    for segment in segments.values():
        result.add_segment_edge(segment=segment, sort=True, copy_segment=True)
    for adj in ref_adjacencies.values():
        result.add_adjacency_edge(adjacency=adj, sort=True, copy_adjacency=True)
    for adj in nov_adjacencies.values():
        result.add_adjacency_edge(adjacency=adj, sort=True, copy_adjacency=True)
    return result


def construct_hiag_inflate_from_haploid_data(hapl_segments, hapl_adjacencies):
    dip_ref_adjacencies = []
    for adj in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, hapl_adjacencies):
        for ph in [Phasing.AA, Phasing.BB]:
            dip_adjacency = adj.get_non_phased_copy()
            dip_adjacency.extra[PHASING] = ph
            dip_ref_adjacencies.append(dip_adjacency)
    dip_nov_adjacencies = []
    for adj in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, hapl_adjacencies):
        for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
            dip_adjacency = adj.get_non_phased_copy()
            dip_adjacency.extra[PHASING] = ph
            dip_nov_adjacencies.append(dip_adjacency)
    dip_segments = []
    for segment in hapl_segments:
        for h in [Haplotype.A, Haplotype.B]:
            dip_segment = segment.get_non_hap_copy()
            dip_segment.extra[HAPLOTYPE] = h
            dip_segments.append(dip_segment)
    for s in dip_segments:
        propagate_haplotype_segment_to_positions(segment=s, inplace=True)
    all_dip_adjacencies = dip_ref_adjacencies + dip_nov_adjacencies
    for a in all_dip_adjacencies:
        propagate_phasing_adjacency_to_positions(adjacency=a, inplace=True)
    hiag = HaplotypeSpecificIntervalAdjacencyGraph(segments=dip_segments, adjacencies=dip_ref_adjacencies + dip_nov_adjacencies)
    hiag.build_graph()
    return hiag


construct_hsiag = construct_hiag

HSIAG = HaplotypeSpecificIntervalAdjacencyGraph
HIAG = HSIAG
