from typing import Iterable, List

import networkx as nx

from rck.core.structures import Adjacency, Position, Strand, AdjacencyType


class ComplexRearrSignature(object):
    def __init__(self, adjacencies: Iterable[Adjacency], ref_locations: Iterable[Adjacency] = None):
        self.adjacencies: List[Adjacency] = list(adjacencies)
        self.ref_adjacencies: List[Adjacency] = list(ref_locations) if ref_locations is not None else self.infer_ref_adjacencies(self.adjacencies)
        self.k = len(self.ref_adjacencies)

    @classmethod
    def infer_ref_adjacencies(cls, adjacencies: Iterable[Adjacency]) -> List[Adjacency]:
        result = set()
        for adjacency in adjacencies:
            for p in [adjacency.position1, adjacency.position2]:
                pr = Position.get_reciprocal(position=p)
                result.add(Adjacency(position1=p, position2=pr, adjacency_type=AdjacencyType.REFERENCE))
        return list(result)


def get_complex_rearrangements_signatures(adjacencies: Iterable[Adjacency]) -> Iterable[ComplexRearrSignature]:
    cr_graph = nx.MultiGraph()
    for adjacency in adjacencies:
        p1 = adjacency.position1.get_non_hap_copy()
        p2 = adjacency.position2.get_non_hap_copy()
        if p1.strand == Strand.REVERSE:
            p1 = Position.get_reciprocal(position=p1)
        if p2.strand == Strand.REVERSE:
            p2 = Position.get_reciprocal(position=p2)
        cr_graph.add_edge(p1, p2, adjacency=adjacency)
    result: List[ComplexRearrSignature] = []
    for cc in nx.connected_component_subgraphs(cr_graph):
        result.append(ComplexRearrSignature(adjacencies=[edge[2]["adjacency"] for edge in cc.edges(data=True)]))
    return result
