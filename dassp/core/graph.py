import networkx as nx

from dassp.core.structures import Segment, Adjacency, AdjacencyType


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
        self.internal_check_consistency()
        self.build_graph()

    def internal_check_consistency(self):
        self.check_consistency(segments=self.segments, adjacencies=self.adjacencies)

    @classmethod
    def check_consistency(cls, segments, adjacencies):
        for a in adjacencies:
            if a.adjacency_type == AdjacencyType.NOVEL:
                continue
            p1, p2 = a.position1, a.position2
            if p1.chromosome != p2.chromosome:
                raise ValueError("Reference adjacency {a} links positions {p1} and {p2} from different chromosomes."
                                 "".format(a=str(a), p1=str(p1), p2=str(p2)))
        segments_extremities = set()
        for s in segments:
            segments_extremities.add(s.start_position)
            segments_extremities.add(s.end_position)
        for a in adjacencies:
            if a.position1 not in segments_extremities or a.position2 not in segments_extremities:
                raise ValueError("Adjacency {a} links positions {p1} and {p2} that do not correspond to segments' extremities positions."
                                 "".format(a=str(a), p1=str(a.position1), p2=str(a.position2)))

    def build_graph(self):
        for s in self.segments:
            u = s.start_position.idx
            v = s.end_position.idx
            if u not in self.graph:
                self.graph.add_node(n=u, object=s.start_position)
            if v not in self.graph:
                self.graph.add_node(n=v, object=s.end_position)
            self.graph.add_edge(u=u, v=v, object=s)
        for a in self.adjacencies:
            self.graph.add_edge(u=a.position1.idx, v=a.position2.idx, object=a)

    def nodes(self, data=True):
        for n, attr in self.graph.nodes(data=True):
            yield node_tuple_based_on_flag(n=n, attr=attr, data=data)

    def edges(self, data=True, nbunch=None):
        return self.graph.edges(nbunch=nbunch, data=data)

    def segment_edges(self, data=True, nbunch=None):
        for u, v, attr in self.edges(nbunch=nbunch, data=True):
            if isinstance(attr["object"], Segment):
                yield edge_tuple_based_on_flag(u, v, attr, data)

    def adjacency_edges(self, data=True, nbunch=None):
        for u, v, attr in self.edges(nbunch=nbunch, data=True):
            if isinstance(attr["object"], Adjacency):
                yield edge_tuple_based_on_flag(u, v, attr, data)

    def ref_adjacency_edges(self, data=True, nbunch=None):
        for u, v, attr in self.adjacency_edges(data=True, nbunch=nbunch):
            if attr["object"].adjacency_type == AdjacencyType.REFERENCE:
                yield edge_tuple_based_on_flag(u, v, attr, data)

    def nov_adjacency_edges(self, data=True, nbunch=None):
        for u, v, attr in self.adjacency_edges(data=True, nbunch=nbunch):
            if attr["object"].adjacency_type == AdjacencyType.NOVEL:
                yield edge_tuple_based_on_flag(u, v, attr, data)

    def get_segment_edge(self, node, data=True):
        segment_edges = list(self.segment_edges(data=True, nbunch=node))
        assert len(segment_edges) == 1
        u = segment_edges[0][0]
        v = segment_edges[0][1]
        attr = segment_edges[0][2]
        return edge_tuple_based_on_flag(u=u, v=v, attr=attr, data=data)


IAG = IntervalAdjacencyGraph
