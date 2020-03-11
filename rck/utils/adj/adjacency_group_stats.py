from collections import defaultdict

from rck.core.structures import AdjacencyGroup


def groups_size_tally(adjacency_groups):
    result = defaultdict(int)
    for ag in adjacency_groups:
        ag: AdjacencyGroup = ag
        result[len(ag.adjacencies_ids)] += 1
    return result
