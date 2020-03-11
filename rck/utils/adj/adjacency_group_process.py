from collections import defaultdict
from copy import deepcopy

import networkx as nx

from rck.core.structures import AdjacencyGroup, AdjacencyGroupType
from rck.core.io import AG_LABELING, EXTERNAL_NA_ID


def refine_labeling_groups_old(adj_groups, gid_suffix="", retain_source_gids=False, iag=None):
    graph = nx.Graph()
    entries_to_adj_groups = defaultdict(list)
    groups_by_ids = {}
    for group in adj_groups:
        groups_by_ids[group.gid] = group
        entries = [(aid, index) for aid, index in zip(group.adjacencies_ids, group.extra.get(AG_LABELING, []))]
        if len(entries) < 2:
            continue
        for entry in entries:
            entries_to_adj_groups[entry].append(group)
        for l, r in zip(entries[:-1], entries[1:]):
            graph.add_edge(l, r)
    result = []
    cnt = 0
    for cc in nx.connected_component_subgraphs(graph):
        entries = list(cc.nodes())
        gids = set()
        for entry in entries:
            groups = entries_to_adj_groups[entry]
            for group in groups:
                gids.add(group.gid)
        groups = [groups_by_ids[gid] for gid in gids]
        aids = [entry[0] for entry in entries]
        indexes = [entry[1] for entry in entries]
        alignments = list(set(group.extra.get("alignment", "") for group in groups if len(group.extra.get("alignment", "")) > 0))
        extra = {
            AG_LABELING: indexes,
        }
        if len(alignments) > 0:
            extra["alignment"] = alignments
        if retain_source_gids:
            extra["source_gids"] = sorted(gids)
        gid = str(cnt)
        if len(gid_suffix) > 0:
            gid += "_" + gid_suffix
        ag = AdjacencyGroup(gid=gid, aids=aids, group_type=AdjacencyGroupType.LABELING, extra=extra)
        result.append(ag)
        cnt += 1
    return result


def refine_labeling_groups_without_iag(adj_groups, gid_suffix="", retain_source_gids=False):
    graph = nx.Graph()
    entries_to_adj_groups = defaultdict(list)
    groups_by_ids = {}
    for group in adj_groups:
        groups_by_ids[group.gid] = group
        entries = [(aid, index) for aid, index in zip(group.adjacencies_ids, group.extra.get(AG_LABELING, []))]
        if len(entries) < 2:
            continue
        for entry in entries:
            entries_to_adj_groups[entry].append(group)
        for l, r in zip(entries[:-1], entries[1:]):
            graph.add_edge(l, r)
    result = []
    cnt = 0
    for cc in nx.connected_component_subgraphs(graph):
        entries = list(cc.nodes())
        gids = set()
        for entry in entries:
            groups = entries_to_adj_groups[entry]
            for group in groups:
                gids.add(group.gid)
        groups = [groups_by_ids[gid] for gid in gids]
        aids = [entry[0] for entry in entries]
        indexes = [entry[1] for entry in entries]
        alignments = list(set(group.extra.get("alignment", "") for group in groups if len(group.extra.get("alignment", "")) > 0))
        extra = {
            AG_LABELING: indexes,
        }
        if len(alignments) > 0:
            extra["alignment"] = alignments
        if retain_source_gids:
            extra["source_gids"] = sorted(gids)
        gid = str(cnt)
        if len(gid_suffix) > 0:
            gid += "_" + gid_suffix
        ag = AdjacencyGroup(gid=gid, aids=aids, group_type=AdjacencyGroupType.LABELING, extra=extra)
        result.append(ag)
        cnt += 1
    return result


def refined_labeling_groups(adj_groups, iag=None, adjacencies=None, gid_suffix="", retain_source_gids=False):
    graph = nx.Graph()
    if adjacencies is not None:
        adjacencies_by_external_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in adjacencies}
    else:
        adjacencies_by_external_ids = {}
    adjacencies_by_positions = defaultdict(list)
    groups_by_ids = {}
    entries_to_adj_groups = defaultdict(list)
    for group in adj_groups:
        groups_by_ids[group.gid] = group
        internal_entries = [(aid, index) for aid, index in zip(group.adjacencies_ids, group.extra.get(AG_LABELING, []))]
        if len(internal_entries) < 2:
            continue
        entries = []
        if iag is not None and adjacencies is not None:
            for aid, index in internal_entries:
                adjacency = adjacencies_by_external_ids[aid]
                position = adjacency.position1 if index == 0 else adjacency.position2
                adjacencies_by_positions[position].append(adjacency)
                entries.append(position)
        else:
            entries = internal_entries
        for entry in entries:
            entries_to_adj_groups[entry].append(group)
        for l, r in zip(entries[:-1], entries[1:]):
            graph.add_edge(l, r)
            if iag is not None and adjacencies is not None:
                l_ref_edges = list(iag.ref_adjacency_edges(nbunch=l, data=False))
                r_ref_edges = list(iag.ref_adjacency_edges(nbunch=r, data=False))
                ref_edges = l_ref_edges + r_ref_edges
                for (u, v) in ref_edges:
                    graph.add_edge(u, v)
    result = []
    cnt = 0
    for cc in nx.connected_component_subgraphs(graph):
        internal_entries = list(cc.nodes())
        gids = set()
        for entry in internal_entries:
            groups = entries_to_adj_groups[entry]
            for group in groups:
                gids.add(group.gid)
        entries = []
        if iag is not None and adjacencies is not None:
            for entry in internal_entries:
                if entry in adjacencies_by_positions:
                    for adjacency in adjacencies_by_positions[entry]:
                        aid = adjacency.extra.get(EXTERNAL_NA_ID, adjacency.stable_id_non_phased)
                        index = 0 if adjacency.position1 == entry else 1
                        entries.append((aid, index))
        else:
            entries = internal_entries
        groups = [groups_by_ids[gid] for gid in gids]
        aids = [entry[0] for entry in entries]
        indexes = [entry[1] for entry in entries]
        alignments = list(set(group.extra.get("alignment", "") for group in groups if len(group.extra.get("alignment", "")) > 0))
        extra = {
            AG_LABELING: indexes,
        }
        if len(alignments) > 0:
            extra["alignment"] = alignments
        if retain_source_gids:
            extra["source_gids"] = sorted(gids)
        gid = str(cnt)
        if len(gid_suffix) > 0:
            gid += "_" + gid_suffix
        ag = AdjacencyGroup(gid=gid, aids=aids, group_type=AdjacencyGroupType.LABELING, extra=extra)
        result.append(ag)
        cnt += 1
    return result


def projected_groups(groups, adjacencies, adjacencies_by_external_ids=None, gid_suffix=""):
    if adjacencies_by_external_ids is None:
        adjacencies_by_external_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in adjacencies}
    result = []
    for group in groups:
        projected = [aid in adjacencies_by_external_ids for aid in group.adjacencies_ids]
        aids = [aid for aid, allowed in zip(group.adjacencies_ids, projected) if allowed]
        if len(aids) == 0:
            continue
        if group.adjacencies is not None:
            adjacencies = [adj for adj, allowed in zip(group.adjacencies, projected) if allowed]
        else:
            adjacencies = None
        ag = deepcopy(group)
        ag.adjacencies_ids = aids
        ag.adjacencies = adjacencies
        if group.group_type == AdjacencyGroupType.LABELING:
            if len(aids) < 2:
                continue
            labeling = [index for index, allowed in zip(group.extra[AG_LABELING], projected) if allowed]
            ag.extra[AG_LABELING] = labeling
        if not all(projected) and len(gid_suffix) > 0 and group.group_type != AdjacencyGroupType.GENERAL:
            ag.gid += "_" + gid_suffix
        result.append(ag)
    return result
