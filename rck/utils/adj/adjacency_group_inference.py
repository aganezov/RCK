import sys
from collections import defaultdict

import pysam

from rck.core.io import EXTERNAL_NA_ID, AG_LABELING
from rck.core.structures import AdjacencyGroup, AdjacencyGroupType, Strand


def infer_sniffles_molecule_groups(adjacencies, extra_rnames_field="rnames", gid_suffix=""):
    reads_to_adjacencies_ids = defaultdict(set)
    for adj in adjacencies:
        read_names = adj.extra.get(extra_rnames_field, "").split(",")
        if len(read_names) == 1 and len(read_names[0]) == 0:
            continue
        for read_name in read_names:
            reads_to_adjacencies_ids[read_name].add(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased))
    result = []
    for cnt, (read_name, aids) in enumerate(reads_to_adjacencies_ids.items()):
        if len(aids) < 2:
            continue
        extra = {"source": read_name}
        gid = str(cnt)
        if len(gid_suffix) > 0:
            gid += "_" + gid_suffix
        ag = AdjacencyGroup(gid=gid, aids=list(aids), group_type=AdjacencyGroupType.MOLECULE, extra=extra)
        result.append(ag)
    return result


def infer_short_nas_labeling_groups(adjacencies, gid_suffix="", max_size=1000, allow_intermediate_same=False, allow_intermediate_tra=False):
    positions_to_adjacencies = defaultdict(list)
    positions_by_chrs = defaultdict(list)
    result = []
    for adj in adjacencies:
        p1, p2 = adj.position1, adj.position2
        positions_to_adjacencies[p1].append(adj)
        positions_to_adjacencies[p2].append(adj)
        positions_by_chrs[p1.chromosome].append(p1)
        positions_by_chrs[p2.chromosome].append(p2)
    positions_by_chr_to_index = {}
    for chr_name in list(positions_by_chrs.keys()):
        positions_by_chrs[chr_name] = sorted(positions_by_chrs[chr_name], key=lambda p: (p.coordinate, p.strand))
        positions_by_chr_to_index[chr_name] = {p: cnt for cnt, p in enumerate(positions_by_chrs[chr_name])}
    processed_adj_ids = set()
    cnt = 0
    for adj in adjacencies:
        aid = adj.stable_id_non_phased
        if aid in processed_adj_ids:
            continue
        p1, p2 = adj.position1, adj.position2
        p1_chr, p2_chr = p1.chromosome, p2.chromosome
        if p1_chr != p2_chr:
            continue
        adj_size = adj.distance_non_hap
        if adj_size > max_size:
            continue
        positions = positions_by_chrs[p1_chr]
        p1_index, p2_index = positions_by_chr_to_index[p1_chr][p1], positions_by_chr_to_index[p1_chr][p2]
        aid = adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased)
        aids = [aid, aid]
        gid = str(cnt)
        if len(gid_suffix) > 0:
            gid += "_" + gid_suffix
        extra = {AG_LABELING: [0, 1]}
        ag = AdjacencyGroup(gid=gid, aids=aids, group_type=AdjacencyGroupType.LABELING, extra=extra)
        if abs(p1_index - p2_index) == 1:
            result.append(ag)
            cnt += 1
        else:
            if not (allow_intermediate_same or allow_intermediate_tra):
                continue
            intermediate_indexes = list(range(p1_index + 1, p2_index))
            has_same = False
            has_tra = False
            allow = False
            for index in intermediate_indexes:
                position = positions[index]
                adjs = positions_to_adjacencies[position]
                for adj in adjs:
                    has_same |= adj.position1.chromosome == adj.position2.chromosome
                    has_tra |= adj.position1.chromosome != adj.position2.chromosome
            allow |= has_same and allow_intermediate_same
            allow |= has_tra and allow_intermediate_tra
            if allow:
                result.append(ag)
                cnt += 1
        processed_adj_ids.add(adj.stable_id_non_phased)
    return result


def get_mode_str(format="bam", input=False):
    result = "r" if input else "w"
    if format == "bam":
        result += "b"
    elif format == "cram":
        result += "c"
    return result


def get_labeling_groups(read_alignments, read_adjacencies, strategy="skip", delta=500):
    result = []
    read_alignments = sorted(read_alignments, key=lambda e: (e.query_alignment_start, e.query_alignment_end))
    processed_positions = set()
    positions_by_chrs = defaultdict(list)
    positions_to_alignments = defaultdict(list)
    alignments_to_positions = defaultdict(list)
    positions_to_adjacencies = defaultdict(list)
    for adj in read_adjacencies:
        p1 = adj.position1
        p2 = adj.position2
        positions_by_chrs[p1.chromosome].append(p1)
        positions_by_chrs[p2.chromosome].append(p2)
        for alignment in read_alignments:
            for p in [p1, p2]:
                if p.chromosome == alignment.reference_name and \
                        (alignment.reference_start - delta <= p.coordinate <= alignment.reference_end + delta):
                    positions_to_alignments[p].append(alignment)
                    alignments_to_positions[alignment].append(p)
        positions_to_adjacencies[p1].append(adj)
        positions_to_adjacencies[p2].append(adj)
    for alignment in alignments_to_positions.keys():
        alignments_to_positions[alignment] = sorted(alignments_to_positions[alignment], key=lambda p: (p.coordinate, p.strand))
    for chr_name in list(positions_by_chrs.keys()):
        positions_by_chrs[chr_name] = sorted(positions_by_chrs[chr_name], key=lambda p: (p.coordinate, p.strand))
    cnt = 0
    for adj in read_adjacencies:
        aid = adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased)
        for p in [adj.position1, adj.position2]:
            if p in processed_positions:
                continue
            processed_positions.add(p)
            alignments = positions_to_alignments[p]
            for alignment in alignments:
                positions_on_alignment = alignments_to_positions[alignment]
                p_index = positions_on_alignment.index(p)
                direction_neighbours = positions_on_alignment[:p_index] if p.strand == Strand.FORWARD else positions_on_alignment[p_index + 1:]
                if len(direction_neighbours) == 0:
                    continue
                else:
                    neighbour = direction_neighbours[-1] if p.strand == Strand.FORWARD else direction_neighbours[0]
                    processed_positions.add(neighbour)
                    neighbour_ids = [adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased) for adj in positions_to_adjacencies[neighbour]]
                    adj_index = 0 if p == adj.position1 else 1
                    neighbour_indexes = [0 if neighbour == adj.position1 else 1 for adj in positions_to_adjacencies[neighbour]]
                    extra = {
                        "alignment": alignment.query_name,
                        AG_LABELING: [adj_index] + neighbour_indexes
                    }
                    ag = AdjacencyGroup(gid=cnt, aids=[aid] + neighbour_ids, group_type=AdjacencyGroupType.LABELING, extra=extra)
                    result.append(ag)
                    cnt += 1
                    break
    return result


def infer_alignment_labeling_groups(adjacencies, alignment_file_name, alignment_format="bam",
                                    extra_rnames_field="rnames", gid_suffix="", inconsistent_traversal_strategy="skip"):
    result = []
    reads_to_adjacencies_ids = defaultdict(set)
    adjacencies_by_aids = {}
    for adj in adjacencies:
        aid = adj.stable_id_non_phased
        read_names = adj.extra.get(extra_rnames_field, "").split(",")
        if len(read_names) == 1 and len(read_names[0]) == 0:
            continue
        for read_name in read_names:
            reads_to_adjacencies_ids[read_name].add(aid)
        adjacencies_by_aids[aid] = adj
    mode = get_mode_str(format=alignment_format, input=True)
    current_read_name = None
    current_entries = []
    cnt = 0
    alignment_1k_counter = 0
    with pysam.AlignmentFile(alignment_file_name, mode) as i_stream:
        if "SO:queryname" not in i_stream.text:
            raise ValueError("Input alignment file needs to be sorted by read (i.e., query) name. It is not.")
        for alignment_cnt, entry in enumerate(i_stream):
            if alignment_cnt / 1000 >= alignment_1k_counter:
                alignment_1k_counter += 1
            if entry.qname != current_read_name:
                if len(current_entries) > 0 and current_read_name in reads_to_adjacencies_ids:
                    adjacencies = [adjacencies_by_aids[aid] for aid in reads_to_adjacencies_ids[current_read_name]]
                    groups = get_labeling_groups(read_alignments=current_entries, read_adjacencies=adjacencies, strategy=inconsistent_traversal_strategy)
                    for group in groups:
                        gid = str(cnt)
                        if len(gid_suffix) > 0:
                            gid += "_" + gid_suffix
                        group.gid = gid
                        cnt += 1
                    result.extend(groups)
                current_read_name = entry.qname
                current_entries = [entry]
            else:
                current_entries.append(entry)
        if len(current_entries) > 0 and current_read_name in reads_to_adjacencies_ids:
            adjacencies = [adjacencies_by_aids[aid] for aid in reads_to_adjacencies_ids[current_read_name]]
            groups = get_labeling_groups(read_alignments=current_entries, read_adjacencies=adjacencies, strategy=inconsistent_traversal_strategy)
            for group in groups:
                gid = str(cnt)
                if len(gid_suffix) > 0:
                    gid += "_" + gid_suffix
                group.gid = gid
                cnt += 1
            result.extend(groups)
    return result
