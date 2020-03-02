from collections import defaultdict
from copy import deepcopy
from typing import Dict, List, Set, Optional, Tuple

import numpy

from rck.core.structures import SegmentCopyNumberProfile, AdjacencyCopyNumberProfile, Adjacency, Haplotype, AdjacencyType, Segment, Strand, Position, \
    refined_scnt_with_adjacencies_and_telomeres


def get_flipped_scnt(segments: List[Segment], scnt: Dict[str, SegmentCopyNumberProfile], flip_prob: float = 0.5) -> Dict[str, SegmentCopyNumberProfile]:
    result = {clone_id: SegmentCopyNumberProfile() for clone_id in scnt.keys()}
    for segment in segments:
        sid = segment.stable_id_non_hap
        flipping = numpy.random.choice([True, False], size=1, p=[flip_prob, 1-flip_prob])[0]
        hap_a, hap_b = (Haplotype.B, Haplotype.A) if flipping else (Haplotype.A, Haplotype.B)
        for clone_id in result.keys():
            r_scnp: SegmentCopyNumberProfile = result[clone_id]
            s_scnp: SegmentCopyNumberProfile = scnt[clone_id]
            r_scnp.set_cn_record(sid=sid, hap=hap_a, cn=s_scnp.get_cn(sid=sid, haplotype=Haplotype.A, default=-1))
            r_scnp.set_cn_record(sid=sid, hap=hap_b, cn=s_scnp.get_cn(sid=sid, haplotype=Haplotype.B, default=-1))
    return result


def get_noisy_scnt(segments: List[Segment], scnt: Dict[str, SegmentCopyNumberProfile], chunk_size: int = 50000) -> Tuple[List[Segment], Dict[str, SegmentCopyNumberProfile]]:
    chrs = {s.chromosome for s in segments}
    chrs_min = {chrom: 3000000000 for chrom in chrs}
    chrx_max = {chrom: 0 for chrom in chrs}
    for s in segments:
        if s.start_coordinate < chrs_min[s.chromosome]:
            chrs_min[s.chromosome] = s.start_coordinate
        if s.end_coordinate > chrx_max[s.chromosome]:
            chrx_max[s.chromosome] = s.end_coordinate
    chunk_segments = []
    for chrom in chrs:
        min_c = chrs_min[chrom]
        max_c = chrx_max[chrom]
        current_left = min_c
        for boundary in range(min_c + chunk_size, max_c, chunk_size):
            s = Segment.from_chromosome_coordinates(chromosome=chrom, start=current_left + 1, end=boundary)
            current_left = s.end_coordinate
            chunk_segments.append(s)
    all_positions = set()
    for tmp_segments in [segments, chunk_segments]:
        for s in tmp_segments:
            all_positions.add(s.start_position)
            all_positions.add(s.end_position)
    refined_segments, refined_scnt, segments_ids_mapping = refined_scnt_with_adjacencies_and_telomeres(segments=segments, scnt=scnt, telomere_positions=all_positions)
    refined_segments_by_chrs = defaultdict(list)
    for s in refined_segments:
        refined_segments_by_chrs[s.chromosome].append(s)
    for chrom in list(refined_segments_by_chrs.keys()):
        refined_segments_by_chrs[chrom] = sorted(refined_segments_by_chrs[chrom], key=lambda s: (s.start_coordinate, s.end_coordinate))
    result = {clone_id: SegmentCopyNumberProfile() for clone_id in scnt.keys()}
    considered_rs_segments = set()
    for s in chunk_segments:
        spanning_refined_segments: List[Segment] = []
        for rs in refined_segments_by_chrs[s.chromosome]:
            if rs.start_coordinate >= s.start_coordinate:
                if rs.end_coordinate <= s.end_coordinate:
                    assert rs not in considered_rs_segments
                    spanning_refined_segments.append(rs)
                    considered_rs_segments.add(rs)
                else:
                    break
        for clone_id in result:
            source_scnp: SegmentCopyNumberProfile = refined_scnt[clone_id]
            target_scnp: SegmentCopyNumberProfile = result[clone_id]
            for hap in [Haplotype.A, Haplotype.B]:
                cn = 0
                for rs in spanning_refined_segments:
                    rs_cn = source_scnp.get_cn(sid=rs.stable_id_non_hap, haplotype=hap, default=-1)
                    assert rs_cn != -1
                    cn += (rs_cn * rs.length)
                fractional_cn = 1.0 * cn / sum(rs.length for rs in spanning_refined_segments)
                cn = round(fractional_cn)
                target_scnp.set_cn_record(sid=s.stable_id_non_hap, hap=hap, cn=cn)
    return chunk_segments, result


def get_unlabeled_adjacencies(adjacencies: List[Adjacency], acnt: Dict[str, AdjacencyCopyNumberProfile],
                              novel_only: bool = True, clones: Optional[Set[str]] = None, present_only: bool = True) -> List[Adjacency]:
    result = adjacencies
    if novel_only:
        result = [a for a in result if a.adjacency_type == AdjacencyType.NOVEL]
    if clones is None:
        clones = set(acnt.keys())
    if present_only:
        tmp_result = []
        for adj in result:
            aid = adj.stable_id_non_phased
            for clone in clones:
                if acnt[clone].get_combined_cn(aid=aid, default=-1) > 0:
                    tmp_result.append(adj)
                    break
        result = tmp_result
    return result


def get_noisy_adjacencies(adjacencies: List[Adjacency], coordinate_noise_prob: bool = 0.5, coordinate_noise_max: int = 50,
                          fp_rate: float = 0.00) -> List[Adjacency]:
    tmp_result = []
    for adj in adjacencies:
        adj = deepcopy(adj)
        adj.extra["or_coordinates"] = (adj.position1.coordinate, adj.position2.coordinate)
        if adj.position1.chromosome != adj.position2.chromosome:
            tmp_result.append(adj)
            continue
        add_noise = numpy.random.choice([True, False], size=1, p=[coordinate_noise_prob, 1 - coordinate_noise_prob])[0]
        if add_noise:
            adj_length = adj.distance_non_hap
            left_shift = numpy.random.randint(low=-1 * coordinate_noise_max - 1, high=coordinate_noise_max + 1)
            right_shift = numpy.random.randint(low=-1 * coordinate_noise_max - 1, high=coordinate_noise_max + 1)
            if left_shift - right_shift < adj_length:
                adj.position1.coordinate += left_shift
                adj.position2.coordinate += right_shift
        tmp_result.append(adj)
    result = tmp_result
    fp_cnt = int(fp_rate * len(result))
    chrs = {adj.position1.chromosome for adj in result}
    chrs_min = {chrom: 3000000000 for chrom in chrs}
    chrs_max = {chrom: 0 for chrom in chrs}
    for adj in adjacencies:
        for p in [adj.position1, adj.position2]:
            if p.coordinate < chrs_min[p.chromosome]:
                chrs_min[p.chromosome] = p.coordinate
            if p.coordinate > chrs_max[p.chromosome]:
                chrs_max[p.chromosome] = p.coordinate
    for _ in range(fp_cnt):
        inter = numpy.random.choice([True, False], size=1, p=[0.01, 0.99])[0]
        chr1 = numpy.random.choice(list(chrs), size=1)[0]
        chr2 = chr1
        coordinate1 = numpy.random.randint(low=chrs_min[chr1], high=chrs_max[chr1])
        strand1 = numpy.random.choice([Strand.FORWARD, Strand.REVERSE], size=1, p=[0.5, 0.5])[0]
        strand2 = numpy.random.choice([Strand.FORWARD, Strand.REVERSE], size=1, p=[0.5, 0.5])[0]
        if inter:
            chr2 = numpy.random.choice(list(chrs), size=1)[0]
            coordinate2 = numpy.random.randint(low=chrs_min[chr2], high=chrs_max[chr2])
        else:
            coordinate2 = numpy.random.randint(low=coordinate1, high=chrs_max[chr1])
        p1 = Position(chromosome=chr1, coordinate=coordinate1, strand=strand1)
        p2 = Position(chromosome=chr2, coordinate=coordinate2, strand=strand2)
        adj = Adjacency(position1=p1, position2=p2, adjacency_type=AdjacencyType.NOVEL, extra={"fake": True})
        result.append(adj)
    return result
