from collections import defaultdict
from enum import Enum

from core.io import EXTERNAL_NA_ID
from dassp.core.structures import PositionCluster, sorted_segments_donot_overlap


class ClusteringStrategy(Enum):
    IterativeSlider = 0


def cluster_positions(positions, sigma=25, strategy=ClusteringStrategy.IterativeSlider):
    if strategy == ClusteringStrategy.IterativeSlider:
        return cluster_positions_iterative_slider(positions=positions, sigma=sigma)
    return []


def cluster_positions_iterative_slider(positions, sigma):
    per_chrom_positions = defaultdict(list)
    for p in positions:
        per_chrom_positions[p.chromosome].append(p)
    for chromosome in per_chrom_positions.keys():
        per_chrom_positions[chromosome] = sorted(per_chrom_positions[chromosome])
    clusters = []
    for chromosome in sorted(per_chrom_positions.keys(), key=lambda entry: (-len(entry), entry)):  # sorting by chromosome name length first, and then by actual name
        positions = per_chrom_positions[chromosome]
        current_position = positions[0]
        current_cluster_pl = [current_position]
        for p in positions[1:]:
            if p.coordinate - current_position.coordinate <= sigma:
                current_cluster_pl.append(p)
            else:
                pc = PositionCluster(positions=current_cluster_pl)
                clusters.append(pc)
                current_cluster_pl = [p]
            current_position = p
        pc = PositionCluster(positions=current_cluster_pl)
        clusters.append(pc)
    return clusters


def positions_aligned(segments_positions, other_positions):
    segments_positions_ids = {p.stable_id_non_hap for p in segments_positions}
    other_positions_ids = {p.stable_id_non_hap for p in other_positions}
    result = other_positions_ids <= segments_positions_ids
    return result


def nas_groups_aligned(nas_groups_ids, nas):
    nas_ids = {na.extre[EXTERNAL_NA_ID] for na in nas}
    for group in nas_groups_ids:
        for external_na_id in group.adjacencies:
            if external_na_id not in nas_ids:
                return False
    return True
