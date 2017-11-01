from collections import defaultdict
from enum import Enum

from dassp.core.structures import PositionCluster


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
    for chromosome in sorted(per_chrom_positions.keys(), key=lambda entry: (-len(entry), entry)):   # sorting by chromosome name length first, and then by actual name
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
