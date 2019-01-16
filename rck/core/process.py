from enum import Enum

from rck.core.io import EXTERNAL_NA_ID


class ClusteringStrategy(Enum):
    IterativeSlider = 0


def positions_aligned(segments_positions, other_positions):
    segments_positions_ids = {p.stable_id_non_hap for p in segments_positions}
    other_positions_ids = {p.stable_id_non_hap for p in other_positions}
    result = other_positions_ids <= segments_positions_ids
    return result


def adj_groups_concur(adj_groups, adjacencies):
    nas_ids = {na.extra.get(EXTERNAL_NA_ID, na.idx) for na in adjacencies}
    for group in adj_groups:
        for aid in group.adjacencies_ids:
            if aid not in nas_ids:
                return False
    return True
