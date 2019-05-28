import itertools
from collections import defaultdict

from rck.core.structures import refined_scnt_with_adjacencies_and_telomeres, refined_scnt, cn_distance_inter_scnt


class CloneCollectionCNDistanceInstance(object):
    def __init__(self, instances1, instances2, mapping1, mapping2):
        self.instances1 = instances1
        self.instances2 = instances2
        self.mapping1 = mapping1
        self.mapping2 = mapping2

    def __str__(self):
        return f'instance 1: ({",".join(self.instances1)}); instance 2: ({",".join(self.instances2)})'


def cn_distance(segments1, scnt1, segments2, scnt2, both_haplotype_specific=False):
    positions_by_chr = defaultdict(set)
    for segments in [segments1, segments2]:
        for segment in segments:
            positions_by_chr[segment.chromosome].add(segment.start_position)
            positions_by_chr[segment.chromosome].add(segment.end_position)
    outermost_positions_per_chromosomes = {}
    for chr_name, positions in positions_by_chr.items():
        outermost_positions_per_chromosomes[chr_name] = {
            "start": min(positions, key=lambda p: p.coordinate),
            "end":   max(positions, key=lambda p: p.coordinate)
        }
    segments1, scnt1, _ = refined_scnt(segments=segments1, scnt=scnt1, merge_fragments=False, fill_gaps=True, extend_outermost=True,
                                       outermost_positions=outermost_positions_per_chromosomes, outermost_positions_margin=0)
    segments2, scnt2, _ = refined_scnt(segments=segments2, scnt=scnt2, merge_fragments=False, fill_gaps=True, extend_outermost=True,
                                       outermost_positions=outermost_positions_per_chromosomes, outermost_positions_margin=0)
    all_positions = set()
    for segments in [segments1, segments2]:
        for segment in segments:
            all_positions.add(segment.start_position)
            all_positions.add(segment.end_position)
    segments1, scnt1, _ = refined_scnt_with_adjacencies_and_telomeres(segments=segments1, scnt=scnt1, telomere_positions=all_positions)
    segments2, scnt2, _ = refined_scnt_with_adjacencies_and_telomeres(segments=segments2, scnt=scnt2, telomere_positions=all_positions)
    clone_ids1, clone_ids2 = list(set(scnt1.keys())), list(set(scnt2.keys()))
    matching_clone_ids_cnt = min((len(clone_ids1), len(clone_ids2)))
    result = {}
    for clone_ids1_instances in itertools.combinations(clone_ids1, matching_clone_ids_cnt):
        for clone_ids2_instances in itertools.permutations(clone_ids2, matching_clone_ids_cnt):
            clone_ids1_mapping = {str(cnt): clone_id for cnt, clone_id in enumerate(clone_ids1_instances)}
            clone_ids2_mapping = {str(cnt): clone_id for cnt, clone_id in enumerate(clone_ids2_instances)}
            tmp_scnt1 = {key: scnt1[value] for key, value in clone_ids1_mapping.items()}
            tmp_scnt2 = {key: scnt2[value] for key, value in clone_ids2_mapping.items()}
            clone_specific_distances = cn_distance_inter_scnt(tensor1=tmp_scnt1, tensor2=tmp_scnt2, segments=segments1)
            case = CloneCollectionCNDistanceInstance(instances1=clone_ids1_instances, instances2=clone_ids2_instances,
                                                     mapping1=clone_ids1_mapping, mapping2=clone_ids2_mapping)
            result[case] = clone_specific_distances
    return result
