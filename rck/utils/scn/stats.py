import itertools
from collections import defaultdict
from typing import Dict, Tuple, List

from rck.core.structures import refined_scnt_with_adjacencies_and_telomeres, refined_scnt, cn_distance_inter_scnt, check_and_fill_segments_to_fragments, \
    segments_to_fragments_reversed, SegmentCopyNumberProfile, Haplotype, Segment

T1 = "t1"
T2 = "t2"


class CloneCollectionCNInstance(object):
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
            "end": max(positions, key=lambda p: p.coordinate)
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
    assert set(segments1) == set(segments2)
    assert len(segments1) == len(set(segments1))
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
            case = CloneCollectionCNInstance(instances1=clone_ids1_instances, instances2=clone_ids2_instances,
                                             mapping1=clone_ids1_mapping, mapping2=clone_ids2_mapping)
            result[case] = clone_specific_distances
    return result


def cn_changed_results(tensor1, tensor2, segments, segments_to_fragments=None, check_clone_ids_match=True) -> Dict[str, Dict[str, SegmentCopyNumberProfile]]:
    if check_clone_ids_match and set(tensor1.keys()) != set(tensor2.keys()):
        raise Exception()
    assert len(sorted(tensor1.keys())) == len(sorted(tensor2.keys()))
    clone_ids_source = sorted(tensor1.keys())
    clone_ids_target = sorted(tensor2.keys())
    sids_to_segments = {segment.stable_id_non_hap: segment for segment in segments}
    segments_to_fragments = check_and_fill_segments_to_fragments(segments=segments, segments_to_fragments=segments_to_fragments)
    fragments_to_segments = segments_to_fragments_reversed(segments_to_fragments=segments_to_fragments)
    result = {clone_id: {T1: SegmentCopyNumberProfile(), T2: SegmentCopyNumberProfile()} for clone_id in clone_ids_target}
    for fid, sids in fragments_to_segments.items():
        fid_value_ab = 0
        fid_value_ba = 0
        for cid_source, cid_target in zip(clone_ids_source, clone_ids_target):
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                cn_dif_B = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                value = ((cn_dif_A + cn_dif_B) * segment_length)
                fid_value_ab += value
        for cid_source, cid_target in zip(clone_ids_source, clone_ids_target):
            for sid in sids:
                segment_length = sids_to_segments[sid].length
                cn_dif_A = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.A, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                cn_dif_B = abs(tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.B, default=0) - tensor2[cid_target].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                value = ((cn_dif_A + cn_dif_B) * segment_length)
                fid_value_ba += value
        if fid_value_ab <= fid_value_ba:
            hap_a = Haplotype.A
            hap_b = Haplotype.B
        else:
            hap_a = Haplotype.B
            hap_b = Haplotype.A
        for cid_source, cid_target in zip(clone_ids_source, clone_ids_target):
            for sid in sids:
                t1_scnp: SegmentCopyNumberProfile = result[cid_source][T1]
                t2_scnp: SegmentCopyNumberProfile = result[cid_target][T2]
                t1_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.A, default=0))
                t1_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=tensor1[cid_source].get_cn(sid=sid, haplotype=Haplotype.B, default=0))
                t2_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=tensor2[cid_source].get_cn(sid=sid, haplotype=hap_a, default=0))
                t2_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=tensor2[cid_source].get_cn(sid=sid, haplotype=hap_b, default=0))
    return result


def cn_change(segments1, scnt1, segments2, scnt2, both_haplotype_specific=False) -> Dict[CloneCollectionCNInstance, Tuple]:
    positions_by_chr = defaultdict(set)
    for segments in [segments1, segments2]:
        for segment in segments:
            positions_by_chr[segment.chromosome].add(segment.start_position)
            positions_by_chr[segment.chromosome].add(segment.end_position)
    outermost_positions_per_chromosomes = {}
    for chr_name, positions in positions_by_chr.items():
        outermost_positions_per_chromosomes[chr_name] = {
            "start": min(positions, key=lambda p: p.coordinate),
            "end": max(positions, key=lambda p: p.coordinate)
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
    result: Dict[CloneCollectionCNInstance, Tuple] = {}
    for clone_ids1_instances in itertools.combinations(clone_ids1, matching_clone_ids_cnt):
        for clone_ids2_instances in itertools.permutations(clone_ids2, matching_clone_ids_cnt):
            clone_ids1_mapping = {str(cnt): clone_id for cnt, clone_id in enumerate(clone_ids1_instances)}
            clone_ids2_mapping = {str(cnt): clone_id for cnt, clone_id in enumerate(clone_ids2_instances)}
            tmp_scnt1 = {key: scnt1[value] for key, value in clone_ids1_mapping.items()}
            tmp_scnt2 = {key: scnt2[value] for key, value in clone_ids2_mapping.items()}
            clone_specific_cn_change_result = cn_changed_results(tensor1=tmp_scnt1, tensor2=tmp_scnt2, segments=segments1)
            clone_specific_distances = cn_distance_inter_scnt(tensor1=tmp_scnt1, tensor2=tmp_scnt2, segments=segments1)
            case = CloneCollectionCNInstance(instances1=clone_ids1_instances, instances2=clone_ids2_instances,
                                             mapping1=clone_ids1_mapping, mapping2=clone_ids2_mapping)
            result[case] = (segments1, clone_specific_distances, clone_specific_cn_change_result)
    return result


def changed_haploid_segments(scnp1: SegmentCopyNumberProfile, scnp2: SegmentCopyNumberProfile, segments: List[Segment]) -> List[Segment]:
    result = []
    for segment in segments:
        sid = segment.stable_id_non_hap
        source_cn = scnp1.get_combined_cn(sid=sid)
        target_cn = scnp2.get_combined_cn(sid=sid, default=0)
        if source_cn != target_cn:
            result.append(segment)
    return result


def changed_segments(scnp1: SegmentCopyNumberProfile, scnp2: SegmentCopyNumberProfile, segments: List[Segment]) -> List[Segment]:
    result = []
    for segment in segments:
        sid = segment.stable_id_non_hap
        for hap in [Haplotype.A, Haplotype.B]:
            source_cn = scnp1.get_cn(sid=sid, haplotype=hap, default=-1)
            target_cn = scnp2.get_cn(sid=sid, haplotype=hap, default=-1)
            assert -1 not in {source_cn, target_cn}
            if source_cn != target_cn:
                result.append(segment)
    return result


def changed_ampl_segments(scnp1: SegmentCopyNumberProfile, scnp2: SegmentCopyNumberProfile, segments: List[Segment]) -> List[Segment]:
    result = []
    for segment in segments:
        sid = segment.stable_id_non_hap
        for hap in [Haplotype.A, Haplotype.B]:
            source_cn = scnp1.get_cn(sid=sid, haplotype=hap, default=-1)
            target_cn = scnp2.get_cn(sid=sid, haplotype=hap, default=-1)
            assert -1 not in {source_cn, target_cn}
            if source_cn > 1 and target_cn <= 1:
                result.append(segment)
    return result


def changed_neutral_segments(scnp1: SegmentCopyNumberProfile, scnp2: SegmentCopyNumberProfile, segments: List[Segment]) -> List[Segment]:
    result = []
    for segment in segments:
        sid = segment.stable_id_non_hap
        for hap in [Haplotype.A, Haplotype.B]:
            source_cn = scnp1.get_cn(sid=sid, haplotype=hap, default=-1)
            target_cn = scnp2.get_cn(sid=sid, haplotype=hap, default=-1)
            assert -1 not in {source_cn, target_cn}
            if source_cn == 1 and target_cn != 1:
                result.append(segment)
    return result


def changed_loss_segments(scnp1: SegmentCopyNumberProfile, scnp2: SegmentCopyNumberProfile, segments: List[Segment]) -> List[Segment]:
    result = []
    for segment in segments:
        sid = segment.stable_id_non_hap
        for hap in [Haplotype.A, Haplotype.B]:
            source_cn = scnp1.get_cn(sid=sid, haplotype=hap, default=-1)
            target_cn = scnp2.get_cn(sid=sid, haplotype=hap, default=-1)
            assert -1 not in {source_cn, target_cn}
            if source_cn == 0 and target_cn > 0:
                result.append(segment)
    return result

