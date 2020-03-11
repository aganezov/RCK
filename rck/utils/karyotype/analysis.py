from collections import defaultdict

from rck.core.io import FALSE_POSITIVE, AG_LABELING
from rck.core.structures import Haplotype, CNBoundaries, AdjacencyGroupType, Phasing


def scnb_violations(scnt, scnb, segments_syncs, segments=None, clone_ids=None, short_circuit=False):
    if clone_ids is None:
        clone_ids = set(scnt.keys()) & set(scnb.keys())
    result = defaultdict(list)
    for clone_id in clone_ids:
        scnp = scnt[clone_id]
        scnbp = scnb[clone_id]
        sids = [s.stable_id_non_hap for s in segments] if segments is not None else scnp.records.keys()
        for sid in sids:
            sync_indicator = segments_syncs[sid]
            cna, cnb = scnp.get_cn(sid=sid, haplotype=Haplotype.A), scnp.get_cn(sid=sid, haplotype=Haplotype.B)
            lower_a = scnbp.get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.LOWER)
            lower_b = scnbp.get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.LOWER)
            upper_a = scnbp.get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.UPPER)
            upper_b = scnbp.get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.UPPER)
            if sync_indicator == 1:
                if (not lower_a <= cna <= upper_a) or (not lower_b <= cnb <= upper_b):
                    result[clone_id].append(sid)
            else:
                if (not lower_b <= cna <= upper_b) or (not lower_a <= cnb <= upper_a):
                    result[clone_id].append(sid)
        if short_circuit and len(result[clone_id]) > 0:
            return result
    return result


def unique_realization_violations(adjacencies, acnt):
    pass


def adjacency_groups_molecule_violations(groups, acnt, clone_ids=None, skip_missing_fp=True, short_circuit=False):
    result = []
    if clone_ids is None:
        clone_ids = sorted(acnt.keys())
    for group in filter(lambda ag: ag.group_type == AdjacencyGroupType.MOLECULE, groups):
        group_fp = group.extra.get(FALSE_POSITIVE, None)
        if group_fp is None:
            if not skip_missing_fp:
                result.append(group)
                if short_circuit:
                    return result
            continue
        group_is_good = False
        for clone_id in clone_ids:
            acnp = acnt[clone_id]
            clone_present = acnp.haploid_adjacencies_present(adjacencies=group.adjacencies)
            inferred_fp = 1 - (len(clone_present) * 1.0 / len(group.adjacencies))
            group_is_good |= inferred_fp <= group_fp
        if not group_is_good:
            result.append(group)
        if short_circuit and len(result) > 0:
            return result
    return result


def adjacency_groups_general_violations(groups, acnt, clone_ids=None, skip_missing_fp=True, short_circuit=False):
    result = []
    if clone_ids is None:
        clone_ids = sorted(acnt.keys())
    for group in filter(lambda ag: ag.group_type == AdjacencyGroupType.GENERAL, groups):
        group_fp = group.extra.get(FALSE_POSITIVE, None)
        if group_fp is None:
            if not skip_missing_fp:
                result.append(group)
                if short_circuit:
                    return result
            continue
        total_present = set()
        for clone_id in clone_ids:
            acnp = acnt[clone_id]
            clone_specific = acnp.haploid_adjacencies_present(adjacencies=group.adjacencies)
            for adjacency in clone_specific:
                total_present.add(adjacency.stable_id_non_phased)
        inferred_fp = 1 - (len(total_present) * 1.0 / len(group.adjacencies_ids))
        if inferred_fp > group_fp:
            result.append(group)
        if short_circuit and len(result) > 0:
            return result
    return result


def adjacency_groups_labeling_violations(groups, acnt, clone_ids=None, short_circuit=False):
    if clone_ids is None:
        clone_ids = sorted(acnt.keys())
    result = []
    for group in filter(lambda ag: ag.group_type == AdjacencyGroupType.LABELING, groups):
        adjacencies = group.adjacencies
        indexes = group.extra[AG_LABELING]
        haplotype_specific_presence = defaultdict(list)
        for hap in [Haplotype.A, Haplotype.B]:
            for adjacency, index in zip(adjacencies, indexes):
                phasings = [Phasing.AA] if hap == Haplotype.A else [Phasing.BB]
                local_phasings = [Phasing.AB, Phasing.BA] if hap == Haplotype.A else [Phasing.BA, Phasing.AB]
                phasings.append(local_phasings[index])
                aid = adjacency.stable_id_non_phased
                cn = 0
                for clone_id in clone_ids:
                    acnp = acnt[clone_id]
                    for ph in phasings:
                        cn += acnp.get_cn(aid=aid, phasing=ph)
                haplotype_specific_presence[hap].append(cn != 0)
        haplotype_specific_presence = dict(haplotype_specific_presence)
        for hap in [Haplotype.A, Haplotype.B]:
            haplotype_specific_presence[hap] = any(haplotype_specific_presence[hap])
        if sum(haplotype_specific_presence.values()) > 1:
            result.append(group)
            if short_circuit:
                return result
    return result


def nas_fp_violations(acnt, fp, adjacencies=None):
    pass
