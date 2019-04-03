from collections import defaultdict
from copy import deepcopy

from rck.core.io import COPY_NUMBER
from rck.core.structures import Strand, Position, sorted_segments_donot_overlap, Haplotype, SegmentCopyNumberProfile


def refined_segments(segments, additional_positions=None, additional_positions_by_chrs=None):
    fragments = deepcopy(segments)
    if additional_positions is None:
        additional_positions = []
    if additional_positions_by_chrs is None:
        additional_positions_by_chrs = defaultdict(list)
    refined_segments = []
    segments_ids_mapping = defaultdict(list)
    source_fragments_by_chrs = defaultdict(list)
    for fragment in fragments:
        source_fragments_by_chrs[fragment.chromosome].append(fragment)
    for position in additional_positions:
        additional_positions_by_chrs[position.chromosome].append(position)
    for chr_name in list(source_fragments_by_chrs.keys()):
        source_fragments_by_chrs[chr_name] = sorted(source_fragments_by_chrs[chr_name], key=lambda s: (s.start_coordinate, s.end_coordinate))
        if not sorted_segments_donot_overlap(segments=source_fragments_by_chrs[chr_name]):
            raise ValueError("Some segments overlap on chromosome {chr_name}.".format(chr_name=chr_name))
    for chr_name in list(additional_positions_by_chrs.keys()):
        additional_positions_by_chrs[chr_name] = sorted(additional_positions_by_chrs[chr_name], key=lambda p: (p.coordinate, p.strand))
    for chr_name in source_fragments_by_chrs:
        chr_fragments = iter(source_fragments_by_chrs[chr_name])
        chr_positions = iter(additional_positions_by_chrs[chr_name])
        current_fragment = next(chr_fragments, None)
        current_position = next(chr_positions, None)
        current_segment = deepcopy(current_fragment)
        refined_segments.append(current_segment)
        new_to_old = []
        while current_fragment is not None and current_position is not None:
            if current_position.coordinate < current_fragment.start_coordinate:
                current_position = next(chr_positions, None)
            elif current_position.coordinate == current_segment.start_coordinate and current_position.strand == Strand.REVERSE:
                current_position = next(chr_positions, None)
            elif current_position.coordinate == current_segment.end_coordinate and current_position.strand == Strand.FORWARD:
                current_position = next(chr_positions, None)
            elif current_position.coordinate <= current_fragment.end_coordinate:
                if current_position.strand == Strand.FORWARD:
                    left_partition_coordinate = current_position.coordinate
                else:
                    left_partition_coordinate = current_position.coordinate - 1
                right_partition_coordinate = left_partition_coordinate + 1

                new_end_position = Position(chromosome=chr_name, coordinate=left_partition_coordinate, strand=Strand.FORWARD)
                new_start_position = Position(chromosome=chr_name, coordinate=right_partition_coordinate, strand=Strand.REVERSE)
                current_segment.end_position = new_end_position
                new_to_old.append(current_segment.stable_id_non_hap)
                current_segment = deepcopy(current_fragment)
                current_segment.start_position = new_start_position
                refined_segments.append(current_segment)
                current_position = next(chr_positions, None)
            elif current_position.coordinate > current_fragment.end_coordinate:
                new_to_old.append(current_segment.stable_id_non_hap)
                current_fragment_id = current_fragment.stable_id_non_hap
                for sid in new_to_old:
                    segments_ids_mapping[current_fragment_id].append(sid)
                    segments_ids_mapping[sid].append(current_fragment_id)
                current_fragment = next(chr_fragments, None)
                current_segment = deepcopy(current_fragment)
                if current_fragment is not None:
                    refined_segments.append(current_segment)
                new_to_old = []
            else:
                raise ValueError("Something went wrong")
        if current_fragment is not None and current_segment is not None:
            new_to_old.append(current_segment.stable_id_non_hap)
            current_fragment_id = current_fragment.stable_id_non_hap
            for sid in new_to_old:
                segments_ids_mapping[current_fragment_id].append(sid)
                segments_ids_mapping[sid].append(current_fragment_id)
        current_fragment = next(chr_fragments, None)
        current_segment = deepcopy(current_fragment)
        while current_fragment is not None:
            sid = current_segment.stable_id_non_hap
            segments_ids_mapping[sid].append(sid)
            refined_segments.append(current_segment)
            current_fragment = next(chr_fragments, None)
            current_segment = deepcopy(current_fragment)
    return refined_segments, segments_ids_mapping


def lies_within(segment, segments, fully=False):
    for s in segments:
        if fully:
            if s.start_coordinate <= segment.start_coordinate and segment.end_coordinate <= s.end_coordinate:
                return True
        else:
            if s.start_coordinate <= segment.start_coordinate <= s.end_coordinate or \
                    s.start_coordinate <= segment.end_coordinate <= s.end_coordinate or \
                    segment.start_coordinate <= s.start_coordinate <= segment.end_coordinate or \
                    segment.start_coordinate <= s.end_coordinate <= segment.end_coordinate:
                return True
    return False


def filter_segments_by_chromosomal_regions(segments, include=None, exclude=None, include_full=True):
    if include is None:
        include = []
    if exclude is None:
        exclude = []
    include_by_chr = defaultdict(list)
    for segment in include:
        include_by_chr[segment.chromosome.lower()].append(segment)
    exclude_by_chr = defaultdict(list)
    for segment in exclude:
        exclude_by_chr[segment.chromosome.lower()].append(segment)
    for chr in list(include_by_chr.keys()):
        include_by_chr[chr] = sorted(include_by_chr[chr], key=lambda s: (s.start_coordinate, s.end_coordinate))
    for chr in list(exclude_by_chr.keys()):
        exclude_by_chr[chr] = sorted(exclude_by_chr[chr], key=lambda s: (s.start_coordinate, s.end_coordinate))
    include_by_chr = dict(include_by_chr)
    exclude_by_chr = dict(exclude_by_chr)
    for segment in segments:
        retain = True
        chromosome = segment.chromosome
        include_segments_chr = include_by_chr.get(chromosome, [])
        if len(include_segments_chr) != 0:
            retain = lies_within(segment=segment, segments=include_segments_chr, fully=include_full)
        elif len(include_by_chr) > 0:
            retain = False
        if not retain:
            continue
        exclude_segments_chr = exclude_by_chr.get(chromosome, [])
        if len(exclude_segments_chr) != 0:
            retain = not lies_within(segment=segment, segments=exclude_segments_chr, fully=include_full)
        if retain:
            yield segment


def iter_haploid_segments(segments, copy=True):
    for segment in segments:
        result = segment
        if copy:
            result = deepcopy(segment)
        if COPY_NUMBER in result.extra:
            cn_entry = result.extra[COPY_NUMBER]
            result_cn_entry = {}
            for clone_id, cn_dict in cn_entry.items():
                new_cn_dict = {Haplotype.A: sum(cn_dict.values())}
                result_cn_entry[clone_id] = new_cn_dict
            result.extra[COPY_NUMBER] = result_cn_entry
        yield result


def haploid_segments(segments, copy=True):
    return list(iter_haploid_segments(segments=segments, copy=copy))


def get_haploid_scnt(segments, scnt):
    result = {clone_id: SegmentCopyNumberProfile() for clone_id in sorted(scnt.keys())}
    for segment in segments:
        sid = segment.stable_id_non_hap
        for clone_id in result.keys():
            result_scnp: SegmentCopyNumberProfile = result[clone_id]
            source_scnp: SegmentCopyNumberProfile = scnt[clone_id]
            total_cn = source_scnp.get_combined_cn(sid=sid)
            result_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=total_cn)
    return result
