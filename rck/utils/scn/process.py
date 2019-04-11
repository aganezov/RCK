from collections import defaultdict
from copy import deepcopy

from rck.core.io import COPY_NUMBER, stringify_adjacency_cn_entry
from rck.core.structures import Strand, Position, sorted_segments_donot_overlap, Haplotype, SegmentCopyNumberProfile, Segment, refined_scnt_with_adjacencies_and_telomeres
from rck.utils.adj.process import REMOVE


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


def filter_segments_by_chromosomal_regions(segments, include=None, exclude=None, include_full=True, exclude_full=False):
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
            retain = not lies_within(segment=segment, segments=exclude_segments_chr, fully=exclude_full)
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


def filter_segments_by_extra(segments, keep_extra_field=None, keep_extra_field_missing_strategy=None, remove_extra_field=None, remove_extra_field_missing_strategy=None):
    for segment in segments:
        if keep_extra_field is not None and len(keep_extra_field) > 0:
            keep = False
            for field, regex_list in keep_extra_field.items():
                for regex in regex_list:
                    if field in segment.extra:
                        data = segment.extra[field]
                        if not isinstance(data, str):
                            if field == COPY_NUMBER:
                                data = stringify_adjacency_cn_entry(entry=data)
                            else:
                                data = str(data)
                        if regex.search(data) is not None:
                            keep = True
                    elif keep_extra_field_missing_strategy == REMOVE:
                        keep = keep or False
            if not keep:
                continue

        if remove_extra_field is not None and len(remove_extra_field) > 0:
            keep = True
            for field, regex_list in remove_extra_field.items():
                for regex in regex_list:
                    if field in segment.extra:
                        data = segment.extra[field]
                        if not isinstance(data, str):
                            if field == COPY_NUMBER:
                                data = stringify_adjacency_cn_entry(entry=data)
                            else:
                                data = str(data)
                        if regex.search(data) is not None:
                            keep = False
                            break
                    elif remove_extra_field_missing_strategy == REMOVE:
                        keep = False
                        break
            if not keep:
                continue
        yield segment


def filter_segments_by_size(segments, min_size=0, max_size=1000000000):
    for segment in segments:
        size = abs(segment.end_coordinate - segment.start_coordinate)
        if min_size <= size <= max_size:
            yield segment


def segment_within_segment(inner_segment, outer_segment):
    start_in = outer_segment.start_coordinate <= inner_segment.start_coordinate <= outer_segment.end_coordinate
    end_in = outer_segment.start_coordinate <= inner_segment.end_coordinate <= outer_segment.end_coordinate
    return start_in and end_in



def get_circa_segments_cna_fractions(segments, scnt, clone_id, window_size=10000000, chr_sizes=None, cna_type="ampl", haploid=False, inverse=False):
    intra_result = defaultdict(list)
    segments_by_chrs = defaultdict(list)
    for segment in segments:
        segments_by_chrs[segment.chromosome].append(segment)
    for chr_name in sorted(segments_by_chrs.keys()):
        segments_by_chrs[chr_name] = sorted(segments_by_chrs[chr_name], key=lambda s: (s.start_coordinate, s.end_coordinate))
    windows_by_chr = defaultdict(list)
    if chr_sizes is None:
        chr_sizes = {}
    for chr_name in set(chr_sizes.keys()) | set(segments_by_chrs.keys()):
        start = 0
        default = segments_by_chrs[chr_name][-1].end_coordinate if chr_name in segments_by_chrs else 0
        end = chr_sizes.get(chr_name, default)
        windows_boundaries = list(range(start, end, window_size))
        if windows_boundaries[-1] != end:
            windows_boundaries.append(end)
        for lb, rb in zip(windows_boundaries[:-1], windows_boundaries[1:]):
            segment = Segment.from_chromosome_coordinates(chromosome=chr_name, start=lb + 1, end=rb)
            windows_by_chr[chr_name].append(segment)
    windows = []
    for ws in windows_by_chr.values():
        chr_name = ws[0].chromosome
        if chr_name in segments_by_chrs and segments_by_chrs[chr_name][0].start_coordinate > ws[0].start_coordinate:
            segments_by_chrs[chr_name][0].start_position.coordinate = ws[0].start_coordinate
        if chr_name in segments_by_chrs and segments_by_chrs[chr_name][-1].end_coordinate < ws[-1].end_coordinate:
            segments_by_chrs[chr_name][-1].end_position.coordinate = ws[-1].end_coordinate
        for w in ws:
            windows.append(w)
    positions = []
    for w in windows:
        if w.chromosome in segments_by_chrs:
            positions.append(w.start_position)
            positions.append(w.end_position)
    r_segments, r_scnt, _ = refined_scnt_with_adjacencies_and_telomeres(segments=segments, scnt=scnt, telomere_positions=positions)
    for chr_name in windows_by_chr.keys():
        chr_windows = iter(windows_by_chr[chr_name])
        current_window = next(chr_windows, None)
        for segment in r_segments:
            if current_window is None:
                break
            if segment.start_coordinate < current_window.start_coordinate:
                continue
            elif segment_within_segment(inner_segment=segment, outer_segment=current_window):
                intra_result[current_window].append(segment)
            else:
                current_window = next(chr_windows, None)
    scnp: SegmentCopyNumberProfile = r_scnt[clone_id]
    result = {}
    for window, segments in intra_result.items():
        cna_fraction = 0
        for segment in segments:
            length_fraction = segment.length / window.length
            if haploid:
                cns = [scnp.get_combined_cn(sid=segment.stable_id_non_hap)]
                base = 2
            else:
                cns = [scnp.get_cn(sid=segment.stable_id_non_hap, haplotype=Haplotype.A), scnp.get_cn(sid=segment.stable_id_non_hap, haplotype=Haplotype.B)]
                base = 1
            amplified = any(map(lambda cn: cn > base, cns))
            deletions = any(map(lambda cn: cn < base, cns))
            if cna_type == "ampl" and amplified:
                cna_fraction += length_fraction
            elif cna_type == "del" and deletions:
                cna_fraction += length_fraction
        if inverse:
            cna_fraction = 1 - cna_fraction
        result[window] = cna_fraction
    return result
