from collections import defaultdict
from copy import deepcopy

from rck.core.structures import Strand, Position, sorted_segments_donot_overlap


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
            elif current_position.coordinate <= current_fragment.end_position:
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
            elif current_position.coordinate > current_fragment.end_position:
                new_to_old.append(current_segment.stable_id_non_hap)
                current_fragment_id = current_fragment.stable_id_non_hap
                for sid in new_to_old:
                    segments_ids_mapping[current_fragment_id].append(sid)
                    segments_ids_mapping[sid].append(current_fragment_id)
                current_fragment = next(chr_fragments, None)
                current_segment = deepcopy(current_fragment)
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




