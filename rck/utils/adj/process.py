import argparse
from collections import defaultdict
import re
from copy import deepcopy, copy

from sortedcontainers import SortedList
import statistics

from rck.core.io import EXTERNAL_NA_ID, COPY_NUMBER, stringify_adjacency_cn_entry, ADJACENCY_TYPE
from rck.core.structures import Strand, Adjacency, Position, Phasing, AdjacencyCopyNumberProfile, Segment
from rck.utils.adj.convert import GUNDEM_PER_SAMPLE_SUPPORT

ORIGIN_IDS = "origin_ids"

KEEP = "keep"
REMOVE = "remove"


class Merger(object):
    def __init__(self, origin_ids_field=ORIGIN_IDS):
        self.origin_ids_field = origin_ids_field
        self.adjs_by_ids = {}
        self.merged_position_pairs = set()
        self.ajds_ids_by_positions = defaultdict(set)
        self.data = defaultdict(lambda: {
            Strand.FORWARD: SortedList(),
            Strand.REVERSE: SortedList(),
        })

    def add_adjacency(self, adjacency, force_sv_type=False, force_na_type=False, max_distance=-1):
        aid = adjacency.extra.get(EXTERNAL_NA_ID, adjacency.idx)
        lp, rp = adjacency.position1, adjacency.position2
        chrl, chrr = lp.chromosome, rp.chromosome
        coordl, coordr = lp.coordinate, rp.coordinate
        strandl, strandr = lp.strand, rp.strand
        lp_neighbours = self.get_closest_coordinates(chromosome=chrl, coordinate=coordl, strand=strandl, max_distance=max_distance)
        rp_neighbours = self.get_closest_coordinates(chromosome=chrr, coordinate=coordr, strand=strandr, max_distance=max_distance)
        if lp_neighbours is None or rp_neighbours is None:
            self.__insert_adjacency(adjacency=adjacency)
            self.adjs_by_ids[aid] = adjacency
            self.add_adjs_id_by_position(aid=aid, chromosome=chrl, strand=strandl, coordinate=coordl)
            self.add_adjs_id_by_position(aid=aid, chromosome=chrr, strand=strandr, coordinate=coordr)
            self.merged_position_pairs.add(((chrl, strandl, coordl), (chrr, strandr, coordr)))
            return
        neighbour_pairs = []
        for ln in lp_neighbours:
            for rn in rp_neighbours:
                if chrl == chrr and strandl == strandr and ln >= rn:
                    continue
                neighbour_pairs.append((ln, rn))
        if len(neighbour_pairs) == 0:
            neighbour_pairs = [(-1, -1)]
            fall_through = True
        else:
            fall_through = False
            neighbour_pairs = sorted(neighbour_pairs, key=lambda e: abs(e[0] - coordl) + abs(e[1] - coordr))
        for lpn, rpn in neighbour_pairs:
            if fall_through:
                break
            lpn_adj_ids = self.get_adjs_ids_by_position(chromosome=chrl, strand=strandl, coordinate=lpn)
            rpn_adj_ids = self.get_adjs_ids_by_position(chromosome=chrr, strand=strandr, coordinate=rpn)
            common_adjs_ids = lpn_adj_ids & rpn_adj_ids
            if len(common_adjs_ids) == 0:
                continue
            # here we merge the stuffffff
            common_adjs = [self.adjs_by_ids[common_aid] for common_aid in common_adjs_ids]
            lp_coordinates = []
            rp_coordinates = []
            aid_flipped = {}
            for adj in common_adjs:
                p1_chr, p2_chr = adj.position1.chromosome, adj.position2.chromosome
                p1_coord, p2_coord = adj.position1.coordinate, adj.position2.coordinate
                p1_strand, p2_strand = adj.position1.strand, adj.position2.strand
                if (p1_chr, p1_strand) == (chrl, strandl) and (p2_chr, p2_strand) == (chrr, strandr):
                    aid_flipped[adj.extra.get(EXTERNAL_NA_ID, adj.idx)] = False
                    lp_coordinates.append(p1_coord)
                    rp_coordinates.append(p2_coord)
                elif (p1_chr, p1_strand) == (chrr, strandr) and (p2_chr, p2_strand) == (chrl, strandl):
                    aid_flipped[adj.extra.get(EXTERNAL_NA_ID, adj.idx)] = True
                    lp_coordinates.append(p2_coord)
                    rp_coordinates.append(p1_coord)
                else:
                    raise Exception("Problems in merging")
            lp_coordinates.append(coordl)
            rp_coordinates.append(coordr)
            new_lp_coordinate = int(statistics.mean(lp_coordinates))
            new_rp_coordinate = int(statistics.mean(rp_coordinates))
            if common_adjs_ids == lpn_adj_ids:
                self.data[chrl][strandl].discard(lpn)
            if common_adjs_ids == rpn_adj_ids:
                self.data[chrr][strandr].discard(rpn)
            if new_lp_coordinate not in self.data[chrl][strandl]:
                self.data[chrl][strandl].add(new_lp_coordinate)
            if new_rp_coordinate not in self.data[chrr][strandr]:
                self.data[chrr][strandr].add(new_rp_coordinate)
            assert ((chrl, strandl, lpn), (chrr, strandr, rpn)) in self.merged_position_pairs or ((chrr, strandr, rpn), (chrl, strandl, lpn)) in self.merged_position_pairs
            self.merged_position_pairs.discard(((chrl, strandl, lpn), (chrr, strandr, rpn)))
            self.merged_position_pairs.discard(((chrr, strandr, rpn), (chrl, strandl, lpn)))
            self.merged_position_pairs.add(((chrl, strandl, new_lp_coordinate), (chrr, strandr, new_rp_coordinate)))
            for common_adj_id in common_adjs_ids:
                # l_coord, r_coord = (rpn, lpn) if aid_flipped[common_adj_id] else (lpn, rpn)
                self.delete_adjs_id_by_position(aid=common_adj_id, chromosome=chrl, coordinate=lpn, strand=strandl)
                self.delete_adjs_id_by_position(aid=common_adj_id, chromosome=chrr, coordinate=rpn, strand=strandr)
            self.adjs_by_ids[aid] = adjacency
            for current_aid in sorted(common_adjs_ids) + [aid]:
                self.add_adjs_id_by_position(aid=current_aid, chromosome=chrl, coordinate=new_lp_coordinate, strand=strandl)
                self.add_adjs_id_by_position(aid=current_aid, chromosome=chrr, coordinate=new_rp_coordinate, strand=strandr)
            break
        else:
            self.__insert_adjacency(adjacency=adjacency)
            self.adjs_by_ids[adjacency.extra.get(EXTERNAL_NA_ID, adjacency.idx)] = adjacency
            self.add_adjs_id_by_position(aid=aid, chromosome=chrl, strand=strandl, coordinate=coordl)
            self.add_adjs_id_by_position(aid=aid, chromosome=chrr, strand=strandr, coordinate=coordr)
            self.merged_position_pairs.add(((chrl, strandl, coordl), (chrr, strandr, coordr)))

    def __insert_adjacency(self, adjacency):
        if adjacency.position1.coordinate not in self.data[adjacency.position1.chromosome][adjacency.position1.strand]:
            self.data[adjacency.position1.chromosome][adjacency.position1.strand].add(adjacency.position1.coordinate)
        if adjacency.position2.coordinate not in self.data[adjacency.position2.chromosome][adjacency.position2.strand]:
            self.data[adjacency.position2.chromosome][adjacency.position2.strand].add(adjacency.position2.coordinate)

    def get_closest_coordinates(self, chromosome, coordinate, strand, max_distance=-1):
        coordinates = self.data[chromosome][strand]
        if len(coordinates) == 0:
            return None
        if coordinate in coordinates:
            return [coordinate]
        left_index = coordinates.bisect_left(coordinate)
        neighbours = []
        for i in reversed(range(0, left_index)):
            other_coordinate = coordinates[i]
            if abs(other_coordinate - coordinate) < max_distance:
                neighbours.append(other_coordinate)
            else:
                break
        for i in range(left_index, len(coordinates)):
            other_coordinate = coordinates[i]
            if abs(other_coordinate - coordinate) < max_distance:
                neighbours.append(other_coordinate)
            else:
                break
        if len(neighbours) == 0:
            return None
        return sorted(neighbours, key=lambda e: abs(e - coordinate))

    def get_closest_coordinate_right(self, chromosome, coordinate, strand, max_distance=1):
        coordinates = self.data[chromosome][strand]
        if coordinate in coordinates:
            return coordinate
        right_coordinate = coordinates.bisect_right(coordinate)
        if len(coordinates) == 0 or right_coordinate == len(coordinates):
            return None
        right_neighbour = coordinates[right_coordinate]
        if max_distance >= 0 and abs(right_neighbour - coordinate) > max_distance:
            return None
        return right_neighbour

    def get_adjs_ids_by_position(self, chromosome, coordinate, strand):
        return self.ajds_ids_by_positions.get((chromosome, strand, coordinate), set())

    def add_adjs_id_by_position(self, aid, chromosome, coordinate, strand):
        self.ajds_ids_by_positions[(chromosome, strand, coordinate)].add(aid)

    def delete_adjs_id_by_position(self, aid, chromosome, coordinate, strand):
        if (chromosome, strand, coordinate) in self.ajds_ids_by_positions:
            self.ajds_ids_by_positions[(chromosome, strand, coordinate)].discard(aid)

    def get_merged_adjacencies(self, merged_template="{cnt}_merged"):
        result = []
        processed_ids = set()
        cnt = 0
        for lp, rp in self.merged_position_pairs:
            lp_chr, lp_strand, lp_coord = lp
            rp_chr, rp_strand, rp_coord = rp
            lp_ids = self.get_adjs_ids_by_position(chromosome=lp_chr, strand=lp_strand, coordinate=lp_coord)
            rp_ids = self.get_adjs_ids_by_position(chromosome=rp_chr, strand=rp_strand, coordinate=rp_coord)
            common_ids = lp_ids & rp_ids
            assert len(common_ids) > 0
            assert len(common_ids & processed_ids) == 0
            processed_ids.update(common_ids)
            pos1 = Position(chromosome=lp_chr, coordinate=lp_coord, strand=lp_strand)
            pos2 = Position(chromosome=rp_chr, coordinate=rp_coord, strand=rp_strand)
            adj = Adjacency(position1=pos1, position2=pos2, extra={
                EXTERNAL_NA_ID: merged_template.format(cnt=cnt),
                self.origin_ids_field: ",".join(common_ids),
            })
            result.append(adj)
            cnt += 1
        return result

    @staticmethod
    def adjacencies_mergeable(adjacency1, adjacency2, max_distance=500, separate_distances=True):
        distance11 = Position.non_hap_distance_strand_specific(pos1=adjacency1.position1, pos2=adjacency2.position1)
        distance12 = Position.non_hap_distance_strand_specific(pos1=adjacency1.position2, pos2=adjacency2.position2)
        distance21 = Position.non_hap_distance_strand_specific(pos1=adjacency1.position2, pos2=adjacency2.position1)
        distance22 = Position.non_hap_distance_strand_specific(pos1=adjacency1.position1, pos2=adjacency2.position2)
        if not separate_distances:
            return min(distance11 + distance12, distance21 + distance22) <= max_distance
        case1 = distance11 <= max_distance and distance12 <= max_distance
        case2 = distance21 <= max_distance and distance22 <= max_distance
        return case1 or case2


ANNOTATE_RETAINED_EXTRA_FIELD = "retained_by"


def filter_adjacencies_by_chromosomal_regions(adjacencies, include=None, exclude=None, include_both=True, exclude_both=False, include_spanning=False, exclude_spanning=False,
                                              annotate_retained=False, annotate_retained_extra_field_prefix=None, annotated_retained_segments_extra_field=None,
                                              annotate_short_circ=False):
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
    annotate_retained_extra_field = ANNOTATE_RETAINED_EXTRA_FIELD
    if annotate_retained_extra_field_prefix is not None:
        annotate_retained_extra_field = "_".join([annotate_retained_extra_field_prefix, annotate_retained_extra_field])
    for adj in adjacencies:
        retain = True
        chr1, chr2 = adj.position1.chromosome, adj.position2.chromosome
        chr1, chr2 = chr1.lower(), chr2.lower()
        include_segments_chr1 = include_by_chr.get(chr1, [])
        include_segments_chr2 = include_by_chr.get(chr2, [])
        annotations_segments = []
        if len(include_segments_chr1) != 0 or len(include_segments_chr2) != 0:
            chr1_segments_in = coordinate_in_segments(coordinate=adj.position1.coordinate, segments=include_segments_chr1, short_circ=annotate_short_circ)
            chr2_segments_in = coordinate_in_segments(coordinate=adj.position2.coordinate, segments=include_segments_chr2, short_circ=annotate_short_circ)
            chr1in = len(chr1_segments_in) > 0
            chr2in = len(chr2_segments_in) > 0
            if include_both:
                retain = chr1in and chr2in
            else:
                retain = chr1in or chr2in
            if retain and annotate_retained:
                annotations_segments.extend(chr1_segments_in)
                annotations_segments.extend(chr2_segments_in)
            if include_spanning and chr1 == chr2 and len(include_segments_chr1) != 0:
                spanned_segments = coordinates_span_segments(coordinate1=adj.position1.coordinate, coordinate2=adj.position2.coordinate, segments=include_segments_chr1,
                                                             partial=False,
                                                             short_circ=annotate_short_circ)
                retain |= len(spanned_segments) > 0
                if retain and annotate_retained:
                    annotations_segments.extend(spanned_segments)
        elif len(include_by_chr) > 0:
            retain = False
        if not retain:
            continue
        exclude_segments_chr1 = exclude_by_chr.get(chr1, [])
        exclude_segments_chr2 = exclude_by_chr.get(chr2, [])
        if len(exclude_segments_chr1) != 0 or len(exclude_segments_chr2) != 0:
            chr1in = coordinate_in_any_segment(coordinate=adj.position1.coordinate, segments=exclude_segments_chr1)
            chr2in = coordinate_in_any_segment(coordinate=adj.position2.coordinate, segments=exclude_segments_chr2)
            if exclude_both:
                retain = not (chr1in and chr2in)
            else:
                retain = not (chr1in or chr2in)
            if exclude_spanning and chr1 == chr2 and len(exclude_segments_chr1) != 0:
                retain &= not coordinates_span_any_segment(coordinate1=adj.position1.coordinate, coordinate2=adj.position2.coordinate, segments=exclude_segments_chr1,
                                                           partial=False)
        if retain:
            if annotate_retained and len(annotations_segments) > 0:
                annotations = set()
                for segment in annotations_segments:
                    annotations.add(segment.extra.get(annotated_retained_segments_extra_field, segment.stable_id_non_hap))
                adj.extra[annotate_retained_extra_field] = sorted(annotations)
            yield adj


def filter_adjacencies_by_size(adjacencies, min_size=0, max_size=1000000000, size_extra_field=None, size_extra_seq_field=None,
                               allow_inter_chr=True, allow_intra_chr=True, size_extra_field_abs=True, allow_self_loops=True):
    for adj in adjacencies:
        adj_size = None
        try:
            adj_size = int(float(adj.extra[size_extra_field]))
            if size_extra_field_abs:
                adj_size = abs(adj_size)
        except (KeyError, ValueError):
            pass
        if adj_size is None:
            try:
                adj_size = len(adj.extra[size_extra_seq_field])
            except (KeyError, ValueError):
                pass
        if adj_size is None:
            adj_size = adj.distance_non_hap
        if adj.position1.chromosome != adj.position2.chromosome:
            if allow_inter_chr:
                yield adj
        elif allow_intra_chr:
            if not allow_self_loops and adj_size == 0:
                continue
            if adj_size < min_size or adj_size > max_size:
                continue
            yield adj


def coordinate_in_segments(coordinate, segments, short_circ=False):
    result = []
    segments = sorted(segments, key=lambda s: (s.start_coordinate, s.end_coordinate))
    for segment in segments:
        if coordinate < segment.start_coordinate:
            break
        if segment.start_coordinate <= coordinate <= segment.end_coordinate:
            result.append(segment)
            if short_circ:
                return result
    return result


def coordinates_span_segments(coordinate1, coordinate2, segments, partial=True, short_circ=False):
    result = []
    coordinate1, coordinate2 = sorted([coordinate1, coordinate2])
    segments = sorted(segments, key=lambda s: (s.start_coordinate, s.end_coordinate))
    for segment in segments:
        if coordinate2 < segment.start_coordinate:
            break
        if partial and (segment.start_coordinate <= coordinate1 <= segment.end_coordinate or segment.start_coordinate <= coordinate2 <= segment.end_coordinate):
            result.append(segment)
            if short_circ:
                return result
        elif coordinate1 <= segment.start_coordinate and coordinate2 >= segment.end_coordinate:
            result.append(segment)
            if short_circ:
                return result
    return result


def coordinate_in_any_segment(coordinate, segments):
    for segment in segments:
        if segment.start_coordinate <= coordinate <= segment.end_coordinate:
            return True
    return False


def coordinates_span_any_segment(coordinate1, coordinate2, segments, partial=True):
    for segment in segments:
        if partial and (segment.start_coordinate <= coordinate1 <= segment.end_coordinate or segment.start_coordinate <= coordinate2 <= segment.end_coordinate):
            return True
        elif coordinate1 <= segment.start_coordinate and coordinate2 >= segment.end_coordinate:
            return True
    return False


def get_shared_nas_parser():
    shared_parser = argparse.ArgumentParser(add_help=False)
    shared_parser.add_argument("--no-append-id", action="store_false", dest="append_id_suffix")
    shared_parser.add_argument("--o-extra-fields", default="all")
    shared_parser.add_argument("--chrs-include", action="append", nargs=1)
    shared_parser.add_argument("--chrs-include-file", type=argparse.FileType("rt"))
    shared_parser.add_argument("--chrs-include-no-both", action="store_false", dest="include_both")
    shared_parser.add_argument("--chrs-include-spanning", action="store_true", dest="include_spanning")
    shared_parser.add_argument("--chrs-exclude", action="append", nargs=1)
    shared_parser.add_argument("--chrs-exclude-file", type=argparse.FileType("rt"))
    shared_parser.add_argument("--chrs-exclude-both", action="store_true", dest="exclude_both")
    shared_parser.add_argument("--chrs-exclude-spanning", action="store_true", dest="exclude_spanning")
    return shared_parser


def get_chromosome_strip_parser():
    chr_strip_parser = argparse.ArgumentParser(add_help=False)
    chr_strip_parser.add_argument("--no-strip-chr", action="store_false", dest="strip_chr")
    return chr_strip_parser


def get_chromosome_regions_dict(string_entries):
    result = defaultdict(list)
    for entry in string_entries:
        if len(entry) == 0:
            continue
        data = entry.split(":")
        assert len(data) in [1, 2]
        chromosome = data[0]
        if len(data) > 1:
            start, stop = map(int, data[1].split("-"))
            result[chromosome].append((start, stop))
        else:
            result[chromosome].append((-1, 3000000000))
    for chromosome in result.keys():
        result[chromosome] = sorted(result[chromosome], key=lambda e: (e[0], e[1]))
    return result


def iter_over_string_entries_from_source(source):
    for line in source:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        yield line


def get_extra_field_regexes(string_entries):
    result = defaultdict(list)
    for entry in string_entries:
        data = entry.split("=")
        data = [data[0], "=".join(data[1:])]
        field_name = data[0]
        regex = re.compile(data[1])
        result[field_name].append(regex)
    return result


def position_in_region(position, region):
    return region[0] <= position.coordinate <= region[1]


def filter_adjacencies_by_extra(adjacencies,
                                keep_extra_field=None, keep_extra_field_missing_strategy=None,
                                remove_extra_field=None, remove_extra_field_missing_strategy=None):
    for adj in adjacencies:
        if keep_extra_field is not None and len(keep_extra_field) > 0:
            keep = False
            for field, regex_list in keep_extra_field.items():
                for regex in regex_list:
                    if field in adj.extra:
                        data = adj.extra[field]
                        if not isinstance(data, str):
                            if field == COPY_NUMBER:
                                data = stringify_adjacency_cn_entry(entry=data)
                            else:
                                data = str(data)
                        if regex.search(data) is not None:
                            keep = True
                    elif field.lower() == ADJACENCY_TYPE.lower():
                        data = str(adj.adjacency_type.value)
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
                    if field in adj.extra:
                        data = adj.extra[field]
                        if not isinstance(data, str):
                            if field == COPY_NUMBER:
                                data = stringify_adjacency_cn_entry(entry=data)
                            else:
                                data = str(data)
                        if regex.search(data) is not None:
                            keep = False
                            break
                    elif field.lower() == ADJACENCY_TYPE.lower():
                        data = str(adj.adjacency_type.value)
                        if regex.search(data) is not None:
                            keep = False
                            break
                    elif remove_extra_field_missing_strategy == REMOVE:
                        keep = False
                        break
            if not keep:
                continue
        yield adj


def processed_gundem2015_adjacencies(adjacencies, sample_names, min_per_sample_cnt, remove_per_sample_cnt_data=True):
    result = []
    for adj in adjacencies:
        for sample_name in sample_names:
            if GUNDEM_PER_SAMPLE_SUPPORT in adj.extra and \
                    sample_name in adj.extra[GUNDEM_PER_SAMPLE_SUPPORT] and \
                    adj.extra[GUNDEM_PER_SAMPLE_SUPPORT][sample_name] >= min_per_sample_cnt:
                result.append(adj)
                break
        if remove_per_sample_cnt_data and GUNDEM_PER_SAMPLE_SUPPORT in adj.extra:
            del adj.extra[GUNDEM_PER_SAMPLE_SUPPORT]
    return result


def refined_adjacencies_reciprocal(novel_adjacencies, max_distance, inplace=False):
    if not inplace:
        novel_adjacencies = deepcopy(novel_adjacencies)
    positions_by_chr = defaultdict(set)
    adjacencies_by_positions = defaultdict(list)
    for na in novel_adjacencies:
        p1_chr = na.position1.chromosome
        p2_chr = na.position2.chromosome
        positions_by_chr[p1_chr].add(copy(na.position1))
        positions_by_chr[p2_chr].add(copy(na.position2))
        adjacencies_by_positions[na.position1].append(na)
        adjacencies_by_positions[na.position2].append(na)
    new_positions_by_chr = {}
    for chr_name in sorted(positions_by_chr.keys()):
        new_positions_by_chr[chr_name] = sorted(positions_by_chr[chr_name], key=lambda p: p.coordinate)
    positions_by_chr = new_positions_by_chr
    for chr_name in sorted(positions_by_chr.keys()):
        positions = positions_by_chr[chr_name]
        if len(positions) < 2:
            continue
        suitable_for_merging = True
        for lp, rp in zip(positions[:-1], positions[1:]):
            if not suitable_for_merging:
                suitable_for_merging = True
                continue
            if positions_are_reciprocal(p1=lp, p2=rp, max_distance=max_distance):
                lp_adjacencies = adjacencies_by_positions[lp]
                rp_adjacencies = adjacencies_by_positions[rp]
                lp_adjacencies_ids = {a.stable_id_non_phased for a in lp_adjacencies}
                rp_adjacencies_ids = {a.stable_id_non_phased for a in rp_adjacencies}
                if len(lp_adjacencies_ids & rp_adjacencies_ids) > 0:
                    continue
                new_p1, new_p2 = merged_reciprocal_positions(p1=lp, p2=rp)
                if new_p1.strand == rp.strand:
                    new_p1, new_p2 = new_p2, new_p1
                for adj in lp_adjacencies:
                    update_position_in_adjacency(adjacency=adj, old_position=lp, new_position=new_p1)
                for adj in rp_adjacencies:
                    update_position_in_adjacency(adjacency=adj, old_position=rp, new_position=new_p2)
                suitable_for_merging = False
    return novel_adjacencies


def iter_haploid_adjacencies(adjacencies, copy=True):
    for adj in adjacencies:
        result = adj
        if copy:
            result = deepcopy(adj)
        if COPY_NUMBER in result.extra:
            cn_entry = result.extra[COPY_NUMBER]
            result_cn_entry = {}
            for clone_id, cn_dict in cn_entry.items():
                new_cn_dict = {Phasing.AA: sum(cn_dict.values())}
                result_cn_entry[clone_id] = new_cn_dict
            result.extra[COPY_NUMBER] = result_cn_entry
        yield result


def haploid_adjacencies(adjacencies, copy=True):
    return list(iter_haploid_adjacencies(adjacencies=adjacencies, copy=copy))


def get_haploid_acnt(adjacencies, acnt):
    result = {clone_id: AdjacencyCopyNumberProfile() for clone_id in sorted(acnt.keys())}
    for adj in adjacencies:
        aid = adj.stable_id_non_phased
        for clone_id in result.keys():
            result_acnp: AdjacencyCopyNumberProfile = result[clone_id]
            source_acnp: AdjacencyCopyNumberProfile = acnt[clone_id]
            total_cn = source_acnp.get_combined_cn(aid=aid)
            result_acnp.set_cn_record(aid=aid, phasing=Phasing.AA, cn=total_cn)
    return result


def positions_are_reciprocal(p1, p2, max_distance):
    if p1.chromosome != p2.chromosome:
        return False
    return abs(p1.coordinate - p2.coordinate) <= max_distance and p1.strand != p2.strand


def merged_reciprocal_positions(p1, p2):
    p1, p2 = deepcopy(p1), deepcopy(p2)
    fsp, rsp = (p1, p2) if p1.strand == Strand.FORWARD else (p2, p1)
    left = min([p1.coordinate, p2.coordinate])
    distance = abs(p1.coordinate - p2.coordinate)
    mid_coordinate = left + int(distance / 2)
    fsp.coordinate = mid_coordinate
    rsp.coordinate = mid_coordinate + 1
    return fsp, rsp


def update_position_in_adjacency(adjacency, old_position, new_position):
    if adjacency.position1 == old_position:
        adjacency.position1 = deepcopy(new_position)
    if adjacency.position2 == old_position:
        adjacency.position2 = deepcopy(new_position)


def element_before_window(window, element, adj_full_cnt=True):
    """
    Returning True in every ambiguous situation, as that simply skips the element
    """
    if isinstance(element, Position):
        return element.coordinate < window.start_coordinate
    elif isinstance(element, Adjacency):
        p1_same_chr = element.position1.chromosome == window.chromosome
        p2_same_chr = element.position2.chromosome == window.chromosome
        if adj_full_cnt:
            if not (p1_same_chr and p2_same_chr):
                return True
            if (p1_same_chr and element.position1.coordinate < window.start_coordinate) or (p2_same_chr and element.position2.coordinate < window.start_coordinate):
                return True
            return False
        else:
            if p1_same_chr and element.position1.coordinate < window.start_coordinate:
                return True
            elif p2_same_chr and element.position2.coordinate < window.start_coordinate:
                return True
            else:
                return False
    return True


def element_in_window(window, element, adj_full_cnt=True):
    if isinstance(element, Position):
        return window.start_coordinate <= element.coordinate <= window.end_coordinate
    elif isinstance(element, Adjacency):
        p1_same_chr = element.position1.chromosome == window.chromosome
        p2_same_chr = element.position2.chromosome == window.chromosome
        p1_in = p1_same_chr and window.start_coordinate <= element.position1.coordinate <= window.end_coordinate
        p2_in = p2_same_chr and window.start_coordinate <= element.position2.coordinate <= window.end_coordinate
        if adj_full_cnt:
            cnt = p1_in and p2_in
        else:
            cnt = p1_in or p2_in
        return cnt
    return False


def get_circa_adj_cnt(adjacencies, window_size=10000000, chr_sizes=None, element="breakend", adj_full_cnt=True):
    result = defaultdict(int)
    adjacencies_ids_by_chrs = defaultdict(set)
    adjacencies_by_ids = {adj.stable_id_non_phased: adj for adj in adjacencies}
    for adj in adjacencies_by_ids.values():
        adjacencies_ids_by_chrs[adj.position1.chromosome].add(adj.stable_id_non_phased)
        adjacencies_ids_by_chrs[adj.position2.chromosome].add(adj.stable_id_non_phased)
    adjacencies_by_chr = defaultdict(list)
    for chr_name in list(adjacencies_ids_by_chrs.keys()):
        sorted_chr_adjacencies = sorted([adjacencies_by_ids[aid] for aid in adjacencies_ids_by_chrs[chr_name]],
                                        key=lambda adj: (adj.position1.coordinate, adj.position2.coordinate))
        adjacencies_by_chr[chr_name] = sorted_chr_adjacencies
    windows_by_chr = defaultdict(list)
    if chr_sizes is None:
        chr_sizes = {}
    for chr_name in set(chr_sizes.keys()) | set(adjacencies_by_chr.keys()):
        start = 0
        default = adjacencies_by_chr[chr_name][-1].position2.coordinate if chr_name in adjacencies_by_chr else 0
        end = chr_sizes.get(chr_name, default)
        windows_boundaries = list(range(start, end, window_size))
        if windows_boundaries[-1] != end:
            windows_boundaries.append(end)
        for lb, rb in zip(windows_boundaries[:-1], windows_boundaries[1:]):
            segment = Segment.from_chromosome_coordinates(chromosome=chr_name, start=lb + 1, end=rb)
            windows_by_chr[chr_name].append(segment)
    # counted_entries = set()
    for chr_name in windows_by_chr.keys():
        chr_windows = iter(windows_by_chr[chr_name])
        current_window = next(chr_windows, None)
        elements = []
        if element == "breakend":
            for adj in adjacencies_by_chr[chr_name]:
                for p in [adj.position1, adj.position2]:
                    if p.chromosome == chr_name:
                        elements.append(p)
        elif element == "adj":
            for adj in adjacencies_by_chr[chr_name]:
                p1_in = adj.position1.chromosome == chr_name
                p2_in = adj.position2.chromosome == chr_name
                if adj_full_cnt:
                    cnt = p1_in and p2_in
                else:
                    cnt = p1_in or p2_in
                if cnt:
                    elements.append(adj)
        for entry in elements:
            if current_window is None:
                break
            # element_id = element.stable_id_non_hap if isinstance(element, Position) else element.stable_id_non_phased
            # if element_id in counted_entries:
            #     continue
            # counted_entries.add(element_id)
            if element_before_window(current_window, entry, adj_full_cnt=adj_full_cnt):
                continue
            elif element_in_window(current_window, entry, adj_full_cnt=adj_full_cnt):
                result[current_window] += 1
            else:  # after window
                current_window = next(chr_windows, None)
    return result


def update_adjacencies(target_adjacencies, source_adjacencies,
                       update_coords=True, update_coord1=True, update_coord2=True,
                       update_strands=True, update_strand1=True, update_strand2=True,
                       extra_exclude=None, extra_include=None, include_missing=True):
    if extra_exclude is None:
        extra_exclude = set()
    if extra_include is None:
        extra_include = set()
    source_adjacencies_by_aids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in source_adjacencies}
    for adj in target_adjacencies:
        source_adj: Adjacency = source_adjacencies_by_aids.get(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased), None)
        if source_adj is None:
            yield adj
        if update_coords:
            if update_coord1:
                adj.position1.coordinate = source_adj.position1.coordinate
            if update_coord2:
                adj.position2.coordinate = source_adj.position2.coordinate
        if update_strands:
            if update_strand1:
                adj.position1.strand = source_adj.position1.strand
            if update_strand2:
                adj.position2.strand = source_adj.position2.strand
        if len(extra_exclude) == 1 and list(extra_exclude)[0] == "all":
            yield adj
        updated_extra = deepcopy(adj.extra)
        for key in adj.extra.keys():
            if key in extra_exclude:
                continue
            if len(extra_include) > 0 and key not in extra_include:
                continue
            updated_extra[key] = source_adj.extra.get(key, updated_extra[key])
        if include_missing:
            for key, value in source_adj.extra.items():
                if key in updated_extra:
                    continue
                updated_extra[key] = value
        adj.extra = updated_extra
        yield adj
