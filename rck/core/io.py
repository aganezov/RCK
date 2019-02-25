# -*- coding: utf-8 -*-
import argparse
import ast
import csv
import sys
import datetime
import itertools
import logging
import os
from collections import defaultdict
from copy import deepcopy
from enum import Enum

from rck.core.structures import AdjacencyCopyNumberProfile, AdjacencyGroup, CNBoundaries, SegmentCopyNumberBoundaries, AdjacencyGroupType
from rck.core.structures import SegmentCopyNumberProfile, Haplotype, AdjacencyType, Phasing
from rck.core.structures import Position, Strand, Adjacency, Segment

csv.field_size_limit(sys.maxsize)

AID = "aid"
GID = "gid"
AIDS = "aids"
AG_TYPE = "agt"
AG_LABELING = "agl"
FALSE_POSITIVE = "fp"

CHR1 = "chr1"
COORD1 = "coord1"
CHR2 = "chr2"
COORD2 = "coord2"
STRAND1 = "strand1"
STRAND2 = "strand2"
EXTRA = "extra"
COPY_NUMBER = "cn"
OLD_SEGMENTS_CNS = "CNs"
OLD_ADJACENCIES_CNS = "CNs"
COPY_NUMBER_BOUNDARIES = "cnb"
ADJACENCY_TYPE = "at"

CHR = "chr"
START = "start"
END = "end"
COORD = "coord"
STRAND = "strand"

VCF_GENOTYPE = "GT"


def get_logging_cli_parser():
    logging_parser = argparse.ArgumentParser(add_help=False)
    logging_parser.add_argument("--log-level", type=int,
                                choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL], default=logging.INFO)
    logging_parser.add_argument("--log-format", default="%(asctime)s - %(name)-15s - %(levelname)-7s - %(message)s")
    logging_parser.add_argument("--log-file", default=None)
    return logging_parser


def get_standard_logger_from_args(args, program_name):
    logger = logging.getLogger(program_name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(args.log_format)
    ch = logging.StreamHandler()
    ch.setLevel(args.log_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    if args.log_file is not None:
        full_path = get_full_path(path=args.log_file)
        fh = logging.FileHandler(full_path)
        fh.setLevel(args.log_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger


def extra_description_is_all(extra):
    return (len(extra) == 1 and extra[0].lower() == "all") or (isinstance(extra, str) and extra.lower() == "all")


def parse_cn_entry(cn_string, cn_separator=";"):
    result = dict()
    clone_specific_entries = cn_string.split(cn_separator)
    for entry in clone_specific_entries:
        entry = entry[1:-1]
        clone_id, cna, cnb = entry.split(",")
        cna, cnb = int(cna), int(cnb)
        result[clone_id] = {Haplotype.A: cna, Haplotype.B: cnb}
    return result


def parse_scn_record(scn_string, separator="\t", cn_separator=";", allow_unit_segments=False):
    data = scn_string.split(separator)
    chr_name = data[0]
    start_coordinate = int(data[1])
    end_coordinate = int(data[2])
    cn_entry = data[3]
    start_position = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
    end_position = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    segment = Segment(start_position=start_position, end_position=end_position, allow_unit_length=allow_unit_segments)
    cns = parse_cn_entry(cn_string=cn_entry, cn_separator=cn_separator)
    return segment, cns


def read_scn_tensor(file_name, clone_ids, separator="\t", cn_separator=";", allow_absent_clones=False, allow_unit_segments=False):
    segments = []
    result = {}
    for clone_id in clone_ids:
        result[clone_id] = SegmentCopyNumberProfile()
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            segment, cns = parse_scn_record(scn_string=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in clone_ids:
                if not allow_absent_clones and clone_id not in cns:
                    raise Exception()
            segments.append(segment)
            sid = segment.stable_id_non_hap
            for clone_id in clone_ids:
                scnp = result[clone_id]
                cna = cns[clone_id][Haplotype.A]
                cnb = cns[clone_id][Haplotype.B]
                scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=cna)
                scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=cnb)
    return segments, result


def get_cn_string_entry(clone_id, cna, cnb):
    return "({clone_id},{cna},{cnb})".format(clone_id=str(clone_id), cna=str(cna), cnb=str(cnb))


def get_scn_string_entry(segment, scnt, clone_ids, separator="\t", cn_separator=";"):
    chr_string = "{chr_name}".format(chr_name=str(segment.chromosome))
    start_coord_string = "{start_coord}".format(start_coord=str(segment.start_position.coordinate))
    end_coord_string = "{end_coord}".format(end_coord=str(segment.end_position.coordinate))
    cn_strings = []
    sid = segment.stable_id_non_hap
    for clone_id in clone_ids:
        scnp = scnt[clone_id]
        cna = scnp.get_cn(sid=sid, haplotype=Haplotype.A, default=0)
        cnb = scnp.get_cn(sid=sid, haplotype=Haplotype.B, default=0)
        cn_string = get_cn_string_entry(clone_id=clone_id, cna=cna, cnb=cnb)
        cn_strings.append(cn_string)
    cn_string = cn_separator.join(cn_strings)
    all_string_entries = [chr_string, start_coord_string, end_coord_string, cn_string]
    return separator.join(all_string_entries)


def write_scn_tensor(file_name, scnt, segments, separator="\t", cn_separator=";"):
    clone_ids = sorted(scnt.keys())
    with open(file_name, "wt") as destination:
        print("chr", "start", "end", "CNs", sep=separator, file=destination)
        segments = sorted(segments, key=lambda s: (s.chromosome, s.start_position.coordinate, s.end_position.coordinate))
        for segment in segments:
            scn_string = get_scn_string_entry(segment=segment, scnt=scnt, clone_ids=clone_ids, separator=separator, cn_separator=cn_separator)
            print(scn_string, file=destination)


def read_scnt_from_file(file_name, clone_ids=None, separator="\t", extra_separator=";", remove_cn_data_from_segs=True):
    with open(file_name, "rt") as source:
        return read_scnt_from_source(source=source, clone_ids=clone_ids, separator=separator, extra_separator=extra_separator,
                                     remove_cn_data_from_segs=remove_cn_data_from_segs)


def read_scnt_from_source(source, clone_ids=None, separator="\t", extra_separator=";", remove_cn_data_from_segs=True):
    segments = read_segments_from_source(source=source, separator=separator, extra_separator=extra_separator)
    for segment in segments:
        if COPY_NUMBER not in segment.extra:
            raise ValueError("Trying to read Segment Copy Number tensor from {source}, but the {cn} entry in the {extra} column (or {extra} column is missing al together)"
                             "".format(source=source, cn=COPY_NUMBER, extra=EXTRA))
    scnt = extract_scnt_from_segments(segments=segments, clone_ids=clone_ids)
    if remove_cn_data_from_segs:
        remove_cn_data_from_segments(segments=segments)
    return segments, scnt


def remove_cn_data_from_segments(segments):
    for segment in segments:
        if COPY_NUMBER in segment.extra:
            del segment.extra[COPY_NUMBER]


def read_scnb_from_file(file_name, clone_ids=None, separator="\t", extra_separator=";", remove_cnb_data_from_segs=True, allow_missing=True):
    with open(file_name, "rt") as source:
        return read_scnb_from_source(source=source, clone_ids=clone_ids, separator=separator, extra_separator=extra_separator,
                                     remove_cnb_data_from_segs=remove_cnb_data_from_segs, allow_missing=allow_missing)


def read_scnb_from_source(source, clone_ids=None, separator="\t", extra_separator=";", remove_cnb_data_from_segs=True, allow_missing=True):
    segments = read_segments_from_source(source=source, separator=separator, extra_separator=extra_separator)
    if not allow_missing:
        for segment in segments:
            if COPY_NUMBER_BOUNDARIES not in segment:
                raise ValueError("Trying to read Segment Copy Number tensor from {source}, but the {cn} entry in the {extra} column (or {extra} column is missing al together)"
                                 "".format(source=source, cn=COPY_NUMBER, extra=EXTRA))
    scnb = extract_scnb_from_segments(segments=segments, clone_ids=clone_ids, allow_missing=allow_missing)
    if remove_cnb_data_from_segs:
        remove_cnb_data_from_segments(segments=segments)
    return segments, scnb


def remove_cnb_data_from_segments(segments):
    for segment in segments:
        if COPY_NUMBER_BOUNDARIES in segment.extra:
            del segment.extra[COPY_NUMBER_BOUNDARIES]


def extract_scnt_from_segments(segments, clone_ids=None):
    if clone_ids is None:
        clone_ids = set()
        for segment in segments:
            clone_ids.update(set(segment.extra[COPY_NUMBER].keys()))
    result = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    for segment in segments:
        cn_data = segment.extra[COPY_NUMBER]
        sid = segment.stable_id_non_hap
        for clone_id in clone_ids:
            if clone_id not in cn_data:
                raise ValueError("Clone {cid} is not present in CN data for segment {segment}".format(cid=clone_id, segment=sid))
            for haplotype in cn_data[clone_id]:
                scnp = result[clone_id]
                scnp.set_cn_record(sid=sid, hap=haplotype, cn=cn_data[clone_id][haplotype])
    return result


def extract_scnb_from_segments(segments, clone_ids=None, allow_missing=True):
    if clone_ids is None:
        clone_ids = set()
        for segment in segments:
            clone_ids.update(set(segment.extra[COPY_NUMBER].keys()))
    result = {clone_id: SegmentCopyNumberBoundaries() for clone_id in clone_ids}
    for segment in segments:
        if allow_missing and COPY_NUMBER_BOUNDARIES not in segment.extra:
            continue
        cnb_data = segment.extra[COPY_NUMBER_BOUNDARIES]
        sid = segment.stable_id_non_hap
        for clone_id in clone_ids:
            if clone_id not in cnb_data:
                raise ValueError("Clone {cid} is not present in CNB data for segment {segment}".format(cid=clone_id, segment=sid))
            for haplotype in cnb_data[clone_id]:
                for boundary_type in cnb_data[clone_id][haplotype]:
                    scnb = result[clone_id]
                    scnb.set_cnb_record(sid=sid, hap=haplotype, boundary_type=boundary_type, value=cnb_data[clone_id][haplotype][boundary_type])
    return result


def read_segments_from_file(file_name, separator="\t", extra_separator=";"):
    with open(file_name, "rt") as source:
        return read_segments_from_source(source=source, separator=separator, extra_separator=extra_separator)


def read_segments_from_source(source, separator="\t", extra_separator=";"):
    return list(stream_segments_from_source(source=source, separator=separator, extra_separator=extra_separator))


def parse_segment_extra_cn_boundaries_string(segment_cn_boundaries_string):
    clone_specific_dict = ast.literal_eval(node_or_string=segment_cn_boundaries_string)
    result = defaultdict(lambda: defaultdict(dict))
    for clone_id, data in clone_specific_dict.items():
        for haplotype_str in data:
            boundaries_dict = data[haplotype_str]
            haplotype = Haplotype[haplotype_str]
            for boundaries_string in boundaries_dict:
                boundary_value = int(boundaries_dict[boundaries_string])
                boundary = CNBoundaries.from_string(string=boundaries_string)
                result[clone_id][haplotype][boundary] = boundary_value
    return dict(result)


def parse_segment_extra_cn_string(segment_cn_string):
    clone_specific_dict = ast.literal_eval(node_or_string=segment_cn_string)
    result = defaultdict(dict)
    for clone_id, data in clone_specific_dict.items():
        for haplotype_str in data:
            cn = int(data[haplotype_str])
            haplotype = Haplotype[haplotype_str]
            result[clone_id][haplotype] = cn
    return dict(result)


def parse_segment_old_cn_string(old_cn_string):
    per_clone_entries = old_cn_string.split(";")
    result = {}
    for per_clone_entry in per_clone_entries:
        per_clone_entry = per_clone_entry[1:-1]
        clone, cna, cnb = per_clone_entry.split(",")
        cna, cnb = int(cna), int(cnb)
        result[clone] = {Haplotype.A: cna, Haplotype.B: cnb}
    return result


def parse_segment_extra(extra, extra_separator=";"):
    result = {}
    extra_entries = extra.split(extra_separator)
    for entry in extra_entries:
        data = entry.split("=")
        if len(data) == 1:
            if len(data[0]) > 0:
                result[data[0]] = True
            continue
        key, value = data
        if key == COPY_NUMBER:
            value = parse_segment_extra_cn_string(segment_cn_string=value)
        elif key == COPY_NUMBER_BOUNDARIES:
            value = parse_segment_extra_cn_boundaries_string(segment_cn_boundaries_string=value)
        result[key] = value
    return result


def stream_segments_from_source(source, separator="\t", extra_separator=";"):
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        chromosome = row[CHR]
        start = int(row[START])
        end = int(row[END])
        extra_dict = {}
        if OLD_SEGMENTS_CNS in row:
            extra_dict[COPY_NUMBER] = parse_segment_old_cn_string(old_cn_string=row[OLD_SEGMENTS_CNS])
        if EXTRA in row:
            supp_extra_dict = parse_segment_extra(extra=row[EXTRA], extra_separator=extra_separator)
            extra_dict.update(supp_extra_dict)
        segment = Segment.from_chromosome_coordinates(chromosome=chromosome, start=start, end=end)
        segment.extra = extra_dict
        yield segment


def write_segments_to_file(file_name, segments, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True):
    with open(file_name, "wt") as destination:
        write_segments_to_destination(destination=destination, segments=segments, separator=separator,
                                      extra=extra, extra_separator=extra_separator, extra_fill=extra_fill, sort_segments=sort_segments)


def write_segments_to_destination(destination, segments, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True):
    if sort_segments:
        segments = sorted(segments, key=lambda s: (s.chromosome, s.start_position.coordinate, s.end_position.coordinate))
    header_exntries = [CHR, START, END]
    if extra is not None:
        header_exntries.append(EXTRA)
    writer = csv.DictWriter(destination, fieldnames=header_exntries, delimiter=separator)
    writer.writeheader()
    for segment in segments:
        data = {}
        data[CHR] = segment.chromosome
        data[START] = segment.start_position.coordinate
        data[END] = segment.end_position.coordinate
        if extra is not None:
            extra_strings = []
            if extra_description_is_all(extra=extra):
                for key, value in segment.extra.items():
                    extra_strings.append("{extra_name}={extra_value}".format(extra_name=key.lower(), extra_value=value if value is not None else extra_fill))
            else:
                for entry in extra:
                    extra_strings.append("{extra_name}={extra_value}".format(extra_name=str(entry).lower(), extra_value=segment.extra.get(entry, extra_fill)))
            extra_string = extra_separator.join(extra_strings)
            data[EXTRA] = extra_string
        writer.writerow(data)


def write_scnt_to_file(file_name, segments, scnt, clone_ids=None, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True, inplace=True):
    with open(file_name, "wt") as destination:
        write_scnt_to_destination(destination=destination, scnt=scnt, segments=segments, clone_ids=clone_ids, separator=separator, extra=extra, extra_separator=extra_separator,
                                  extra_fill=extra_fill, sort_segments=sort_segments, inplace=inplace)


def iter_segments_scnt_dummy(segments, scnt, clone_ids=None, inplace=True):
    if clone_ids is None:
        clone_ids = sorted(scnt.keys())
    for segment in segments:
        if not inplace:
            segment = deepcopy(segment)
        if COPY_NUMBER in segment.extra:
            yield segment
            continue
        sid = segment.stable_id_non_hap
        segment.extra[COPY_NUMBER] = {}
        for clone_id in clone_ids:
            clone_specific_entries = {}
            for haplotype in [Haplotype.A, Haplotype.B]:
                clone_specific_entries[haplotype.value] = scnt[clone_id].get_cn(sid=sid, haplotype=haplotype, default=0)
            segment.extra[COPY_NUMBER][clone_id] = dict(clone_specific_entries)
        yield segment


def write_scnt_to_destination(destination, segments, scnt, clone_ids=None, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True, inplace=True):
    if clone_ids is None:
        clone_ids = sorted(scnt.keys())
    segments = iter_segments_scnt_dummy(segments=segments, scnt=scnt, clone_ids=clone_ids, inplace=inplace)
    if extra is None:
        extra = [COPY_NUMBER]
    elif not extra_description_is_all(extra) and COPY_NUMBER not in extra:
        extra.append(COPY_NUMBER)
    write_segments_to_destination(destination=destination, segments=segments, separator=separator, extra=extra, extra_separator=extra_separator, extra_fill=extra_fill,
                                  sort_segments=sort_segments)


def iter_segments_scnb_dummy(segments, scnb, clone_ids=None, inplace=True):
    if clone_ids is None:
        clone_ids = sorted(scnb.keys())
    for segment in segments:
        if not inplace:
            segment = deepcopy(segment)
        if COPY_NUMBER_BOUNDARIES in segment.extra:
            yield segment
            continue
        sid = segment.stable_id_non_hap
        segment.extra[COPY_NUMBER_BOUNDARIES] = {}
        for clone_id in clone_ids:
            clone_specific_entries = defaultdict(dict)
            for haplotype in [Haplotype.A, Haplotype.B]:
                for boundary_type in [CNBoundaries.LOWER, CNBoundaries.UPPER]:
                    clone_specific_entries[haplotype.value][boundary_type.value] = scnb[clone_id].get_cnb(sid=sid, hap=haplotype, boundary_type=boundary_type)
            segment.extra[COPY_NUMBER_BOUNDARIES][clone_id] = dict(clone_specific_entries)
        yield segment


def write_scnb_to_file(file_name, segments, scnb, clone_ids=None, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True, inplace=True):
    with open(file_name, "wt") as destination:
        write_scnb_to_destination(destination=destination, segments=segments, scnb=scnb, clone_ids=clone_ids, separator=separator, extra=extra, extra_fill=extra_fill,
                                  sort_segments=sort_segments, inplace=inplace)


def write_scnb_to_destination(destination, segments, scnb, clone_ids=None, separator="\t", extra="all", extra_separator=";", extra_fill="", sort_segments=True, inplace=True):
    if clone_ids is None:
        clone_ids = sorted(scnb.keys())
    segments = iter_segments_scnb_dummy(segments=segments, scnb=scnb, clone_ids=clone_ids, inplace=inplace)
    if extra is None:
        extra = [COPY_NUMBER_BOUNDARIES]
    elif not extra_description_is_all(extra) and COPY_NUMBER_BOUNDARIES not in extra:
        extra.append(COPY_NUMBER_BOUNDARIES)
    write_segments_to_destination(destination=destination, segments=segments, separator=separator, extra=extra, extra_separator=extra_separator, extra_fill=extra_fill,
                                  sort_segments=sort_segments)


EXTERNAL_NA_ID = AID
SVTYPE = "svtype"


def parse_adjacency_extra_cn_string(adjacency_cn_string):
    clone_specific_dict = ast.literal_eval(node_or_string=adjacency_cn_string)
    result = defaultdict(dict)
    for clone_id, data in clone_specific_dict.items():
        for phasing_str in data:
            cn = int(data[phasing_str])
            phasing = Phasing[phasing_str]
            result[clone_id][phasing] = cn
    return result


def parse_adjacency_old_cn_string(old_cn_string):
    clone_specific_data = old_cn_string.split(";")
    result = {}
    for clone_specific_string in clone_specific_data:
        clone_specific_string = clone_specific_string[1:-1]
        clone, cnaa, cnab, cnba, cnbb = clone_specific_string.split(",")
        cnaa, cnab, cnba, cnbb = int(cnaa), int(cnab), int(cnba), int(cnbb)
        result[clone] = {
            Phasing.AA: cnaa,
            Phasing.AB: cnab,
            Phasing.BA: cnba,
            Phasing.BB: cnbb
        }
    return result


def parse_adjacency_extra(extra, extra_separator=";"):
    result = {}
    extra_entries = extra.split(extra_separator)
    for entry in extra_entries:
        key, value = entry.split("=")
        if key == COPY_NUMBER:
            value = parse_adjacency_extra_cn_string(adjacency_cn_string=value)
        result[key] = value
    return result


def read_adjacencies_from_file(file_name, extra_separator=";", separator="\t"):
    with open(file_name, "rt") as source:
        return read_adjacencies_from_source(source=source, extra_separator=extra_separator, separator=separator)


def read_adjacencies_from_source(source, extra_separator=";", separator="\t"):
    return list(stream_adjacencies_from_source(source=source, extra_separator=extra_separator, separator=separator))


def stream_adjacencies_from_source(source, extra_separator=";", separator="\t"):
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        chr1 = row[CHR1]
        chr2 = row[CHR2]
        coord1 = int(row[COORD1])
        coord2 = int(row[COORD2])
        strand1 = Strand.from_pm_string(string=row[STRAND1])
        strand2 = Strand.from_pm_string(string=row[STRAND2])
        extra_dict = {EXTERNAL_NA_ID: row[EXTERNAL_NA_ID]}
        if OLD_ADJACENCIES_CNS in row:
            extra_dict[COPY_NUMBER] = parse_adjacency_old_cn_string(old_cn_string=row[OLD_ADJACENCIES_CNS])
        if EXTRA in row:
            supp_extra_dict = parse_adjacency_extra(extra=row[EXTRA], extra_separator=extra_separator)
            extra_dict.update(supp_extra_dict)
        pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
        pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
        adjacency = Adjacency(position1=pos1, position2=pos2, extra=extra_dict)
        if ADJACENCY_TYPE in extra_dict:
            adjacency.adjacency_type = AdjacencyType.from_name(name=extra_dict[ADJACENCY_TYPE])
            del extra_dict[ADJACENCY_TYPE]
        else:
            adjacency.adjacency_type = AdjacencyType.REFERENCE if extra_dict[EXTERNAL_NA_ID].startswith("r") else AdjacencyType.NOVEL
        yield adjacency


def write_adjacencies_to_file(file_name, adjacencies, extra="all", extra_fill="", extra_separator=";", separator="\t", sort_adjacencies=True):
    if sort_adjacencies:
        adjacencies = sorted(adjacencies, key=lambda a: (a.position1.chromosome, a.position1.coordinate, a.position2.chromosome, a.position2.coordinate))
    with open(file_name, "wt") as destination:
        write_adjacencies_to_destination(destination=destination, adjacencies=adjacencies,
                                         extra=extra, extra_fill=extra_fill, extra_separator=extra_separator, separator=separator,
                                         sort_adjacencies=False)


def write_adjacencies_to_destination(destination, adjacencies, extra="all", extra_fill="", extra_separator=";", separator="\t", sort_adjacencies=True):
    if sort_adjacencies:
        adjacencies = sorted(adjacencies, key=lambda a: (a.position1.chromosome, a.position1.coordinate, a.position2.chromosome, a.position2.coordinate))
    header_entries = [AID, CHR1, COORD1, STRAND1, CHR2, COORD2, STRAND2]
    if extra is not None:
        header_entries.append(EXTRA)
    writer = csv.DictWriter(destination, fieldnames=header_entries, delimiter=separator)
    writer.writeheader()
    for adjacency in adjacencies:
        data = {}
        data[AID] = adjacency.extra.get(EXTERNAL_NA_ID, adjacency.idx)
        data[CHR1] = str(adjacency.position1.chromosome)
        data[COORD1] = str(adjacency.position1.coordinate)
        data[STRAND1] = str(adjacency.position1.strand)
        data[CHR2] = str(adjacency.position2.chromosome)
        data[COORD2] = str(adjacency.position2.coordinate)
        data[STRAND2] = str(adjacency.position2.strand)
        if extra is not None:
            extra_strings = []
            if extra_description_is_all(extra=extra):
                for key, value in adjacency.extra.items():
                    if isinstance(value, (list, tuple)):
                        value = ",".join(map(str, value))
                    extra_strings.append("{extra_name}={extra_value}".format(extra_name=str(key).lower(), extra_value=value if value is not None else extra_fill))
                extra_strings.append("{extra_name}={extra_value}".format(extra_name=ADJACENCY_TYPE.lower(), extra_value=adjacency.adjacency_type.value))
            else:
                for entry in extra:
                    value = adjacency.extra.get(entry, extra_fill)
                    if isinstance(value, (list, tuple)):
                        value = ",".join(map(str, value))
                    extra_strings.append("{extra_name}={extra_value}".format(extra_name=str(entry).lower(), extra_value=value))
            extra_string = extra_separator.join(extra_strings)
            data[EXTRA] = extra_string
        writer.writerow(data)


def suitable_extra_entry(entry, extra="all"):
    if extra_description_is_all(extra=extra):
        return True
    return str(entry).lower() in extra


VCF_SAMPLE_EXTRA_FIELDS = {
    COPY_NUMBER.lower(),
    COPY_NUMBER_BOUNDARIES.lower(),
    VCF_GENOTYPE.lower()
}


def is_info_extra_entry(entry):
    return str(entry).lower() not in VCF_SAMPLE_EXTRA_FIELDS


def is_format_extra_entry(entry):
    return not is_info_extra_entry(entry=entry)


def get_vcf_info_string(adjacency, extra_fields, extra_fill=""):
    result = []
    result.append("CHR2={chromosome}".format(chromosome=adjacency.position2.chromosome))
    result.append("END={coordinate}".format(coordinate=adjacency.position2.coordinate))
    result.append("STRANDS={s1}{s2}".format(s1=str(adjacency.position1.strand), s2=str(adjacency.position2.strand)))
    extra_lower = {key.lower(): value for key, value in adjacency.extra.items() if key.upper() not in {"CHR2", "END", "STRANDS"}}
    for key in extra_fields:
        if key.upper() in {"CHR2", "END", "STRANDS"}:
            continue
        value = extra_lower.get(key.lower(), extra_fill)
        if isinstance(value, (list, tuple)):
            value = ",".join(map(str, value))
        result.append("{entry_id}={value}".format(entry_id=key, value=value))
    return ";".join(result)


def get_vcf_format_string(adjacency, clone_id, format_fields, extra_fill=""):
    result = []
    if COPY_NUMBER in adjacency.extra:
        total_cn = sum([adjacency.extra[COPY_NUMBER][clone_id][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]])
        genotype_value = "0/0" if total_cn == 0 else "./."
    else:
        genotype_value = "./."
    extra_lower = {key.lower(): value for key, value in adjacency.extra.items()}
    for entry_id in format_fields:
        if entry_id.lower() == VCF_GENOTYPE.lower():
            value = genotype_value
        elif entry_id.lower() == COPY_NUMBER.lower() and COPY_NUMBER in adjacency.extra:
            value = ",".join([str(adjacency.extra[COPY_NUMBER][clone_id][ph]) for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]])
        else:
            value = str(extra_lower.get(entry_id.lower(), extra_fill))
        result.append(value)
    return ":".join(result)


class VCFOutputFormat(Enum):
    SNIFFLES = "Sniffles"
    VCFProper = "VCF"


def write_adjacencies_to_vcf_sniffles_file(file_name, adjacencies, extra="all", extra_fill=".", sort_adjacencies=True):
    with open(file_name, "wt") as destination:
        write_adjacencies_to_vcf_sniffles_destination(destination=destination, adjacencies=adjacencies, extra=extra, extra_fill=extra_fill, sort_adjacencies=sort_adjacencies)


def write_adjacencies_to_vcf_sniffles_destination(destination, adjacencies, extra="all", extra_fill=".", sort_adjacencies=True,
                                                  dummy_clone="dummy_clone", clone_suffix=""):
    if sort_adjacencies:
        adjacencies = sorted(adjacencies, key=lambda a: (a.position1.chromosome, a.position1.coordinate, a.position2.chromosome, a.position2.coordinate))
    adjacencies = list(adjacencies)
    today = datetime.datetime.today()
    print("##fileformat=VCFv4.1", file=destination)
    print("##fileDate={year}{month:0>2}{day:0>2}".format(year=today.year, month=today.month, day=today.day), file=destination)
    print("##source=RCK", file=destination)
    chromosomes = set()
    svtypes = set()
    extra_fields = set()
    for adj in adjacencies:
        chromosomes.add(adj.position1.chromosome)
        chromosomes.add(adj.position2.chromosome)
        if SVTYPE in adj.extra:
            svtypes.add(adj.extra[SVTYPE])
        for key in adj.extra.keys():
            if suitable_extra_entry(entry=key, extra=extra):
                extra_fields.add(key)
    chromosomes = sorted(chromosomes)
    svtypes = sorted(svtypes)
    for chromosome in chromosomes:
        print("##contig=<ID={chromosome}>".format(chromosome=chromosome), file=destination)
    for svtype in svtypes:
        print("##ALT=<ID={svtype},Description=\"{descr}\">".format(svtype=svtype.upper(), descr=svtype), file=destination)
    print("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"chromosome for second position of the reported adjacency\">", file=destination)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"coordinate for second position of the reported adjacency\">", file=destination)
    print("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"strands _x__y_ where _x_ corresponds tot he strand for first position, and _y_ for the \">", file=destination)
    clone_ids = None
    extra_info_fields_and_numbers = {}
    extra_format_fields_and_numbers = {"GT": "1"}
    for adjacency in adjacencies:
        for key, value in adjacency.extra.items():
            if key not in extra_fields:
                continue
            if is_info_extra_entry(entry=key):
                entry_id = str(key).upper()
                entry_number = "." if isinstance(value, (list, tuple)) else "1"
                if entry_id in extra_info_fields_and_numbers and extra_info_fields_and_numbers[entry_id] != entry_number:
                    entry_number = "1"
                extra_info_fields_and_numbers[entry_id] = entry_number
            elif is_format_extra_entry(entry=key):
                entry_id = str(key).upper()
                entry_number = "1"
                extra_format_fields_and_numbers[entry_id] = entry_number
                if isinstance(value, dict):
                    if clone_ids is None:
                        clone_ids = set(value.keys())
                    else:
                        clone_ids &= set(value.keys())
    if clone_ids is None:
        clone_ids = [dummy_clone]
    clone_names = [clone_id + clone_suffix for clone_id in clone_ids]
    clone_ids = sorted(clone_ids)

    for extra_entry_id, extra_entry_number in extra_info_fields_and_numbers.items():
        if extra_entry_id in ["CHR2", "END", "STRANDS"]:
            continue
        print("##INFO=<ID={entry_id},Number={entry_number},Type=String,Description=\"\">".format(entry_id=extra_entry_id, entry_number=extra_entry_number), file=destination)
    for extra_entry_id, extra_entry_number in extra_format_fields_and_numbers.items():
        print("##FORMAT=<ID={entry_id},Number={entry_number},Type=String,Description=\"\">".format(entry_id=extra_entry_id, entry_number=extra_entry_number), file=destination)
    for clone_id in clone_names:
        print("##SAMPLE=<ID={sample_id}>".format(sample_id=clone_id), file=destination)

    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for clone_id in clone_names:
        header.append(clone_id)
    print("\t".join(header), file=destination)

    for adjacency in adjacencies:
        print(adjacency.position1.chromosome, end="\t", file=destination)  # CHROM
        print(adjacency.position1.coordinate, end="\t", file=destination)  # POS
        print(adjacency.extra.get(EXTERNAL_NA_ID, adjacency.idx), end="\t", file=destination)  # ID
        print("N", end="\t", file=destination)  # REF
        print("<{svtype}>".format(svtype=adjacency.extra.get(SVTYPE, "BND")), end="\t", file=destination)  # ALT
        print(".", end="\t", file=destination)  # QUAL
        print("PASS", end="\t", file=destination)  # FILTER
        print(get_vcf_info_string(adjacency=adjacency, extra_fields=sorted(extra_info_fields_and_numbers.keys()), extra_fill=extra_fill), end="\t", file=destination)  # INFO
        print(":".join(extra_format_fields_and_numbers.keys()), end="\t", file=destination)  # FORMAT
        for clone_id in clone_ids:
            print(get_vcf_format_string(adjacency=adjacency, clone_id=clone_id, format_fields=sorted(extra_format_fields_and_numbers.keys())), end="\t", file=destination)  # SAMPLE
        print(file=destination)


def write_acnt_to_file(file_name, acnt, adjacencies, clone_ids=None,
                       extra="all", extra_fill="", extra_separator=";", separator="\t", sort_adjacencies=True, output_reference=True, inplace=True,
                       mix_reference_and_novel=False):
    with open(file_name, "wt") as destination:
        write_acnt_to_destination(destination=destination, acnt=acnt, adjacencies=adjacencies, clone_ids=clone_ids, extra=extra, extra_fill=extra_fill,
                                  extra_separator=extra_separator, separator=separator, sort_adjacencies=sort_adjacencies, output_reference=output_reference, inplace=inplace,
                                  mix_reference_and_novel=mix_reference_and_novel)


def iter_adjacencies_acnt_dummy(adjacencies, acnt, clone_ids=None, inplace=True):
    if clone_ids is None:
        clone_ids = sorted(acnt.keys())
    for adjacency in adjacencies:
        if not inplace:
            adjacency = deepcopy(adjacency)
        if COPY_NUMBER in adjacency.extra:
            yield adjacency
            continue
        aid = adjacency.stable_id_non_phased
        adjacency.extra[COPY_NUMBER] = {}
        for clone_id in clone_ids:
            clone_specific_entries = {}
            for phasing in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                clone_specific_entries[str(phasing)] = acnt[clone_id].get_cn(aid=aid, phasing=phasing, default=0)
            adjacency.extra[COPY_NUMBER][clone_id] = clone_specific_entries
        yield adjacency


def write_acnt_to_destination(destination, acnt, adjacencies, clone_ids=None,
                              extra="all", extra_fill="", extra_separator=";", separator="\t", sort_adjacencies=True, output_reference=True, inplace=False,
                              mix_reference_and_novel=False):
    entries = list(adjacencies)
    if sort_adjacencies:
        entries = sorted(entries, key=lambda a: (a.position1.chromosome, a.position1.coordinate, a.position2.chromosome, a.position2.coordinate))
    if not mix_reference_and_novel:
        ref_adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, adjacencies))
        nov_adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, adjacencies))
        entries = nov_adjacencies + ref_adjacencies
    if not output_reference:
        entries = [adj for adj in entries if adj.adjacency_type == AdjacencyType.NOVEL]
    if clone_ids is None:
        clone_ids = sorted(acnt.keys())
    entries = iter_adjacencies_acnt_dummy(adjacencies=entries, acnt=acnt, clone_ids=clone_ids, inplace=inplace)
    if extra is None:
        extra = [COPY_NUMBER]
    elif not extra_description_is_all(extra=extra) and COPY_NUMBER not in extra:
        extra.append(COPY_NUMBER)
    write_adjacencies_to_destination(destination=destination, adjacencies=entries,
                                     extra=extra, extra_fill=extra_fill, extra_separator=extra_separator, separator=separator, sort_adjacencies=False)


def read_acnt_from_source(source, clone_ids=None, extra_separator=";", separator="\t", remove_cn_data_from_adj=True):
    adjacencies = read_adjacencies_from_source(source=source, extra_separator=extra_separator, separator=separator)
    for adj in adjacencies:
        if COPY_NUMBER not in adj.extra:
            raise ValueError("Trying to read Adjacency Copy Number tensor from {source}, but the {cn} entry in the {extra} column (or {extra} column is missing al together)"
                             "".format(source=source, cn=COPY_NUMBER, extra=EXTRA))
    acnt = extract_acnt_from_adjacencies(adjacencies=adjacencies, clone_ids=clone_ids)
    if remove_cn_data_from_adj:
        remove_cn_data_from_adjacencies(adjacencies=adjacencies)
    return adjacencies, acnt


def remove_cn_data_from_adjacencies(adjacencies):
    for adj in adjacencies:
        if COPY_NUMBER in adj.extra:
            del adj.extra[COPY_NUMBER]


def extract_acnt_from_adjacencies(adjacencies, clone_ids=None):
    if clone_ids is None:
        clone_ids = set()
        for adj in adjacencies:
            clone_ids.update(set(adj.extra[COPY_NUMBER].keys()))
    clone_ids = sorted(clone_ids)
    acnt = {clone_id: AdjacencyCopyNumberProfile() for clone_id in clone_ids}
    for adj in adjacencies:
        cn_data = adj.extra[COPY_NUMBER]
        aid = adj.stable_id_non_phased
        for clone_id in clone_ids:
            if clone_id not in cn_data:
                raise ValueError("Clone {cid} is not present in CN data for adjacency {adj}".format(cid=clone_id, adj=aid))
            for phasing in cn_data[clone_id]:
                acnp = acnt[clone_id]
                acnp._set_cn_record(aid=aid, phasing=phasing, cn=cn_data[clone_id][phasing])
    return acnt


def read_acnt_from_file(file_name, clone_ids=None, extra_separator=";", separator="\t", remove_cn_data_from_adj=True):
    with open(file_name, "rt") as source:
        return read_acnt_from_source(source=source, clone_ids=clone_ids, extra_separator=extra_separator, separator=separator, remove_cn_data_from_adj=remove_cn_data_from_adj)


def parse_acn_string(string_entry, separator="\t", cn_separator=";"):
    cns = {}
    data = string_entry.split(separator)
    naid = data[0]
    chr1 = data[1]
    coord1 = int(data[2])
    strand1 = Strand.from_pm_string(string=data[3])
    chr2 = data[4]
    coord2 = int(data[5])
    strand2 = Strand.from_pm_string(string=data[6])
    cns_string_entry = data[7]
    p1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
    p2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
    adj_type = AdjacencyType.REFERENCE if "r" in naid else AdjacencyType.NOVEL
    adjacency = Adjacency(position1=p1, position2=p2, adjacency_type=adj_type, extra={EXTERNAL_NA_ID: naid})
    cn_data = cns_string_entry.split(cn_separator)
    for entry in cn_data:
        entry = entry[1:-1]
        clone_id, cnAA, cnAB, cnBA, cnBB = entry.split(",")
        cnAA, cnAB, cnBA, cnBB = int(cnAA), int(cnAB), int(cnBA), int(cnBB)
        cns[clone_id] = {Phasing.AA: cnAA, Phasing.AB: cnAB, Phasing.BA: cnBA, Phasing.BB: cnBB}
    return adjacency, cns


def read_acnt_old(file_name, clone_ids, separator="\t", cn_separator=";"):
    clone_ids = sorted(clone_ids)
    result = {clone_id: AdjacencyCopyNumberProfile() for clone_id in clone_ids}
    adjacencies = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            adjacency, cns = parse_acn_string(string_entry=line, separator=separator, cn_separator=cn_separator)
            if EXTERNAL_NA_ID in adjacency.extra and adjacency.extra[EXTERNAL_NA_ID].startswith("r"):
                adjacency.adjacency_type = AdjacencyType.REFERENCE
            else:
                adjacency.adjacency_type = AdjacencyType.NOVEL
            aid = adjacency.stable_id_non_phased
            adjacencies.append(adjacency)
            for clone_id in clone_ids:
                if clone_id not in cns:
                    raise Exception()
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    acnp = result[clone_id]
                    acnp._set_cn_record(aid=aid, phasing=ph, cn=cns[clone_id][ph])
    return adjacencies, result


def get_all_clone_ids_from_dasck_scnt_file(file_name, separator="\t", cn_separator=";", allow_unit_segments=False):
    clone_ids = set()
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            segment, cns = parse_scn_record(scn_string=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in cns.keys():
                clone_ids.add(clone_id)
    return sorted(clone_ids)


def parser_position_string_entry(position_string, separator):
    data = position_string.split(separator)
    chr_name = data[0]
    coordinate = int(data[1])
    strand = Strand.from_pm_string(string=data[2])
    return Position(chromosome=chr_name, coordinate=coordinate, strand=strand)


def read_positions(file_name, separator="\t"):
    result = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            position = parser_position_string_entry(position_string=line, separator=separator)
            result.append(position)
    return result


def read_positions_from_file(file_name, separator="\t", extra_separator=";"):
    with open(file_name, "rt") as source:
        return read_positions_from_source(source=source, separator=separator, extra_separator=extra_separator)


def read_positions_from_source(source, separator="\t", extra_separator=";"):
    return list(stream_positions_from_source(source=source, separator=separator, extra_separator=extra_separator))


def stream_positions_from_source(source, separator="\t", extra_separator=";"):
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        chromosome = row[CHR]
        coordinate = int(row[COORD])
        strand = Strand.from_pm_string(string=row[STRAND])
        yield Position(chromosome=chromosome, coordinate=coordinate, strand=strand)


def read_adjacency_groups_from_file(file_name, default_group_type=AdjacencyGroupType.MOLECULE, separator="\t", aids_separator=",", extra_separator=";"):
    with open(file_name, "rt") as source:
        return read_adjacency_groups_from_source(source=source, default_group_type=default_group_type, separator=separator,
                                                 aids_separator=aids_separator, extra_separator=extra_separator)


def read_adjacency_groups_from_source(source, default_group_type=AdjacencyGroupType.MOLECULE, separator="\t", aids_separator=",", extra_separator=";"):
    return list(stream_adjacency_groups_from_source(source=source, default_group_type=default_group_type, separator=separator,
                                                    aids_separator=aids_separator, extra_separator=extra_separator))


def stream_adjacency_groups_from_source(source, default_group_type=AdjacencyGroupType.MOLECULE, separator="\t", aids_separator=",", extra_separator=";"):
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        gid = row[GID]
        aids = row[AIDS].split(aids_separator)
        extra_dict = {}
        if EXTRA in row:
            extra_dict = parse_adj_groups_extra(extra_string=row[EXTRA], extra_separator=extra_separator)
        group_type = extra_dict.get(AG_TYPE, default_group_type)
        yield AdjacencyGroup(gid=gid, aids=aids, group_type=group_type, extra=extra_dict)


def write_positions_to_file(file_name, positions, extra="all", extra_fill="", separator="\t", extra_separator=";", sort=True):
    with open(file_name, "wt") as destination:
        write_positions_to_destinations(destination=destination, positions=positions, extra=extra, extra_fill=extra_fill, separator=separator, extra_separator=extra_separator,
                                        sort=sort)


def write_positions_to_destinations(destination, positions, extra="all", extra_fill="", separator="\t", extra_separator=";", sort=True):
    header_entries = [CHR, COORD, STRAND]
    if sort:
        positions = sorted(positions, key=lambda p: (p.chromosome, p.coordinate, p.strand))
    writer = csv.DictWriter(destination, fieldnames=header_entries, delimiter=separator)
    writer.writeheader()
    for position in positions:
        data = {}
        data[CHR] = position.chromosome
        data[COORD] = position.coordinate
        data[STRAND] = position.strand
        writer.writerow(data)


def parse_adj_group_extra_labeling_string(extra_labeling_string):
    return list(map(int, extra_labeling_string.split(",")))


def parse_adj_groups_extra(extra_string, extra_separator=";"):
    data = extra_string.split(extra_separator)
    result = {}
    for entry in data:
        key, value = entry.split("=")
        key = key.lower()
        if key == AG_TYPE:
            value = AdjacencyGroupType.from_string(string=value)
        elif key == FALSE_POSITIVE:
            value = float(value)
        elif key == AG_LABELING:
            value = list(map(int, value.split(",")))
        result[key] = value
    return result


def write_adjacency_groups_to_file(file_name, adjacency_groups, separator="\t", aids_separator=",",
                                   extra="all", extra_separator=";", extra_fill=""):
    with open(file_name, "wt") as destination:
        write_adjacency_groups_to_destination(destination=destination, adjacency_groups=adjacency_groups, separator=separator,
                                              aids_separator=aids_separator, extra=extra, extra_separator=extra_separator, extra_fill=extra_fill)


def write_adjacency_groups_to_destination(destination, adjacency_groups, separator="\t", aids_separator=",",
                                          extra="all", extra_separator=";", extra_fill=""):
    header_entries = [GID, AIDS, EXTRA]
    if extra is None:
        extra = [AG_TYPE]
    # if extra is not None:
    #     header_entries.append(EXTRA)
    writer = csv.DictWriter(destination, fieldnames=header_entries, delimiter=separator)
    writer.writeheader()
    for ag in adjacency_groups:
        data = {}
        data[GID] = ag.gid
        data[AIDS] = aids_separator.join(map(str, ag.adjacencies_ids))
        extra_strings = ["{ag_type_string}={ag_type}".format(ag_type_string=AG_TYPE, ag_type=ag.group_type.value)]
        if extra_description_is_all(extra=extra):
            for key, value in ag.extra.items():
                if isinstance(value, (list, tuple)):
                    value = ",".join(map(str, value))
                elif isinstance(value, AdjacencyGroupType):
                    value = value.value
                extra_strings.append("{extra_name}={extra_value}".format(extra_name=str(key).lower(), extra_value=value if value is not None else extra_fill))
        else:
            for entry in extra:
                if entry == AG_TYPE:
                    continue
                value = ag.extra.get(entry, extra_fill)
                if isinstance(value, (list, tuple)):
                    value = ",".join(map(str, value))
                extra_strings.append("{extra_name}={extra_value}".format(extra_name=str(entry).lower(), extra_value=value))
        agt_present = False
        for extra_string in extra_strings:
            agt_present |= "{ag_type_string}=".format(ag_type_string=AG_TYPE) in extra_string
        if not agt_present:
            extra_strings.append("{ag_type_string}={ag_type}".format(ag_type_string=AG_TYPE, ag_type=ag.group_type.value))
        extra_string = extra_separator.join(extra_strings)
        data[EXTRA] = extra_string
        writer.writerow(data)


def get_segment_string_entry(segment, separator="\t"):
    chr_name_entry = "{chr_name}".format(chr_name=str(segment.start_position.chromosome))
    start_entry = "{start_coord}".format(start_coord=str(segment.start_position.coordinate))
    end_entry = "{end_coord}".format(end_coord=str(segment.end_position.coordinate))
    return separator.join([chr_name_entry, start_entry, end_entry])


def write_segments(file_name, segments, separator="\t"):
    with open(file_name, "wt") as destination:
        segments = sorted(segments, key=lambda s: (s.start_position.chromosome, s.start_position.coordinate, s.end_position.coordinate))
        print("chr", "start", "end", sep=separator, file=destination)
        for segment in segments:
            string_entry = get_segment_string_entry(segment=segment, separator=separator)
            print(string_entry, file=destination)


def parse_segment_string_entry(string_entry, separator="\t", allow_unit_segments=False):
    data = string_entry.split(separator)
    chr_name = data[0]
    start_corodinate = int(data[1])
    end_coordinate = int(data[2])
    sp = Position(chromosome=chr_name, coordinate=start_corodinate, strand=Strand.REVERSE)
    ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    return Segment(start_position=sp, end_position=ep, allow_unit_length=allow_unit_segments)


def read_segments(file_name, separator="\t", allow_unit_segments=False):
    result = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            segment = parse_segment_string_entry(string_entry=line, separator=separator, allow_unit_segments=allow_unit_segments)
            result.append(segment)
    return result


def get_cn_boundaries_string_entry(segment, scn_boundaries, clone_ids, cn_separator=";"):
    result = []
    sid = segment.stable_id_non_hap
    for clone_id in clone_ids:
        lower_a, upper_a = scn_boundaries[clone_id][sid][Haplotype.A]
        lower_b, upper_b = scn_boundaries[clone_id][sid][Haplotype.B]
        string_a = "{lower}-{upper}".format(lower=lower_a, upper=upper_a)
        string_b = "{lower}-{upper}".format(lower=lower_b, upper=upper_b)
        string_entry = "({clone_id},{boundaries_a},{boundaries_b})".format(clone_id=str(clone_id), boundaries_a=str(string_a), boundaries_b=str(string_b))
        result.append(string_entry)
    return cn_separator.join(result)


def write_scn_boundaries(file_name, segments, scn_boundaries, separator="\t", cn_separator=";"):
    clone_ids = sorted(scn_boundaries.keys())
    with open(file_name, "wt") as destination:
        print("chr", "start", "end", "CNBoundaries", sep=separator, file=destination)
        for s in segments:
            string_entries = []
            string_entries.append("{chr_name}".format(chr_name=str(s.start_position.chromosome)))
            string_entries.append("{start}".format(start=str(s.start_position.coordinate)))
            string_entries.append("{end}".format(end=str(s.end_position.coordinate)))
            string_entries.append(get_cn_boundaries_string_entry(segment=s, scn_boundaries=scn_boundaries, clone_ids=clone_ids, cn_separator=cn_separator))
            string_entry = separator.join(string_entries)
            print(string_entry, file=destination)


def parse_cn_boundaries_entry_string(string_entry):
    data = string_entry.split("-")
    return int(data[0]), int(data[1])


def parse_cn_boundaries_string(string_entry, cn_separator=";"):
    result = {}
    data = string_entry.split(cn_separator)
    for entry in data:
        entry = entry[1:-1]
        entry_data = entry.split(",")
        clone_id = entry_data[0]
        boundaries_a_string = entry_data[1]
        boundaries_b_string = entry_data[2]
        boundaries_a = parse_cn_boundaries_entry_string(string_entry=boundaries_a_string)
        boundaries_b = parse_cn_boundaries_entry_string(string_entry=boundaries_b_string)
        result[clone_id] = {Haplotype.A: boundaries_a, Haplotype.B: boundaries_b}
    return result


def parse_scn_boundaries_string_entry(string_entry, separator="\t", cn_separator=";", allow_unit_segments=False):
    data = string_entry.split(separator)
    chr_name = data[0]
    start_coordinate = int(data[1])
    end_coordinate = int(data[2])
    cn_boundaries_string = data[3]
    sp = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
    ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    segment = Segment(start_position=sp, end_position=ep, allow_unit_length=allow_unit_segments)
    boundaries = parse_cn_boundaries_string(string_entry=cn_boundaries_string, cn_separator=cn_separator)
    return segment, boundaries


def read_scn_boundaries(file_name, segments, clone_ids, separator="\t", cn_separator=";", allow_unit_segments=False):
    result = {clone_id: {s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in segments} for clone_id in clone_ids}
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            segment, boundaries = parse_scn_boundaries_string_entry(string_entry=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in clone_ids:
                result[clone_id][segment.stable_id_non_hap][Haplotype.A] = boundaries[clone_id][Haplotype.A]
                result[clone_id][segment.stable_id_non_hap][Haplotype.B] = boundaries[clone_id][Haplotype.B]
    return result


def get_full_path(path):
    return os.path.abspath(os.path.expanduser(path))
