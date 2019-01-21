import csv
import math

from rck.core.structures import SegmentCopyNumberProfile, Segment, Haplotype
from rck.utils.adj.convert import strip_chr

BATTENBERG_SAMPLE_NAME = "sample"
BATTENBERG_CHROMOSOME = "chr"
BATTENBERG_START_POSITION = "startpos"
BATTENBERG_END_POSITION = "endpos"
BATTENBERG_CLONE1_CN_A = "nMaj1_A"
BATTENBERG_CLONE1_CN_B = "nMin1_A"
BATTENBERG_CLONE2_CN_A = "nMaj2_A"
BATTENBERG_CLONE2_CN_B = "nMin2_A"


def battenberg_force_non_negativity(cn):
    return cn if cn >= 0 else 0


def battenberg_get_subclonal_cn(subclonal_cn_string, clonal_cn_int):
    if subclonal_cn_string == "NA":
        return clonal_cn_int
    return battenberg_force_non_negativity(int(subclonal_cn_string))


def battenberg_get_scnt_from_battenberg_file(file_name, sample_name, separator="\t", chr_strip=True):
    with open(file_name, "rt") as source:
        return get_scnt_from_battenberg_source(source=source, sample_name=sample_name, separator=separator, chr_strip=chr_strip)


def get_scnt_from_battenberg_source(source, sample_name, separator="\t", chr_strip=True):
    clone1_name = "1"
    clone2_name = "2"
    scnt = {clone1_name: SegmentCopyNumberProfile(), clone2_name: SegmentCopyNumberProfile()}
    segments = []
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        if BATTENBERG_SAMPLE_NAME in row and row[BATTENBERG_SAMPLE_NAME] != sample_name:
            continue
        start_coordinate = int(row[BATTENBERG_START_POSITION])
        end_coordinate = int(row[BATTENBERG_END_POSITION])
        chromosome = row[BATTENBERG_CHROMOSOME]
        if chr_strip:
            chromosome = strip_chr(chr_string=chromosome)
        segment = Segment.from_chromosome_coordinates(chromosome=chromosome, start=start_coordinate, end=end_coordinate)
        clone1_scnp = scnt[clone1_name]
        clone2_scnp = scnt[clone2_name]
        cn1a = battenberg_force_non_negativity(int(row[BATTENBERG_CLONE1_CN_A]))
        cn1b = battenberg_force_non_negativity(int(row[BATTENBERG_CLONE1_CN_B]))
        clone1_scnp.set_cn_record_for_segment(segment=segment, cn=cn1a, haplotype=Haplotype.A)
        clone1_scnp.set_cn_record_for_segment(segment=segment, cn=cn1b, haplotype=Haplotype.B)
        cn2a = battenberg_get_subclonal_cn(subclonal_cn_string=row[BATTENBERG_CLONE2_CN_A], clonal_cn_int=cn1a)
        cn2b = battenberg_get_subclonal_cn(subclonal_cn_string=row[BATTENBERG_CLONE2_CN_B], clonal_cn_int=cn1b)
        clone2_scnp.set_cn_record_for_segment(segment=segment, cn=cn2a, haplotype=Haplotype.A)
        clone2_scnp.set_cn_record_for_segment(segment=segment, cn=cn2b, haplotype=Haplotype.B)
        segments.append(segment)
    return segments, scnt


def hatchet_get_clone_ids_from_file(file_name, sample_name, separator="\t", min_usage=0.01):
    result = set()
    candidates = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            data = line.split(separator)
            clone_data = data[6:]
            if line_cnt == 0:
                total_clone_cnt = int(len(clone_data) / 2)
                candidates = [str(cnt) for cnt in range(1, total_clone_cnt + 1)]
            if line.startswith("#"):
                continue
            sample = data[3]
            if sample != sample_name:
                continue
            for candidate_clone_id, clone_usage_str in zip(candidates, clone_data[1::2]):
                clone_usage = float(clone_usage_str)
                if clone_usage < min_usage:
                    continue
                result.add(candidate_clone_id)
            if sorted(result) == candidates:
                return sorted(result)
    return sorted(result)


def get_scnt_from_hatchet_file(file_name, sample_name, separator="\t", clone_ids=None, min_usage=0.01, chr_strip=True):
    if clone_ids is None:
        clone_ids = hatchet_get_clone_ids_from_file(file_name=file_name, sample_name=sample_name, separator=separator, min_usage=min_usage)
    with open(file_name, "rt") as source:
        return get_scnt_from_hatchet_source(source=source, separator=separator, clone_ids=clone_ids, chr_strip=chr_strip)


def get_scnt_from_hatchet_source(source, clone_ids, separator="\t", chr_strip=True):
    scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    segments = []
    clone_id_mappings = {}
    for line_cnt, line in enumerate(source):
        line = line.strip()
        data = line.split(separator)
        clone_data = data[6:]
        if line_cnt == 0:
            total_clone_cnt = int(len(clone_data) / 2)
            candidates = [str(cnt) for cnt in range(1, total_clone_cnt + 1)]
            for position_cnt, candidate in enumerate(candidates):
                if candidate in clone_ids:
                    clone_id_mappings[candidate] = position_cnt
        clone_cn_strs = clone_data[::2]
        if line.startswith("#") or len(line) == 0:
            continue
        sample_name = data[3]
        if sample_name != sample_name:
            continue
        chromosome = data[0]
        if chr_strip:
            chromosome = strip_chr(chr_string=chromosome)
        start_coord = int(data[1])
        end_coord = int(data[2]) - 1
        segment = Segment.from_chromosome_coordinates(chromosome=chromosome, start=start_coord, end=end_coord)
        segments.append(segment)
        fid = segment.stable_id_non_hap
        for clone_id in clone_ids:
            cns_str = clone_cn_strs[clone_id_mappings[clone_id]]
            data = cns_str.split("|")
            cna = int(data[0])
            cnb = int(data[1])
            scnt[clone_id].set_cn_record(sid=fid, hap=Haplotype.A, cn=cna)
            scnt[clone_id].set_cn_record(sid=fid, hap=Haplotype.B, cn=cnb)
    return segments, scnt


REMIXT_CHROMOSOME = "chromosome"
REMIXT_START_POSITION = "start"
REMIXT_END_POSITION = "end"
REMIXT_CLONE1_CN_A = "major_1"
REMIXT_CLONE1_CN_B = "minor_1"
REMIXT_CLONE2_CN_A = "major_2"
REMIXT_CLONE2_CN_B = "minor_2"


def get_scnt_from_remixt_file(file_name, separator="\t", chr_strip=True):
    with open(file_name, "rt") as source:
        return get_scnt_from_remixt_source(source=source, separator=separator, chr_strip=chr_strip)


def get_scnt_from_remixt_source(source, separator="\t", chr_strip=True):
    segments = []
    clone1_id = "1"
    clone2_id = "2"
    scnt = {clone1_id: SegmentCopyNumberProfile(), clone2_id: SegmentCopyNumberProfile()}
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        chromosome = row[REMIXT_CHROMOSOME]
        if chr_strip:
            chromosome = strip_chr(chr_string=chromosome)
        start_coordinate = int(row[REMIXT_START_POSITION])
        end_coordinate = int(row[REMIXT_END_POSITION]) - 1
        segment = Segment.from_chromosome_coordinates(chromosome=chromosome, start=start_coordinate, end=end_coordinate)
        segments.append(segment)
        sid = segment.stable_id_non_hap
        clone_1_cn_a = int(row[REMIXT_CLONE1_CN_A])
        clone_1_cn_b = int(row[REMIXT_CLONE1_CN_B])
        clone_2_cn_a = int(row[REMIXT_CLONE2_CN_A])
        clone_2_cn_b = int(row[REMIXT_CLONE2_CN_A])
        clone1_scnp = scnt[clone1_id]
        clone2_scnp = scnt[clone2_id]
        clone1_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=clone_1_cn_a)
        clone1_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=clone_1_cn_b)
        clone2_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=clone_2_cn_a)
        clone2_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=clone_2_cn_b)
    return segments, scnt


TITAN_CHROMOSOME = "Chromosome"
TITAN_START_POSITION = "Start"
TITAN_END_POSITION = "End"
TITAN_MAJOR_CN = "MajorCN"
TITAN_MINOR_CN = "MinorCN"
TITAN_CLONE_ID = "Clonal_Cluster"
TITAN_CORRECTED_CN = "Corrected_Copy_Number"
TITAN_SAMPLE_NAME = "Sample"


def titan_get_clone_ids_from_file(file_name, sample_name, separator="\t"):
    with open(file_name, "rt") as source:
        result = set()
        reader = csv.DictReader(source, delimiter=separator)
        for row in reader:
            if row[TITAN_SAMPLE_NAME] != sample_name:
                continue
            clone_id = row[TITAN_CLONE_ID]
            if clone_id != "NA":
                result.add(clone_id)
        return sorted(result)


def get_scnt_from_titan_file(file_name, sample_name, clone_ids=None, separator="\t", corrected_cn_fix="None", chr_strip=True):
    if clone_ids is None:
        clone_ids = titan_get_clone_ids_from_file(file_name=file_name, sample_name=sample_name, separator=separator)
    with open(file_name, "rt") as source:
        return get_scnt_from_titan_source(source=source, sample_name=sample_name, clone_ids=clone_ids, separator=separator, corrected_cn_fix=corrected_cn_fix, chr_strip=chr_strip)


def get_scnt_from_titan_source(source, sample_name, clone_ids, separator="\t", corrected_cn_fix="None", chr_strip=True):
    scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    segments = []
    reader = csv.DictReader(source, delimiter=separator)
    for row in reader:
        if row[TITAN_SAMPLE_NAME] != sample_name:
            continue
        chromosome = row[TITAN_CHROMOSOME]
        if chr_strip:
            chromosome = strip_chr(chr_string=chromosome)
        segment = Segment.from_chromosome_coordinates(chromosome=chromosome, start=int(row[TITAN_START_POSITION]), end=int(row[TITAN_END_POSITION]))
        sid = segment.stable_id_non_hap
        segments.append(segment)
        major_cn, minor_cn = int(row[TITAN_MAJOR_CN]), int(row[TITAN_MINOR_CN])
        if minor_cn > major_cn:
            minor_cn, major_cn = major_cn, minor_cn
        titan_clone_id = row[TITAN_CLONE_ID]
        corrected_cn = int(row[TITAN_CORRECTED_CN])
        for clone_id in clone_ids:
            scnp = scnt[clone_id]
            if titan_clone_id == clone_id:
                if major_cn + minor_cn != corrected_cn and corrected_cn_fix != "None":
                    diff = corrected_cn - major_cn - minor_cn
                    ###
                    # initialize as 0 when corrected_cn_fix strategy does not match any known, yet is not "None"
                    ###
                    major_cn_addition = 0
                    minor_cn_addition = 0
                    if corrected_cn_fix == "equal":
                        major_cn_addition = int(math.ceil(diff / 2))
                        minor_cn_addition = diff - major_cn_addition
                    elif corrected_cn_fix == "relative-dist":
                        relative_relation = minor_cn * 1.0 / major_cn
                        major_cn_addition = int(math.ceil(diff / (1 + relative_relation)))
                        minor_cn_addition = diff - major_cn_addition
                    major_cn += major_cn_addition
                    minor_cn += minor_cn_addition
                scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=major_cn)
                scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=minor_cn)
            else:
                scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=1)
                scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=1)
    return segments, scnt
