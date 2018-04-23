import csv

from collections import defaultdict

from dassp.core.structures import Position, Strand, Adjacency, AdjacencyType

SAMPLE_NAME = "patient"
CHR1 = "chr1"
CHR2 = "chr2"
COORDINATE1 = "pos1"
COORDINATE2 = "pos2"
STRAND1 = "str1"
STRAND2 = "str2"
SAMPLES_LIST = "samples"
SAMPLES_COUNTS = "tumour-counts"


def flip_strand(strand_string):
    if strand_string == "-":
        return "+"
    elif strand_string == "+":
        return "-"
    raise Exception()


def read_breakpoints_into_novel_adjacencies(file_name, sample_basename, sample_name_suffix,
                                            min_cnt=1, flip_second_strand=False, min_distance=1000, remove_adjs_sharing_exact_positions=True):
    result = {}
    positions_to_adjacencies = defaultdict(list)
    with open(file_name, "rt") as csv_source:
        reader = csv.DictReader(csv_source)
        for row in reader:
            if row[SAMPLE_NAME] != sample_basename:
                continue
            samples_list = row[SAMPLES_LIST].split("/")
            samples_cnts = list(map(int, row[SAMPLES_COUNTS].split("/")))
            sample_index = samples_list.index(sample_name_suffix)
            if sample_index < 0:
                continue
            sample_cnt = samples_cnts[sample_index]
            if sample_cnt < min_cnt:
                continue
            position1 = Position(chromosome=row[CHR1], coordinate=int(row[COORDINATE1]), strand=Strand.from_pm_string(row[STRAND1]))
            strand2_string = flip_strand(row[STRAND2]) if flip_second_strand else row[STRAND2]
            position2 = Position(chromosome=row[CHR2], coordinate=int(row[COORDINATE2]), strand=Strand.from_pm_string(strand2_string))
            if min_distance > 0 and position1.chromosome == position2.chromosome:
                left_position = min([position1.coordinate, position2.coordinate])
                right_position = max([position1.coordinate, position2.coordinate])
                if right_position - left_position < min_distance:
                    continue
            novel_adjacency = Adjacency(position1=position1, position2=position2, adjacency_type=AdjacencyType.NOVEL)
            positions_to_adjacencies[position1].append(novel_adjacency)
            positions_to_adjacencies[position2].append(novel_adjacency)
            result[novel_adjacency.stable_id_non_phased] = novel_adjacency
    if remove_adjs_sharing_exact_positions:
        for position, adjacencies in positions_to_adjacencies.items():
            if len(adjacencies) > 1:
                for adj in adjacencies:
                    if adj.stable_id_non_phased in result:
                        del result[adj.stable_id_non_phased]
    return result
