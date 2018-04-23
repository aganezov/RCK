import argparse
import csv

import os

import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from dassp.core.structures import Position, Strand, Adjacency, AdjacencyType
from dassp.core.io import EXTERNAL_NA_ID


def allows_field(field_value, allowed_list):
    if len(allowed_list) == 0:
        return True
    return field_value in allowed_list


def flip_strand_string(strand_string):
    if strand_string == "-":
        return "+"
    return "-"


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Converter from Gundem2015 breakpoint prediction results to DASCK novel adjacencies input format")
    parser.add_argument('gundem15_na_file')
    parser.add_argument("--patients", default="")
    parser.add_argument("--delimiter", default="\t")
    parser.add_argument("--sample", default="")
    parser.add_argument("--min-sample-cnt", type=int, default=1)
    parser.add_argument("--no-flip-second-strand", action="store_false", dest="flip_second_strand")
    parser.add_argument("--brass2", choices=['all', 'brass2', 'no-brass2'], default='all')
    parser.add_argument("--na-id-label", choices=['naid', 'prediction_id'], default='prediction_id')
    parser.add_argument('dasck_na_file')
    args = parser.parse_args()
    patients = args.patients.split(",") if len(args.patients) > 0 else []
    samples = args.sample.split(",") if len(args.sample) > 0 else []
    brass2 = []
    if args.brass2 == 'brass2':
        brass2 = ['brass2']
    elif args.brass2 == 'no-brass2':
        brass2 = ['no-brass2']
    gundem15_na_file = os.path.expanduser(args.gundem15_na_file)
    gundem15_na_file = os.path.abspath(gundem15_na_file)
    na_id = 1
    adjacencies = []
    with open(gundem15_na_file, "rt") as source:
        reader = csv.DictReader(source, delimiter=args.delimiter)
        for line in reader:
            patient = line['patient']
            if not allows_field(field_value=patient, allowed_list=patients):
                continue
            brass2_flag = line['brass2']
            if not allows_field(field_value=brass2_flag, allowed_list=brass2):
                continue
            chr1_name = line['chr1']
            chr1_coordinate = int(line['pos1'])
            str1 = Strand.from_pm_string(string=line['str1'])
            chr2_name = line['chr2']
            chr2_coordinate = int(line['pos2'])
            str2 = line['str2']
            if args.flip_second_strand:
                str2 = flip_strand_string(strand_string=str2)
            str2 = Strand.from_pm_string(string=str2)
            samples_entries = line['samples'].split("/")
            sample_counts = list(map(lambda e: int(e), line['tumour-counts'].split("/")))
            allowed = False
            for sample_name, sample_cnt in zip(samples_entries, sample_counts):
                if allows_field(sample_name, samples) and sample_cnt >= args.min_sample_cnt:
                    allowed = True
            if not allowed:
                continue
            p1 = Position(chromosome=chr1_name, coordinate=chr1_coordinate, strand=str1)
            p2 = Position(chromosome=chr2_name, coordinate=chr2_coordinate, strand=str2)
            adjacency = Adjacency(position1=p1, position2=p2, adjacency_type=AdjacencyType.NOVEL, extra={EXTERNAL_NA_ID: na_id})
            adjacencies.append(adjacency)
            na_id += 1
    with open(args.dasck_na_file, "wt") as destination:
        print(args.na_id_label, "chromosome_1", "position_1", "strand_1", "chromosome_2", "position_2", "strand_2", sep="\t", file=destination)
        for adjacency in adjacencies:
            na_strings_entries = []
            na_strings_entries.append("{na_id}".format(na_id=adjacency.extra[EXTERNAL_NA_ID]))
            na_strings_entries.append("{chr1}".format(chr1=adjacency.position1.chromosome))
            na_strings_entries.append("{coord1}".format(coord1=adjacency.position1.coordinate))
            na_strings_entries.append("{strand1}".format(strand1=adjacency.position1.strand))
            na_strings_entries.append("{chr2}".format(chr2=adjacency.position2.chromosome))
            na_strings_entries.append("{coord2}".format(coord2=adjacency.position2.coordinate))
            na_strings_entries.append("{strand2}".format(strand2=adjacency.position2.strand))
            na_string_entry = "\t".join(na_strings_entries)
            print(na_string_entry, file=destination)
