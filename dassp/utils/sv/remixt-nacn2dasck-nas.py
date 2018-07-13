import argparse
import csv

import os

import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)


from dassp.core.io import EXTERNAL_NA_ID, write_adjacencies
from dassp.core.structures import Strand, Position, Adjacency, AdjacencyType


def allows_field(field_value, allowed_list):
    if len(allowed_list) == 0:
        return True
    return field_value in allowed_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Converter from ReMixT breakpoint copy number predictions to DASCK novel adjacencies input (loosing copy numbers)")
    parser.add_argument("remixt_nacn_file")
    parser.add_argument("--clone-ids", choices=["1", "2", "1,2"], default="1,2")
    parser.add_argument("--no-skip-absent", action="store_false", dest="skip_absent")
    parser.add_argument("--no-remixt-na-correction", action="store_false", dest="remixt_correction")
    parser.add_argument("dasck_na_file")
    args = parser.parse_args()
    clone_ids = args.clone_ids.split(",")
    adjacencies = []
    remixt_nacn_file = os.path.expanduser(args.remixt_nacn_file)
    remixt_nacn_file = os.path.abspath(remixt_nacn_file)
    with open(remixt_nacn_file, "rt") as source:
        reader = csv.DictReader(source, delimiter="\t")
        for line in reader:
            naid = line["prediction_id"]
            chr1_name = line["chromosome_1"]
            chr1_coordinate = int(line["position_1"])
            chr1_strand = Strand.from_pm_string(line["strand_1"])
            chr2_name = line["chromosome_2"]
            chr2_coordinate = int(line["position_2"])
            chr2_strand = Strand.from_pm_string(line["strand_2"])
            clone1_cn = int(line["cn_1"])
            clone2_cn = int(line["cn_2"])
            skip = True
            for clone_cn, clone_name in zip([clone1_cn, clone2_cn], ["1", "2"]):
                if clone_name in clone_ids and clone_cn > 0:
                    skip = False
            skip = skip and args.skip_absent
            if skip:
                continue
            if args.remixt_correction and chr1_strand == Strand.FORWARD:
                chr1_coordinate -= 1
            if args.remixt_correction and chr2_strand == Strand.FORWARD:
                chr2_coordinate -= 1
            position1 = Position(chromosome=chr1_name, coordinate=chr1_coordinate, strand=chr1_strand)
            position2 = Position(chromosome=chr2_name, coordinate=chr2_coordinate, strand=chr2_strand)
            adjacency = Adjacency(position1=position1, position2=position2,
                                  adjacency_type=AdjacencyType.NOVEL, extra={EXTERNAL_NA_ID: naid})
            adjacencies.append(adjacency)
    dasck_na_file = os.path.expanduser(args.dasck_na_file)
    dasck_na_file = os.path.abspath(dasck_na_file)
    write_adjacencies(file_name=dasck_na_file, adjacencies=adjacencies)
