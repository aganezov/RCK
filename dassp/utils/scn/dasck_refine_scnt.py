import argparse
import csv

import os

import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from dassp.core.io import write_scn_tensor, get_all_clone_ids_from_dasck_scnt_file, read_scn_tensor
from dassp.core.structures import SegmentCopyNumberProfile, Position, Strand, Segment, Haplotype, refined_scnt

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Refiner of DASCK formatted segment copy number tensor")
    parser.add_argument('dasck_scn_file_input')
    parser.add_argument("--no-allow-missed-clones", action="store_false", dest="allow_missed_clones")
    parser.add_argument("--clone-ids", default="")
    parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    parser.add_argument("--max-merge-gap", type=int, default=1000000)
    parser.add_argument("--no-fill-gaps", action="store_false", dest="fill_gaps")
    parser.add_argument("--max-fill-gap", type=int, default=1000000)
    parser.add_argument('dasck_scn_file_output')
    args = parser.parse_args()
    clone_ids = args.clone_ids.split(",") if len(args.clone_ids) > 0 else []
    dasck_scnt_file_path = os.path.expanduser(args.dasck_scn_file_input)
    dasck_scnt_file_path = os.path.abspath(dasck_scnt_file_path)
    if len(clone_ids) == 0:
        clone_ids = get_all_clone_ids_from_dasck_scnt_file(file_name=dasck_scnt_file_path)
    segments, scnt = read_scn_tensor(file_name=dasck_scnt_file_path, clone_ids=clone_ids)
    if args.merge_fragments or args.fill_gaps:
        segments, scnt = refined_scnt(segments=segments, tensor=scnt,
                                      merge_fragments=args.merge_fragments, max_merge_gap=args.max_merge_gap,
                                      fill_gaps=args.fill_gaps, max_fill_gap=args.max_fill_gap)
    write_scn_tensor(file_name=args.dasck_scn_file_output, scnt=scnt, segments=segments)
