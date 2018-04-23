import argparse

import os

import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)


from dassp.core.io import read_novel_adjacencies, write_adjacencies
from dassp.core.structures import refined_nas, violates_is, removed_short_nas

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Refiner of DASCK formatted novel adjacencies")
    parser.add_argument("dasck_nas_input")
    parser.add_argument("--no-remove-short", action="store_false", dest="remove_short")
    parser.add_argument("--short-min-size", type=int, default=500)
    parser.add_argument("--no-reciprocal-merging", action="store_false", dest="merge_reciprocal")
    parser.add_argument("--max-reciprocal-window", type=int, default=50)
    parser.add_argument("--no-check-isv", action="store_false", dest="check_is_violations")
    parser.add_argument("dasck_nas_output")
    args = parser.parse_args()
    dasck_nas_input_filepath = os.path.expanduser(args.dasck_nas_input)
    dasck_nas_input_filepath = os.path.abspath(dasck_nas_input_filepath)
    nas = read_novel_adjacencies(file_name=dasck_nas_input_filepath)
    if args.remove_short:
        nas = removed_short_nas(novel_adjacencies=nas, min_size=args.short_min_size)
    if args.merge_reciprocal:
        nas = refined_nas(novel_adjacencies=nas, window_size=args.max_reciprocal_window)
    if args.check_is_violations and violates_is(novel_adjacencies=nas):
        raise Exception()
    write_adjacencies(file_name=args.dasck_nas_output, adjacencies=nas)
