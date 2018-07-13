import argparse
import csv

import os

import sys

from core.structures import LengthSpreadRelationships

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)


from dassp.core.io import read_scn_tensor, get_all_clone_ids_from_dasck_scnt_file, write_scn_boundaries
from dassp.core.structures import get_scn_boundaries, SCNBoundariesStrategies


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Creating boundaries for a DASCK formatted segment copy number tensor")
    parser.add_argument("dasck_input_scnt_file")
    parser.add_argument("--boundaries-strategy", choices=["fixed", "uniform-spread", "length-spread", "uniform-min-max"], default="length-spread")
    parser.add_argument("--uniform-spread-size", type=int, default=1)
    parser.add_argument("--length-spread-relation", choices=["dummy"], default="dummy")
    parser.add_argument("--uniform-min", type=int, default=0)
    parser.add_argument("--uniform-max", type=int, default=10)
    parser.add_argument("--min-allow-zero-for-positive", type=int, default=-1)
    parser.add_argument("--max-allow-zero-for-positive", type=int, default=1000000000)
    parser.add_argument("--min-allow-positive-for-zero", type=int, default=-1)
    parser.add_argument("--max-allow-positive-for-zero", type=int, default=1000000000)
    parser.add_argument("--clone-ids", default="")
    parser.add_argument("--is-male", action="store_false", dest="is_female")
    parser.add_argument("--no-allow-unit-segments", dest="allow_unit_segments", action="store_false")
    parser.add_argument("scn_boundaries_file")
    args = parser.parse_args()
    clone_ids = args.clone_ids.split(",") if len(args.clone_ids) > 0 else []
    dasck_input_scnt_filepath = os.path.expanduser(args.dasck_input_scnt_file)
    dasck_input_scnt_filepath = os.path.abspath(dasck_input_scnt_filepath)
    strategy = SCNBoundariesStrategies.from_str_value(string_value=args.boundaries_strategy)
    length_spread_relation = LengthSpreadRelationships.DUMMY
    if len(clone_ids) == 0:
        clone_ids = get_all_clone_ids_from_dasck_scnt_file(file_name=dasck_input_scnt_filepath, allow_unit_segments=args.allow_unit_segments)
    segments, scnt = read_scn_tensor(file_name=dasck_input_scnt_filepath, clone_ids=clone_ids, allow_unit_segments=args.allow_unit_segments)
    boundaries = get_scn_boundaries(segments=segments, scnt=scnt, strategy=strategy,
                                    min_allow_zero_for_positive=args.min_allow_zero_for_positive,
                                    max_allow_zero_for_positive=args.max_allow_zero_for_positive,
                                    min_allow_positive_for_zero=args.min_allow_positive_for_zero,
                                    max_allow_positive_for_zero=args.max_allow_positive_for_zero,
                                    uniform_spread_size=args.uniform_spread_size,
                                    length_spread_relation=length_spread_relation,
                                    uniform_min=args.uniform_min,
                                    uniform_max=args.uniform_max,
                                    is_female=args.is_female)
    write_scn_boundaries(file_name=args.scn_boundaries_file, segments=segments, scn_boundaries=boundaries)
