import argparse

import os

import sys


current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import read_scnt_from_source, extract_scnb_from_segments, write_scnb_to_destination
from rck.core.structures import SCNBoundariesStrategies, LengthSpreadRelationships, SegmentCopyNumberBoundaries


def main():
    parser = argparse.ArgumentParser("Creating boundaries for a RCK formatted segment copy number tensor")
    parser.add_argument('--version', action='version', version=rck.version)
    parser.add_argument("scnt", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("--bnd-strategy", choices=[strategy.value for strategy in SCNBoundariesStrategies], type=SCNBoundariesStrategies.from_string,
                        default=SCNBoundariesStrategies.UNIFORM_MIN_MAX.value)
    parser.add_argument("--uniform-spread-size", type=int, default=1)
    parser.add_argument("--length-spread-relation", choices=[rel.value for rel in LengthSpreadRelationships], type=LengthSpreadRelationships.from_string,
                        default=LengthSpreadRelationships.DUMMY.value)
    parser.add_argument("--uniform-min", type=int, default=0)
    parser.add_argument("--uniform-max", type=int, default=10)
    parser.add_argument("--missing-only", action="store_true", dest="missing_only")
    parser.add_argument("--min-allow-zero-for-positive", type=int, default=-1)
    parser.add_argument("--max-allow-zero-for-positive", type=int, default=1000000000)
    parser.add_argument("--min-allow-positive-for-zero", type=int, default=-1)
    parser.add_argument("--max-allow-positive-for-zero", type=int, default=1000000000)
    parser.add_argument("--clone-ids", default="")
    parser.add_argument("--is-male", action="store_false", dest="is_female")
    parser.add_argument("--no-allow-unit-segments", dest="allow_unit_segments", action="store_false")
    parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--o-with-scnt", action="store_true", dest="output_scnt")
    parser.add_argument("--separator", default="\t")
    args = parser.parse_args()
    clone_ids = args.clone_ids.split(",") if len(args.clone_ids) > 0 else None
    segments, scnt = read_scnt_from_source(source=args.scnt, clone_ids=clone_ids)
    if clone_ids is None:
        clone_ids = sorted(scnt.keys())
    try:
        scnb = extract_scnb_from_segments(segments=segments, clone_ids=clone_ids)
    except ValueError:
        scnb = {clone_id: SegmentCopyNumberBoundaries() for clone_id in clone_ids}
    for clone_id in clone_ids:
        scnb[clone_id].fill(segments=segments, scnp=scnt[clone_id], missing_only=args.missing_only, strategy=args.bnd_strategy,
                            min_allow_zero_for_positive=args.min_allow_zero_for_positive,
                            max_allow_zero_for_positive=args.max_allow_zero_for_positive,
                            min_allow_positive_for_zero=args.min_allow_positive_for_zero,
                            max_allow_positive_for_zero=args.max_allow_positive_for_zero,
                            uniform_spread_size=args.uniform_spread_size,
                            length_spread_relation=args.length_spread_relation,
                            uniform_min=args.uniform_min,
                            uniform_max=args.uniform_max,
                            is_female=args.is_female)
    write_scnb_to_destination(destination=args.output, segments=segments, scnb=scnb, clone_ids=clone_ids)


if __name__ == "__main__":
    main()
