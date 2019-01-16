import argparse
import os

import sys
from copy import deepcopy

current_file_level = 2
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import read_adjacencies_from_file, \
    get_logging_cli_parser, get_standard_logger_from_args, get_full_path, read_scnt_from_file, read_positions_from_source, \
    write_segments_to_file, write_scnt_to_file
from rck.core.structures import refined_scnt, refined_scnt_with_adjacencies_and_telomeres


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-input-refine", parents=[get_logging_cli_parser()])
    parser.add_argument("--version", action="version", version=rck.version)
    parser.add_argument("--scnt", required=True)
    parser.add_argument("--adjacencies", required=True)
    parser.add_argument("--clone-ids", default=None)
    parser.add_argument("--scnt-separator", default="\t")
    parser.add_argument("--adjacencies-separator", default="\t")
    parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    parser.add_argument("--fragments-max-merge-gap", type=int, default=1000000000)
    parser.add_argument("--no-fill-gaps-fragments", action="store_false", dest="fill_gaps_fragments")
    parser.add_argument("--fragments-max-fill-gap", type=int, default=1000000000)
    parser.add_argument("--no-allow-unit-segments", action="store_false", dest="allow_unit_segments")
    parser.add_argument("--telomere-positions", type=argparse.FileType("rt"))
    parser.add_argument("--telomere-positions-separator", default="\t")
    parser.add_argument("--output-scnt", required=True)
    parser.add_argument("--output-fragments", required=True)
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-input-refine")
    clone_ids = args.clone_ids.split(",") if args.clone_ids is not None else None
    scnt_file = get_full_path(args.scnt_file)
    adj_file = get_full_path(args.adj)
    segments, scnt = read_scnt_from_file(file_name=scnt_file, clone_ids=clone_ids, separator=args.scnt_separator)
    clone_ids = sorted(set(scnt.keys()))
    segments, scnt, segments_ids_mapping = refined_scnt(segments=segments, scnt=scnt,
                                                        merge_fragments=args.merge_fragments, max_merge_gap=args.fragments_max_merge_gap,
                                                        fill_gaps=args.fill_gaps_fragments, max_fill_gap=args.fragments_max_fill_gap)

    adjacencies = read_adjacencies_from_file(file_name=adj_file, separator=args.adjacencies_separator)
    if args.telomere_positions is not None:
        telomere_positions = read_positions_from_source(source=args.telomere_positions, separator=args.telomere_positions_separator)
    else:
        telomere_positions = []
    fragments = deepcopy(segments)
    segments, scnt = refined_scnt_with_adjacencies_and_telomeres(segments=segments, scnt=scnt, adjacencies=adjacencies, telomere_positions=telomere_positions)
    refined_scnt_file = os.path.expanduser(args.refined_scnt_file)
    refined_scnt_file = os.path.abspath(refined_scnt_file)
    fragments_file = get_full_path(path=args.output_fragments)

    write_segments_to_file(file_name=fragments_file, segments=fragments)
    write_scnt_to_file(file_name=refined_scnt_file, scnt=scnt, segments=segments)


if __name__ == "__main__":
    main()
