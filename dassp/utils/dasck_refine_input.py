import argparse
import os

import sys
from copy import deepcopy

current_file_level = 2
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)


from dassp.core.io import read_scn_tensor, read_novel_adjacencies, get_all_clone_ids_from_dasck_scnt_file, read_positions, write_segments, write_scn_tensor
from dassp.core.structures import refined_scnt, removed_short_nas, refined_nas, violates_is, refined_scnt_with_nas_and_telomeres

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scnt_file")
    parser.add_argument("nas_file")
    parser.add_argument("--clone-ids", default="")
    parser.add_argument("--no-remove-short-nas", action="store_false", dest="remove_short_nas")
    parser.add_argument("--nas-short-min-size", type=int, default=500)
    parser.add_argument("--no-reciprocal-merging-nas", action="store_false", dest="merge_reciprocal_nas")
    parser.add_argument("--nas-max-reciprocal-window", type=int, default=50)
    parser.add_argument("--no-check-isv", action="store_false", dest="check_is_violations")
    parser.add_argument("--no-allow-missed-clones", action="store_false", dest="allow_missed_clones")
    parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    parser.add_argument("--fragments-max-merge-gap", type=int, default=1000000)
    parser.add_argument("--no-fill-gaps-fragments", action="store_false", dest="fill_gaps_fragments")
    parser.add_argument("--fragments-max-fill-gap", type=int, default=1000000)
    parser.add_argument("--telomeres-file", default="")
    parser.add_argument("--tel-human-centromeres", action="store_true", dest="tel_use_centromeres")
    parser.add_argument("--ex-allow-unit-segments", action="store_true", dest="allow_unit_segments")
    parser.add_argument("refined_scnt_file")
    parser.add_argument("refined_fragments_file")
    args = parser.parse_args()
    clone_ids = args.clone_ids.split(",") if len(args.clone_ids) > 0 else []
    dasck_scnt_input_file_path = os.path.expanduser(args.scnt_file)
    dasck_scnt_input_file_path = os.path.abspath(dasck_scnt_input_file_path)
    if len(clone_ids) == 0:
        clone_ids = get_all_clone_ids_from_dasck_scnt_file(file_name=dasck_scnt_input_file_path,
                                                           allow_unit_segments=args.allow_unit_segments)
    segments, scnt = read_scn_tensor(file_name=dasck_scnt_input_file_path, clone_ids=clone_ids,
                                     allow_unit_segments=args.allow_unit_segments)
    if args.merge_fragments or args.fill_gaps_fragments:
        segments, scnt = refined_scnt(segments=segments, tensor=scnt,
                                      merge_fragments=args.merge_fragments, max_merge_gap=args.fragments_max_merge_gap,
                                      fill_gaps=args.fill_gaps_fragments, max_fill_gap=args.fragments_max_fill_gap)
    dasck_nas_input_filepath = os.path.expanduser(args.nas_file)
    dasck_nas_input_filepath = os.path.abspath(dasck_nas_input_filepath)
    nas = read_novel_adjacencies(file_name=dasck_nas_input_filepath)
    if args.remove_short_nas:
        nas = removed_short_nas(novel_adjacencies=nas, min_size=args.nas_short_min_size)
    if args.merge_reciprocal_nas:
        nas = refined_nas(novel_adjacencies=nas, window_size=args.nas_max_reciprocal_window)
    if args.check_is_violations and violates_is(novel_adjacencies=nas):
        raise Exception()
    if args.telomeres_file != "":
        dasck_telomeres_input_filepath = os.path.expanduser(args.telomeres_file)
        dasck_telomeres_input_filepath = os.path.abspath(dasck_telomeres_input_filepath)
        telomeres = read_positions(file_name=dasck_telomeres_input_filepath)
    else:
        telomeres = None
    fragments = deepcopy(segments)
    segments, scnt = refined_scnt_with_nas_and_telomeres(segments=segments, scnt=scnt, novel_adjacencies=nas, telomeres=telomeres)
    refined_scnt_file = os.path.expanduser(args.refined_scnt_file)
    refined_scnt_file = os.path.abspath(refined_scnt_file)
    fragments_file = os.path.expanduser(args.refined_fragments_file)
    fragments_file = os.path.abspath(fragments_file)
    write_segments(file_name=fragments_file, segments=fragments)
    write_scn_tensor(file_name=refined_scnt_file, scnt=scnt, segments=segments)
