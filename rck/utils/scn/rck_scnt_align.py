import argparse

import os

import sys
from copy import deepcopy

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from rck.core.io import read_scn_tensor, get_all_clone_ids_from_dasck_scnt_file, write_scn_tensor
from rck.core.structures import align_scnts


def main():
    parser = argparse.ArgumentParser("Align multiple Segment Copy Number Tensors")
    parser.add_argument("scnt", nargs="+")
    parser.add_argument("--output-suffix", default="aligned")
    parser.add_argument("--no-allow-unit-segments", action="store_false", dest="allow_unit_segments")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()
    scnt_files = {}
    for path in args.scnt:
        full_path = os.path.abspath(os.path.expanduser(path))
        name = os.path.splitext(os.path.basename(full_path))[0]
        if name.endswith(".scnt"):
            name = name[:-5]
        if name.endswith("."):
            name = name[:-1]
        scnt_files[name] = full_path
    scnts_by_name = {}
    segments_by_name = {}
    clone_ids = {}
    for name, path in scnt_files.items():
        clone_ids[name] = get_all_clone_ids_from_dasck_scnt_file(file_name=scnt_files[name], allow_unit_segments=args.allow_unit_segments)
        segments, scnt = read_scn_tensor(file_name=scnt_files[name], clone_ids=clone_ids[name], allow_unit_segments=args.allow_unit_segments)
        scnts_by_name[name] = scnt
        segments_by_name[name] = segments
    if len(scnts_by_name.values()) == 1:
        aligned_segments_by_name, aligned_scnts_by_name = deepcopy(segments_by_name), deepcopy(scnts_by_name)
    else:
        aligned_segments_by_name, aligned_scnts_by_name = align_scnts(segments_by_sample_names=segments_by_name, scnts_by_sample_names=scnts_by_name)
    result_base_names = {}
    cnt = 0
    for name in sorted(scnt_files.keys()):
        new_name = name
        if name in result_base_names:
            new_name = name + str(cnt)
            cnt += 1
        new_name = new_name + "." + args.output_suffix
        result_base_names[name] = new_name
    output_dir = args.output_dir if args.output_dir != "" else os.getcwd()
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    for name, new_name in result_base_names.items():
        scnt = aligned_scnts_by_name[name]
        segments = aligned_segments_by_name[name]
        scnt_path = os.path.join(output_dir, new_name + ".scnt.tsv")
        write_scn_tensor(file_name=scnt_path, scnt=scnt, segments=segments)


if __name__ == "__main__":
    main()
