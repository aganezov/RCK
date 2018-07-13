import argparse
import csv

import os

import sys

current_file_level = 4
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)


from dassp.core.structures import SegmentCopyNumberProfile, Position, Strand, Segment, Haplotype, refined_scnt
from dassp.core.io import write_scn_tensor


def get_clone_ids(hatchet_file_path, sample, separator="\t", min_usage=0.01):
    result = set()
    candidates = []
    with open(hatchet_file_path, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            data = line.split(separator)
            clone_data = data[6:]
            if line_cnt == 0:
                total_clone_cnt = int(len(clone_data) / 2)
                candidates = [str(cnt) for cnt in range(1, total_clone_cnt + 1)]
            if line.startswith("#"):
                continue
            sample_name = data[3]
            if sample_name != sample:
                continue
            for candidate_clone_id, clone_usage_str in zip(candidates, clone_data[1::2]):
                clone_usage = float(clone_usage_str)
                if clone_usage < min_usage:
                    continue
                result.add(candidate_clone_id)
            if sorted(result) == candidates:
                return sorted(result)
    return sorted(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hatchet_scnt")
    parser.add_argument("--sample", type=str, required=True)
    parser.add_argument("--clone_ids", type=str, default="")
    parser.add_argument("--min-usage", type=float, default=0.01)
    parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    parser.add_argument("--max-merge-gap", type=int, default=10000000)
    parser.add_argument("--no-fill-gaps", action="store_false", dest="fill_gaps")
    parser.add_argument("--max-fill-gap", type=int, default=10000000)
    parser.add_argument("dasck_scnt")
    args = parser.parse_args()
    hatchet_file_path = os.path.abspath(os.path.expanduser(args.hatchet_scnt))
    dasck_file_path = os.path.abspath(os.path.expanduser(args.dasck_scnt))
    clone_ids = get_clone_ids(hatchet_file_path=hatchet_file_path, sample=args.sample)
    scnt = {clone_id: SegmentCopyNumberProfile() for clone_id in clone_ids}
    segments = []
    with open(hatchet_file_path, "rt") as source:
        clone_id_mappings = {}
        for line_cnt, line in enumerate(source):
            line = line.strip()
            data = line.split("\t")
            clone_data = data[6:]
            if line_cnt == 0:
                total_clone_cnt = int(len(clone_data) / 2)
                candidates = [str(cnt) for cnt in range(1, total_clone_cnt + 1)]
                for position_cnt, candidate in enumerate(candidates):
                    if candidate in clone_ids:
                        clone_id_mappings[candidate] = position_cnt
            clone_cn_strs = clone_data[::2]
            if line.startswith("#") or len(line) == 0:
                continue
            sample_name = data[3]
            if sample_name != args.sample:
                continue
            chr_name = data[0]
            start_coord = int(data[1])
            end_coord = int(data[2]) - 1
            sp = Position(chromosome=chr_name, coordinate=start_coord, strand=Strand.REVERSE)
            ep = Position(chromosome=chr_name, coordinate=end_coord, strand=Strand.FORWARD)
            fragment = Segment(start_position=sp, end_position=ep)
            segments.append(fragment)
            fid = fragment.stable_id_non_hap
            for clone_id in clone_ids:
                cns_str = clone_cn_strs[clone_id_mappings[clone_id]]
                data = cns_str.split("|")
                cna = int(data[0])
                cnb = int(data[1])
                scnt[clone_id].set_cn_record(sid=fid, hap=Haplotype.A, cn=cna)
                scnt[clone_id].set_cn_record(sid=fid, hap=Haplotype.B, cn=cnb)
    if args.merge_fragments or args.fill_gaps:
        segments, scnt = refined_scnt(segments=segments, tensor=scnt,
                                      merge_fragments=args.merge_fragments, max_merge_gap=args.max_merge_gap,
                                      fill_gaps=args.fill_gaps, max_fill_gap=args.max_fill_gap)
    write_scn_tensor(file_name=dasck_file_path, scnt=scnt, segments=segments)
