import argparse
import csv

import sys
import os

import sys

current_file_level = 4
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from dassp.core.io import write_scn_tensor
from dassp.core.structures import SegmentCopyNumberProfile, Position, Strand, Segment, Haplotype, refined_segments_and_tensor

CLONE1_CN_A = "nMaj1_A"
CLONE1_CN_B = "nMin1_A"
CLONE2_CN_A = "nMaj2_A"
CLONE2_CN_B = "nMin2_A"
START_POSITION = "startpos"
END_POSITION = "endpos"
SAMPLE_NAME = "sample"
CHROMOSOME = "chr"


def force_non_negativity(cn):
    return cn if cn >= 0 else 0


def get_subclonal_cn(subclonal_cn_string, clonal_cn_int):
    if subclonal_cn_string == "NA":
        return clonal_cn_int
    return force_non_negativity(int(subclonal_cn_string))


def read_battenberg_scv_file(file_name, sample_name=None):
    clone1_name = '0'
    clone2_name = '1'
    clone_specific_fcnp = {clone1_name: SegmentCopyNumberProfile(), clone2_name: SegmentCopyNumberProfile()}
    hapl_fragments = {}
    with open(file_name, "rt") as csv_source:
        reader = csv.DictReader(csv_source)
        for row in reader:
            if sample_name is not None and row[SAMPLE_NAME] != sample_name:
                continue
            start_position = Position(chromosome=row[CHROMOSOME], coordinate=int(row[START_POSITION]), strand=Strand.REVERSE)
            end_position = Position(chromosome=row[CHROMOSOME], coordinate=int(row[END_POSITION]), strand=Strand.FORWARD)
            fragment = Segment(start_position=start_position, end_position=end_position)
            clone1_scnp = clone_specific_fcnp[clone1_name]
            clone2_scnp = clone_specific_fcnp[clone2_name]
            cn1a = force_non_negativity(int(row[CLONE1_CN_A]))
            cn1b = force_non_negativity(int(row[CLONE1_CN_B]))
            clone1_scnp.set_cn_record_for_segment(segment=fragment, cn=cn1a, haplotype=Haplotype.A)
            clone1_scnp.set_cn_record_for_segment(segment=fragment, cn=cn1b, haplotype=Haplotype.B)
            cn2a = get_subclonal_cn(subclonal_cn_string=row[CLONE2_CN_A], clonal_cn_int=cn1a)
            cn2b = get_subclonal_cn(subclonal_cn_string=row[CLONE2_CN_B], clonal_cn_int=cn1b)
            clone2_scnp.set_cn_record_for_segment(segment=fragment, cn=cn2a, haplotype=Haplotype.A)
            clone2_scnp.set_cn_record_for_segment(segment=fragment, cn=cn2b, haplotype=Haplotype.B)
            hapl_fragments[fragment.stable_id_non_hap] = fragment
    return hapl_fragments, clone_specific_fcnp


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Converter from Battenberg to DASCK segment copy number results")
    parser.add_argument('battenberg_scn_file')
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    parser.add_argument("--max-merge-gap", type=int, default=1000000)
    parser.add_argument("--no-fill-gaps", action="store_false", dest="fill_gaps")
    parser.add_argument("--max-fill-gap", type=int, default=1000000)
    parser.add_argument('dasck_scn_file')
    args = parser.parse_args()
    clone1_id = '1'
    clone2_id = '2'
    sample_name = args.sample
    scnt = {clone1_id: SegmentCopyNumberProfile(), clone2_id: SegmentCopyNumberProfile()}
    segments = []
    battenberg_scn_file_path = os.path.expanduser(args.battenberg_scn_file)
    battenberg_scn_file_path = os.path.abspath(battenberg_scn_file_path)
    with open(battenberg_scn_file_path, "rt") as source:
        reader = csv.DictReader(source)
        for row in reader:
            if row[SAMPLE_NAME] != sample_name:
                continue
            chr_name = row[CHROMOSOME]
            start_coordinate = int(row[START_POSITION])
            end_coordinate = int(row[END_POSITION])
            sp = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
            ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
            fragment = Segment(start_position=sp, end_position=ep)
            clone1_scnp = scnt[clone1_id]
            clone2_scnp = scnt[clone2_id]
            cn1a = force_non_negativity(cn=int(row[CLONE1_CN_A]))
            cn1b = force_non_negativity(cn=int(row[CLONE1_CN_B]))
            cn2a = get_subclonal_cn(subclonal_cn_string=row[CLONE2_CN_A], clonal_cn_int=cn1a)
            cn2b = get_subclonal_cn(subclonal_cn_string=row[CLONE2_CN_B], clonal_cn_int=cn1b)
            fid = fragment.stable_id_non_hap
            clone1_scnp.set_cn_record(sid=fid, hap=Haplotype.A, cn=cn1a)
            clone1_scnp.set_cn_record(sid=fid, hap=Haplotype.B, cn=cn1b)
            clone2_scnp.set_cn_record(sid=fid, hap=Haplotype.A, cn=cn2a)
            clone2_scnp.set_cn_record(sid=fid, hap=Haplotype.B, cn=cn2b)
            segments.append(fragment)
    if args.merge_fragments or args.fill_gaps:
        segments, scnt = refined_segments_and_tensor(segments=segments, tensor=scnt,
                                                     merge_fragments=args.merge_fragments, max_merge_gap=args.max_merge_gap,
                                                     fill_gaps=args.fill_gaps, max_fill_gap=args.max_fill_gap)
    write_scn_tensor(file_name=args.dasck_scn_file, scnt=scnt, segments=segments)
