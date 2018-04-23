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
from dassp.core.structures import SegmentCopyNumberProfile, Position, Strand, Segment, Haplotype

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Converter from ReMixT to DASCK segment copy number results")
    parser.add_argument('remixt_scn_file')
    parser.add_argument('dasck_scn_file')
    args = parser.parse_args()
    clone1_id = '1'
    clone2_id = '2'
    scnt = {clone1_id: SegmentCopyNumberProfile(), clone2_id: SegmentCopyNumberProfile()}
    segments = []
    remixt_scn_file_path = os.path.expanduser(args.remixt_scn_file)
    with open(remixt_scn_file_path, "rt") as source:
        reader = csv.DictReader(source, delimiter='\t')
        for line in reader:
            chr_name = line['chromosome']
            start_coordinate = int(line['start']) - 1
            end_coordinate = int(line['end'])
            sp = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
            ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
            segment = Segment(start_position=sp, end_position=ep)
            segments.append(segment)
            sid = segment.stable_id_non_hap
            clone_1_cn_a = int(line['major_1'])
            clone_1_cn_b = int(line['minor_1'])
            clone_2_cn_a = int(line['major_2'])
            clone_2_cn_b = int(line['minor_2'])
            clone1_scnp = scnt[clone1_id]
            clone2_scnp = scnt[clone2_id]
            clone1_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=clone_1_cn_a)
            clone1_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=clone_1_cn_b)
            clone2_scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=clone_2_cn_a)
            clone2_scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=clone_2_cn_b)
    write_scn_tensor(file_name=args.dasck_scn_file, scnt=scnt, segments=segments)
