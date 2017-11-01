# -*- coding: utf-8 -*-
import vcf

from core.structures import SegmentCopyNumberRecord
from dassp.core.structures import Position, Strand, Adjacency, Segment


class VCFNovelAdjacencyReader(object):
    pass


class PCAWGVCFNovelAdjacencyReader(VCFNovelAdjacencyReader):
    def __init__(self, file_path):
        self.file_path = file_path
        self.reader = vcf.Reader(filename=file_path)

    @classmethod
    def get_pcawg_adjacency_id_from_position_id(cls, position_id, full=False):
        if full:
            return position_id[:-2]
        return position_id[7:-2]

    @classmethod
    def parse_position_record(cls, position_record):
        extra = {
            "self_id": position_record.ID,
            "mate_id": position_record.INFO["MATEID"],
            "sv_class": position_record.INFO["SVCLASS"],
            "sv_type": position_record.INFO["SVTYPE"],
            "mate_chrom": position_record.INFO["MATECHROM"]
        }
        return Position(chromosome=position_record.CHROM,
                        coordinate=position_record.POS,
                        strand=Strand.from_pm_string(position_record.INFO["STRAND"]),
                        extra=extra)

    def parse(self):
        na_positions = self.parse_pcawg_vcf_position_file()
        n_adjacencies = self.get_adjacencies_from_positions(na_positions)
        return na_positions, n_adjacencies

    def parse_pcawg_vcf_position_file(self):
        positions = []
        for record in vcf.Reader(filename=self.file_path):
            position = self.parse_position_record(position_record=record)
            positions.append(position)
        return positions

    @classmethod
    def get_adjacencies_from_positions(cls, na_positions):
        na_positions_by_ids = {}
        for position in na_positions:
            na_positions_by_ids[position.extra["self_id"]] = position
        n_adjacencies = []
        processed_positions = set()
        for position in na_positions:
            mate_position = na_positions_by_ids[position.extra["mate_id"]]
            if (position.extra["self_id"] in processed_positions and mate_position.extra["self_id"] not in processed_positions) or \
                    (position.extra["self_id"] not in processed_positions and mate_position.extra["self_id"] in processed_positions):
                raise ValueError("Mate positions {p1id} and {p2id} have not been processed together (i.e., one was processed without the other). "
                                 "May indicate a problem with the position ids".format(p1id=position.extra["self_id"], p2id=mate_position.extra["self_id"]))
            if position.extra["self_id"] in processed_positions and mate_position.extra["self_id"] in processed_positions:
                continue
            if position.extra["self_id"] != mate_position.extra["mate_id"] or position.extra["mate_id"] != mate_position.extra["self_id"]:
                raise ValueError("SELFID<->MATEID relationship is not symmetric for positions {p1id} and {p2id}".format(p1id=position.extra["self_id"],
                                                                                                                        p2id=mate_position.extra["self_id"]))
            adjacency_id = cls.get_pcawg_adjacency_id_from_position_id(position_id=position.extra["self_id"])
            if adjacency_id != cls.get_pcawg_adjacency_id_from_position_id(position_id=mate_position.extra["self_id"]):
                raise ValueError("Adjacency ID (i.e., similar parts of positoins ids for mated positions) is not the same when extracted from mated positions")
            adjacency = Adjacency(position1=position, position2=mate_position, idx=adjacency_id)
            n_adjacencies.append(adjacency)
            processed_positions.add(position.extra["self_id"])
            processed_positions.add(mate_position.extra["self_id"])
        return n_adjacencies


class BattenbergSegmentReader(object):
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_battenberg_file(self):
        records = []
        with open(self.file_path, "rt") as source:
            for cnt, line in enumerate(source):
                if cnt == 0:
                    continue
                record = self.parse(battenberg_record=line)
                records.append(record)
        return records

    @classmethod
    def parse(cls, battenberg_record, separator="\t"):
        data = battenberg_record.split(separator)
        chromosome = data[1]
        start_position = Position(chromosome=chromosome, coordinate=int(data[2]), strand=Strand.REVERSE)
        end_position = Position(chromosome=chromosome, coordinate=int(data[3]), strand=Strand.FORWARD)
        maj_a_cn = int(data[8])
        maj_b_cn = int(data[9])
        min_a_cn = int(data[11]) if data[11] != "NA" else None
        min_b_cn = int(data[12]) if data[12] != "NA" else None
        maj_frac = float(data[10])
        min_frac = float(data[13]) if data[13] != "NA" else float(0)
        extra = {"battenberg_id": data[0], "maj_frac": maj_frac, "min_frac": min_frac}
        segment = Segment(start_position=start_position, end_position=end_position)
        return SegmentCopyNumberRecord(segment=segment, maj_a_cn=maj_a_cn, maj_b_cn=maj_b_cn, min_a_cn=min_a_cn, min_b_cn=min_b_cn, extra=extra)
