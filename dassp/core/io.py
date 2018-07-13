# -*- coding: utf-8 -*-
import vcf

from dassp.core.structures import AdjacencyCopyNumberProfile
from dassp.core.structures import SegmentCopyNumberRecord, SegmentCopyNumberProfile, Haplotype, AdjacencyType, Phasing
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


def parse_cn_entry(cn_string, cn_separator=";"):
    result = dict()
    clone_specific_entries = cn_string.split(cn_separator)
    for entry in clone_specific_entries:
        entry = entry[1:-1]
        clone_id, cna, cnb = entry.split(",")
        cna, cnb = int(cna), int(cnb)
        result[clone_id] = {Haplotype.A: cna, Haplotype.B: cnb}
    return result


def parse_scn_record(scn_string, separator="\t", cn_separator=";", allow_unit_segments=False):
    data = scn_string.split(separator)
    chr_name = data[0]
    start_coordinate = int(data[1])
    end_coordinate = int(data[2])
    cn_entry = data[3]
    start_position = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
    end_position = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    segment = Segment(start_position=start_position, end_position=end_position, allow_unit_length=allow_unit_segments)
    cns = parse_cn_entry(cn_string=cn_entry, cn_separator=cn_separator)
    return segment, cns


def read_scn_tensor(file_name, clone_ids, separator="\t", cn_separator=";", allow_absent_clones=False, allow_unit_segments=False):
    segments = []
    result = {}
    for clone_id in clone_ids:
        result[clone_id] = SegmentCopyNumberProfile()
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            segment, cns = parse_scn_record(scn_string=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in clone_ids:
                if not allow_absent_clones and clone_id not in cns:
                    raise Exception()
            segments.append(segment)
            sid = segment.stable_id_non_hap
            for clone_id in clone_ids:
                scnp = result[clone_id]
                cna = cns[clone_id][Haplotype.A]
                cnb = cns[clone_id][Haplotype.B]
                scnp.set_cn_record(sid=sid, hap=Haplotype.A, cn=cna)
                scnp.set_cn_record(sid=sid, hap=Haplotype.B, cn=cnb)
    return segments, result


def get_cn_string_entry(clone_id, cna, cnb):
    return "({clone_id},{cna},{cnb})".format(clone_id=str(clone_id), cna=str(cna), cnb=str(cnb))


def get_scn_string_entry(segment, scnt, clone_ids, separator="\t", cn_separator=";"):
    chr_string = "{chr_name}".format(chr_name=str(segment.chromosome))
    start_coord_string = "{start_coord}".format(start_coord=str(segment.start_position.coordinate))
    end_coord_string = "{end_coord}".format(end_coord=str(segment.end_position.coordinate))
    cn_strings = []
    sid = segment.stable_id_non_hap
    for clone_id in clone_ids:
        scnp = scnt[clone_id]
        cna = scnp.get_cn(sid=sid, haplotype=Haplotype.A, default=0)
        cnb = scnp.get_cn(sid=sid, haplotype=Haplotype.B, default=0)
        cn_string = get_cn_string_entry(clone_id=clone_id, cna=cna, cnb=cnb)
        cn_strings.append(cn_string)
    cn_string = cn_separator.join(cn_strings)
    all_string_entries = [chr_string, start_coord_string, end_coord_string, cn_string]
    return separator.join(all_string_entries)


def write_scn_tensor(file_name, scnt, segments, separator="\t", cn_separator=";"):
    clone_ids = sorted(scnt.keys())
    with open(file_name, "wt") as destination:
        print("chr", "start", "end", "CNs", sep=separator, file=destination)
        segments = sorted(segments, key=lambda s: (s.chromosome, s.start_position.coordinate, s.end_position.coordinate))
        for segment in segments:
            scn_string = get_scn_string_entry(segment=segment, scnt=scnt, clone_ids=clone_ids, separator=separator, cn_separator=cn_separator)
            print(scn_string, file=destination)


EXTERNAL_NA_ID = "external_na_id"


def parse_na_string_entry(na_string, separator="\t"):
    data = na_string.split(separator)
    na_id = data[0]
    na_chr1_name = data[1]
    na_coord1 = int(data[2])
    na_strand1 = Strand.from_pm_string(string=data[3])
    na_chr2_name = data[4]
    na_coord2 = int(data[5])
    na_strand2 = Strand.from_pm_string(string=data[6])
    pos1 = Position(chromosome=na_chr1_name, coordinate=na_coord1, strand=na_strand1)
    pos2 = Position(chromosome=na_chr2_name, coordinate=na_coord2, strand=na_strand2)
    na = Adjacency(position1=pos1, position2=pos2, adjacency_type=AdjacencyType.NOVEL, extra={EXTERNAL_NA_ID: na_id})
    return na


def read_novel_adjacencies(file_name, separator="\t"):
    result = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            na = parse_na_string_entry(na_string=line, separator=separator)
            result.append(na)
    return result


def get_adjacency_string_entry(adjacency, separator="\t"):
    na_strings_entries = []
    na_strings_entries.append("{aid}".format(aid=adjacency.extra.get(EXTERNAL_NA_ID, "NA")))
    na_strings_entries.append("{chr1}".format(chr1=adjacency.position1.chromosome))
    na_strings_entries.append("{coord1}".format(coord1=adjacency.position1.coordinate))
    na_strings_entries.append("{strand1}".format(strand1=adjacency.position1.strand))
    na_strings_entries.append("{chr2}".format(chr2=adjacency.position2.chromosome))
    na_strings_entries.append("{coord2}".format(coord2=adjacency.position2.coordinate))
    na_strings_entries.append("{strand2}".format(strand2=adjacency.position2.strand))
    return separator.join(na_strings_entries)


def get_adjacency_cn_string_entry(adjacency, acnt, clone_id):
    result = [str(clone_id)]
    aid = adjacency.stable_id_non_phased
    for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
        result.append(str(acnt[clone_id].get_cn(aid=aid, phasing=ph, default=0)))
    return "(" + ",".join(result) + ")"


def get_acnt_string_entry(adjacency, acnt, clone_ids, separator="\t", cn_separator=";"):
    result = get_adjacency_string_entry(adjacency=adjacency, separator=separator)
    entries = []
    for clone_id in sorted(clone_ids):
        entries.append(get_adjacency_cn_string_entry(adjacency=adjacency, acnt=acnt, clone_id=clone_id))
    cns_entry = cn_separator.join(entries)
    return separator.join([result, cns_entry])


def write_adjacencies(file_name, adjacencies, separator="\t"):
    adjacencies = sorted(adjacencies, key=lambda a: (a.position1.chromosome, a.position1.coordinate, a.position2.chromosome, a.position2.coordinate))
    with open(file_name, "wt") as destination:
        print("aid", "chr1", "coord1", "strand1", "chr2", "coord2", "strand2", sep="\t", file=destination)
        for adjacency in adjacencies:
            na_string_entry = get_adjacency_string_entry(adjacency=adjacency, separator=separator)
            print(na_string_entry, file=destination)


def write_acnt(file_name, acnt, adjacencies, separator="\t", cn_separator=";", output_reference=True):
    clone_ids = sorted(acnt.keys())
    with open(file_name, "wt") as destination:
        ref_adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, adjacencies))
        nov_adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, adjacencies))
        print("aid", "chr1", "coord1", "strand1", "chr2", "coord2", "strand2", "CNs", file=destination, sep=separator)
        entries = nov_adjacencies
        if output_reference:
            entries += ref_adjacencies
        for adjacency in entries:
            na_string = get_acnt_string_entry(adjacency=adjacency, acnt=acnt, clone_ids=clone_ids, separator=separator, cn_separator=cn_separator)
            print(na_string, file=destination)


def parse_acn_string(string_entry, separator="\t", cn_separator=";"):
    cns = {}
    data = string_entry.split(separator)
    naid = data[0]
    chr1 = data[1]
    coord1 = int(data[2])
    strand1 = Strand.from_pm_string(string=data[3])
    chr2 = data[4]
    coord2 = int(data[5])
    strand2 = Strand.from_pm_string(string=data[6])
    cns_string_entry = data[7]
    p1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
    p2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
    adj_type = AdjacencyType.REFERENCE if "r" in naid else AdjacencyType.NOVEL
    adjacency = Adjacency(position1=p1, position2=p2, adjacency_type=adj_type, extra={EXTERNAL_NA_ID: naid})
    cn_data = cns_string_entry.split(cn_separator)
    for entry in cn_data:
        entry = entry[1:-1]
        clone_id, cnAA, cnAB, cnBA, cnBB = entry.split(",")
        cnAA, cnAB, cnBA, cnBB = int(cnAA), int(cnAB), int(cnBA), int(cnBB)
        cns[clone_id] = {Phasing.AA: cnAA, Phasing.AB: cnAB, Phasing.BA: cnBA, Phasing.BB: cnBB}
    return adjacency, cns


def read_acnt(file_name, clone_ids, separator="\t", cn_separator=";"):
    clone_ids = sorted(clone_ids)
    result = {clone_id: AdjacencyCopyNumberProfile() for clone_id in clone_ids}
    adjacencies = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            adjacency, cns = parse_acn_string(string_entry=line, separator=separator, cn_separator=cn_separator)
            if EXTERNAL_NA_ID in adjacency.extra and adjacency.extra[EXTERNAL_NA_ID].startswith("r"):
                adjacency.adjacency_type = AdjacencyType.REFERENCE
            else:
                adjacency.adjacency_type = AdjacencyType.NOVEL
            aid = adjacency.stable_id_non_phased
            adjacencies.append(adjacency)
            for clone_id in clone_ids:
                if clone_id not in cns:
                    raise Exception()
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    acnp = result[clone_id]
                    acnp._set_cn_record(aid=aid, phasing=ph, cn=cns[clone_id][ph])
    return adjacencies, result


def get_all_clone_ids_from_dasck_scnt_file(file_name, separator="\t", cn_separator=";", allow_unit_segments=False):
    clone_ids = set()
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            segment, cns = parse_scn_record(scn_string=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in cns.keys():
                clone_ids.add(clone_id)
    return sorted(clone_ids)


def parser_position_string_entry(position_string, separator):
    data = position_string.split(separator)
    chr_name = data[0]
    coordinate = int(data[1])
    strand = Strand.from_pm_string(string=data[2])
    return Position(chromosome=chr_name, coordinate=coordinate, strand=strand)


def read_positions(file_name, separator="\t"):
    result = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if len(line) == 0 or line.startswith("#") or line_cnt == 0:
                continue
            position = parser_position_string_entry(position_string=line, separator=separator)
            result.append(position)
    return result


def get_segment_string_entry(segment, separator="\t"):
    chr_name_entry = "{chr_name}".format(chr_name=str(segment.start_position.chromosome))
    start_entry = "{start_coord}".format(start_coord=str(segment.start_position.coordinate))
    end_entry = "{end_coord}".format(end_coord=str(segment.end_position.coordinate))
    return separator.join([chr_name_entry, start_entry, end_entry])


def write_segments(file_name, segments, separator="\t"):
    with open(file_name, "wt") as destination:
        segments = sorted(segments, key=lambda s: (s.start_position.chromosome, s.start_position.coordinate, s.end_position.coordinate))
        print("chr", "start", "end", sep=separator, file=destination)
        for segment in segments:
            string_entry = get_segment_string_entry(segment=segment, separator=separator)
            print(string_entry, file=destination)


def parse_segment_string_entry(string_entry, separator="\t", allow_unit_segments=False):
    data = string_entry.split(separator)
    chr_name = data[0]
    start_corodinate = int(data[1])
    end_coordinate = int(data[2])
    sp = Position(chromosome=chr_name, coordinate=start_corodinate, strand=Strand.REVERSE)
    ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    return Segment(start_position=sp, end_position=ep, allow_unit_length=allow_unit_segments)


def read_segments(file_name, separator="\t", allow_unit_segments=False):
    result = []
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            segment = parse_segment_string_entry(string_entry=line, separator=separator, allow_unit_segments=allow_unit_segments)
            result.append(segment)
    return result


def get_cn_boundaries_string_entry(segment, scn_boundaries, clone_ids, cn_separator=";"):
    result = []
    sid = segment.stable_id_non_hap
    for clone_id in clone_ids:
        lower_a, upper_a = scn_boundaries[clone_id][sid][Haplotype.A]
        lower_b, upper_b = scn_boundaries[clone_id][sid][Haplotype.B]
        string_a = "{lower}-{upper}".format(lower=lower_a, upper=upper_a)
        string_b = "{lower}-{upper}".format(lower=lower_b, upper=upper_b)
        string_entry = "({clone_id},{boundaries_a},{boundaries_b})".format(clone_id=str(clone_id), boundaries_a=str(string_a), boundaries_b=str(string_b))
        result.append(string_entry)
    return cn_separator.join(result)


def write_scn_boundaries(file_name, segments, scn_boundaries, separator="\t", cn_separator=";"):
    clone_ids = sorted(scn_boundaries.keys())
    with open(file_name, "wt") as destination:
        print("chr", "start", "end", "CNBoundaries", sep=separator, file=destination)
        for s in segments:
            string_entries = []
            string_entries.append("{chr_name}".format(chr_name=str(s.start_position.chromosome)))
            string_entries.append("{start}".format(start=str(s.start_position.coordinate)))
            string_entries.append("{end}".format(end=str(s.end_position.coordinate)))
            string_entries.append(get_cn_boundaries_string_entry(segment=s, scn_boundaries=scn_boundaries, clone_ids=clone_ids, cn_separator=cn_separator))
            string_entry = separator.join(string_entries)
            print(string_entry, file=destination)


def parse_cn_boundaries_entry_string(string_entry):
    data = string_entry.split("-")
    return int(data[0]), int(data[1])


def parse_cn_boundaries_string(string_entry, cn_separator=";"):
    result = {}
    data = string_entry.split(cn_separator)
    for entry in data:
        entry = entry[1:-1]
        entry_data = entry.split(",")
        clone_id = entry_data[0]
        boundaries_a_string = entry_data[1]
        boundaries_b_string = entry_data[2]
        boundaries_a = parse_cn_boundaries_entry_string(string_entry=boundaries_a_string)
        boundaries_b = parse_cn_boundaries_entry_string(string_entry=boundaries_b_string)
        result[clone_id] = {Haplotype.A: boundaries_a, Haplotype.B: boundaries_b}
    return result


def parse_scn_boundaries_string_entry(string_entry, separator="\t", cn_separator=";", allow_unit_segments=False):
    data = string_entry.split(separator)
    chr_name = data[0]
    start_coordinate = int(data[1])
    end_coordinate = int(data[2])
    cn_boundaries_string = data[3]
    sp = Position(chromosome=chr_name, coordinate=start_coordinate, strand=Strand.REVERSE)
    ep = Position(chromosome=chr_name, coordinate=end_coordinate, strand=Strand.FORWARD)
    segment = Segment(start_position=sp, end_position=ep, allow_unit_length=allow_unit_segments)
    boundaries = parse_cn_boundaries_string(string_entry=cn_boundaries_string, cn_separator=cn_separator)
    return segment, boundaries


def read_scn_boundaries(file_name, segments, clone_ids, separator="\t", cn_separator=";", allow_unit_segments=False):
    result = {clone_id: {s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in segments} for clone_id in clone_ids}
    with open(file_name, "rt") as source:
        for line_cnt, line in enumerate(source):
            line = line.strip()
            if line_cnt == 0 or len(line) == 0 or line.startswith("#"):
                continue
            segment, boundaries = parse_scn_boundaries_string_entry(string_entry=line, separator=separator, cn_separator=cn_separator, allow_unit_segments=allow_unit_segments)
            for clone_id in clone_ids:
                result[clone_id][segment.stable_id_non_hap][Haplotype.A] = boundaries[clone_id][Haplotype.A]
                result[clone_id][segment.stable_id_non_hap][Haplotype.B] = boundaries[clone_id][Haplotype.B]
    return result
