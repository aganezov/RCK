import csv
import re
from collections import defaultdict
from copy import deepcopy
from enum import Enum

import vcf

from rck.core.io import EXTERNAL_NA_ID, SVTYPE, write_adjacencies_to_destination, COPY_NUMBER, parse_segment_extra
from rck.core.structures import Strand, Position, Adjacency, AdjacencyType, Phasing, Segment

VCF_FILTER = "vcf_filter"
VCF_FILTER2 = "vcf_filter2"
STRIP_CHR = "strip_CHR"
APPEND_ID = "append_id"
EXTRA_FIELDS = "extra_fields"
ID_SUFFIX = "id_suffix"
SVTYPE2 = "svtype2"
EXTERNAL_NA_ID2 = "external_na_id2"
SVLEN = "svlen"
SUPPORT_READ_NAMES = "support_read_names"


def get_most_supported_strands(strands_list):
    result_by_support = {}
    for entry in strands_list:
        strands, support = entry.split(":")
        support = int(support)
        result_by_support[support] = strands
    max_support = max(result_by_support.keys())
    return result_by_support[max_support]


def strip_chr(chr_string):
    if chr_string.startswith("chr"):
        return chr_string[3:]
    return chr_string


def get_strand_string_from_alt_breakend_string(alt_breakend_string):
    if alt_breakend_string.startswith("]") or alt_breakend_string.endswith("]"):
        return "+"
    return "-"


def get_strand_string_from_alt_breakend(alt_breakend):
    return get_strand_string_from_alt_breakend_string(alt_breakend_string=str(alt_breakend))


def update_dict_with_vcf_info(target, vcf_record):
    for key, value in vcf_record:
        target[str(key).lower()] = value
    return target


def get_inv_chrs_from_vcf_record(vcf_record):
    chr11 = str(vcf_record.CHROM)
    chr12 = chr11
    chr21 = chr11
    chr22 = chr11
    return chr11, chr12, chr21, chr22


def get_inv_na1_chrs_from_vcf_record(vcf_record):
    return str(vcf_record.CHROM), str(vcf_record.CHROM)


def get_inv_na2_chrs_from_vcf_record(vcf_record):
    return get_inv_na1_chrs_from_vcf_record(vcf_record=vcf_record)


def get_inv_na1_coords_from_vcf_record(vcf_record):
    return int(vcf_record.POS) - 1, int(vcf_record.INFO["END"])


def get_inv_na2_coords_from_vcf_record(vcf_record):
    a, b = get_inv_na1_coords_from_vcf_record(vcf_record=vcf_record)
    return a + 1, b + 1


def get_inv_na1_strands_from_vcf_record(vcf_record):
    return Strand.FORWARD, Strand.FORWARD


def get_inv_na2_strands_from_vcf_record(vcf_record):
    return Strand.REVERSE, Strand.REVERSE


def get_del_na_chrs_from_vcf_record(vcf_record):
    return get_inv_na1_chrs_from_vcf_record(vcf_record=vcf_record)


def get_del_na_coords_from_vcf_record(vcf_record):
    return int(vcf_record.POS) - 1, int(vcf_record.INFO["END"]) + 1


def get_del_na_strands_from_vcf_record(vcf_record):
    return Strand.FORWARD, Strand.REVERSE


def get_dup_na_chrs_from_vcf_record(vcf_record):
    return get_inv_na1_chrs_from_vcf_record(vcf_record=vcf_record)


def get_dup_na_coords_from_vcf_record(vcf_record):
    return int(vcf_record.POS), int(vcf_record.INFO["END"])


def get_dup_na_strands_from_vcf_record(vcf_record):
    return Strand.REVERSE, Strand.FORWARD


def get_strand_from_alt_breakend(vcf_recrod):
    alt_string = str(vcf_recrod.ALT[0])
    if alt_string.startswith("]") or alt_string.endswith("]"):
        return Strand.FORWARD
    return Strand.REVERSE


def get_vcf_records_from_file(vcf_file_name):
    with open(vcf_file_name, "rt") as source:
        return get_vcf_records_from_source(source=source)


def get_vcf_records_from_source(source):
    result = []
    vcf_reader = vcf.Reader(source)
    if "RNAMES" in vcf_reader.infos:
        vcf_reader.infos["RNAMES"] = vcf_reader.infos["RNAMES"]._replace(num='.')
    for record in vcf_reader:
        result.append(record)
    return result


def get_vcf_records_by_ids(vcf_records):
    return {r.ID: r for r in vcf_records}


def get_nas_from_lumpy_vcf_file(lumpy_vcf_file):
    vcf_records = get_vcf_records_from_file(vcf_file_name=lumpy_vcf_file)
    return get_nas_from_lumpy_vcf_records(lumpy_vcf_records=vcf_records)


def update_nas_ids(nas_by_ids_defaultdict, setup):
    result = {}
    for sv_id, nas in nas_by_ids_defaultdict.items():
        if len(nas) > 1:
            for cnt, na in enumerate(nas, start=1):
                na.extra[EXTERNAL_NA_ID] += "_{cnt}".format(cnt=cnt)
        if setup.get(APPEND_ID, True):
            for na in nas:
                if not na.extra[EXTERNAL_NA_ID].endswith(setup.get(ID_SUFFIX, "")):
                    na.extra[EXTERNAL_NA_ID] += "_{suffix}".format(suffix=setup.get(ID_SUFFIX, ""))
        for na in nas:
            na.extra[EXTERNAL_NA_ID] = na.extra[EXTERNAL_NA_ID].replace(":", "_")
            result[na.extra[EXTERNAL_NA_ID]] = na
    return result


PASS_FILTER_STR = "PASS"


def get_string_vcf_filter(record):
    result = ",".join(map(str, record.FILTER)) if isinstance(record.FILTER, list) else str(record.FILTER)
    if len(result) == 0:
        result = PASS_FILTER_STR
    return result


class StandardizedSVType(Enum):
    INS = "INS"
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    TRA = "TRA"


def get_standardize_sv_type(adjacency: Adjacency):
    # generic for all callers
    if adjacency.position1.chromosome != adjacency.position2.chromosome:
        return StandardizedSVType.TRA
    strands = adjacency.position1.strand, adjacency.position2.strand
    # generic for all callers
    if strands == (Strand.REVERSE, Strand.FORWARD):
        return StandardizedSVType.DUP
    # generic for all callers
    if strands[0] == strands[1]:
        return StandardizedSVType.INV
    lower_extra = {key.lower(): value for key, value in adjacency.extra.items()}
    if "or_svtype" in lower_extra:
        ###
        # covers Sniffles, PBSV, Delly, ...
        ###
        if "ins" in lower_extra["or_svtype"].lower():
            return StandardizedSVType.INS
        ###
        # covers Sniffles, PBSV, Delly, ...
        ###
        if "del" in lower_extra["or_svtype"].lower():
            return StandardizedSVType.DEL
    ###
    # more generic approach if both REF and ALT fields are present in the adjacency data
    ###
    if "alt" in lower_extra and "ref" in lower_extra and all(map(lambda e: str(e).isalpha() and str(e) != "None", lower_extra["alt"])) and lower_extra["ref"].isalpha():
        if isinstance(lower_extra["alt"], list):
            alt = lower_extra["alt"][0]
        else:
            alt = lower_extra["alt"]
        if len(alt) > len(lower_extra["ref"]):
            return StandardizedSVType.INS
        else:
            return StandardizedSVType.DEL
    ###
    # last resort based on positive/negative svlen
    ###
    if SVLEN.lower() in lower_extra:
        if isinstance(lower_extra[SVLEN.lower()], list):
            length = int(float(lower_extra[SVLEN.lower()][0]))
        else:
            length = int(float(lower_extra[SVLEN.lower()]))
        if length > 0:
            return StandardizedSVType.INS
    return StandardizedSVType.DEL


def update_adjacencies_svtype(adjacencies):
    for adj in adjacencies:
        adj_lower_extra = {key.lower(): value for key, value in adj.extra.items()}
        or_key = "or_svtype"
        if SVTYPE.lower() in adj_lower_extra:
            adj_lower_extra[or_key] = adj_lower_extra["svtype"]
        adj.extra = adj_lower_extra
        adj.extra[SVTYPE.lower()] = get_standardize_sv_type(adjacency=adj).value


DUPLICATED_ENTRIES_EXTRA = {
    "CHR2",
    "END",
    "STRANDS"
}


def clear_duplicated_entries_extra(extra):
    to_remove = []
    for key, value in extra.items():
        if key.upper() in DUPLICATED_ENTRIES_EXTRA:
            to_remove.append(key)
    for key in to_remove:
        extra.pop(key)
    return extra


def get_nas_from_lumpy_vcf_records(lumpy_vcf_records, setup=None):
    if setup is None:
        setup = {}
    records_by_ids = get_vcf_records_by_ids(vcf_records=lumpy_vcf_records)
    nas_by_ids = defaultdict(list)
    processed_vcf_ids = set()
    for record in records_by_ids.values():
        base_sv_id = str(record.ID).split("_")[0]
        if base_sv_id in processed_vcf_ids:
            continue
        sv_type = record.INFO["SVTYPE"]
        extra = {}
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        if sv_type in ["DUP", "DEL", "BND"]:
            if sv_type in ["DUP", "DEL"]:
                if "_" in str(record.ID):
                    raise ValueError("Non standard id {sv_id} for a same-chromosome SV of type {sv_type}".format(sv_id=record.ID, sv_type=sv_type))
                extra.update(record.INFO)
                chr1, chr2 = get_dup_na_chrs_from_vcf_record(vcf_record=record)
                if sv_type == "DUP":
                    strand1, strand2 = get_dup_na_strands_from_vcf_record(vcf_record=record)
                    coord1, coord2 = get_dup_na_coords_from_vcf_record(vcf_record=record)
                else:
                    strand1, strand2 = get_del_na_strands_from_vcf_record(vcf_record=record)
                    coord1, coord2 = get_del_na_coords_from_vcf_record(vcf_record=record)
                add_record_ref_alt_to_extra(sv_type, record, extra)
            else:
                if "_" not in str(record.ID):
                    raise ValueError("Non standard id {sv_id} for a BND SV".format(sv_id=record.ID))
                extra.update(record.INFO)
                chr1 = str(record.CHROM)
                coord1 = int(record.POS)
                mate_vcf_entry_id_list = record.INFO["MATEID"]
                mate_vcf_entry_id = mate_vcf_entry_id_list[0]
                mate_record = records_by_ids[mate_vcf_entry_id]
                extra.update(mate_record.INFO)
                extra[VCF_FILTER2] = get_string_vcf_filter(record=mate_record)
                strand1 = get_strand_from_alt_breakend(vcf_recrod=mate_record)
                strand2 = get_strand_from_alt_breakend(vcf_recrod=record)
                chr2 = str(mate_record.CHROM)
                coord2 = int(mate_record.POS)
            if setup.get(STRIP_CHR, True):
                chr1 = strip_chr(chr_string=chr1)
                chr2 = strip_chr(chr_string=chr2)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = base_sv_id
            extra = clear_duplicated_entries_extra(extra=extra)
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            nas_by_ids[base_sv_id].append(na)
            processed_vcf_ids.add(base_sv_id)
        elif sv_type == "INV":
            if "_" in str(record.ID):
                raise Exception("Non standard id {sv_id} for a same-chromosome SV of type {sv_type}".format(sv_id=record.ID, sv_type=sv_type))
            extra1 = deepcopy(record.INFO)
            extra1[VCF_FILTER] = get_string_vcf_filter(record=record)
            extra2 = deepcopy(record.INFO)
            extra2[VCF_FILTER] = get_string_vcf_filter(record=record)
            chr11, chr12 = get_inv_na1_chrs_from_vcf_record(vcf_record=record)
            chr21, chr22 = get_inv_na2_chrs_from_vcf_record(vcf_record=record)
            coord11, coord12 = get_inv_na1_coords_from_vcf_record(vcf_record=record)
            coord21, coord22 = get_inv_na2_coords_from_vcf_record(vcf_record=record)
            strand11, strand12 = get_inv_na1_strands_from_vcf_record(vcf_record=record)
            strand21, strand22 = get_inv_na2_strands_from_vcf_record(vcf_record=record)
            if setup.get(STRIP_CHR, True):
                chr11 = strip_chr(chr_string=chr11)
                chr12 = strip_chr(chr_string=chr12)
                chr21 = strip_chr(chr_string=chr21)
                chr22 = strip_chr(chr_string=chr22)
            pos11 = Position(chromosome=chr11, coordinate=coord11, strand=strand11)
            pos12 = Position(chromosome=chr12, coordinate=coord12, strand=strand12)
            pos21 = Position(chromosome=chr21, coordinate=coord21, strand=strand21)
            pos22 = Position(chromosome=chr22, coordinate=coord22, strand=strand22)
            extra1[EXTERNAL_NA_ID] = base_sv_id
            extra2[EXTERNAL_NA_ID] = base_sv_id
            extra1 = clear_duplicated_entries_extra(extra=extra1)
            extra2 = clear_duplicated_entries_extra(extra=extra2)
            na1 = Adjacency(position1=pos11, position2=pos12, extra=extra1)
            na2 = Adjacency(position1=pos21, position2=pos22, extra=extra2)
            nas_by_ids[base_sv_id].append(na1)
            nas_by_ids[base_sv_id].append(na2)
            processed_vcf_ids.add(base_sv_id)
        else:
            raise Exception("Unknown SVTYPE {sv_type}".format(sv_type=sv_type))
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return list(nas_by_ids.values())


def get_nas_from_longranger_vcf_records(longranger_vcf_records, setup=None):
    if setup is None:
        setup = {}
    records_by_ids = get_vcf_records_by_ids(vcf_records=longranger_vcf_records)
    nas_by_ids = defaultdict(list)
    processed_vcf_ids = set()
    for record in records_by_ids.values():
        sv_id_entries = str(record.ID).split("_")
        sv_id = str(record.ID)
        if sv_id in processed_vcf_ids:
            continue
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        sv_type = record.INFO["SVTYPE"]
        if len(sv_id_entries) <= 2:
            mate_present = False
        else:
            mate_present = True
        assert mate_present == ("MATEID" in record.INFO)
        if mate_present:
            base_sv_id = "_".join(sv_id_entries[:-1])
            for mate_vcf_id in record.INFO["MATEID"]:
                mate_record = records_by_ids[mate_vcf_id]
                extra.update(deepcopy(mate_record.INFO))
                extra[VCF_FILTER2] = get_string_vcf_filter(record=mate_record)
                mate_breakend = record.ALT[0]
                record_breakend = mate_record.ALT[0]
                chr1 = str(record.CHROM)
                chr2 = str(mate_record.CHROM)
                strand1 = Strand.from_pm_string(get_strand_string_from_alt_breakend(record_breakend))
                strand2 = Strand.from_pm_string(get_strand_string_from_alt_breakend(mate_breakend))
                coord1 = int(record.POS)
                coord2 = int(mate_record.POS)
                if setup.get(STRIP_CHR, True):
                    chr1 = strip_chr(chr_string=chr1)
                    chr2 = strip_chr(chr_string=chr2)
                pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
                pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
                extra[EXTERNAL_NA_ID] = base_sv_id
                na = Adjacency(position1=pos1, position2=pos2, extra=extra)
                nas_by_ids[base_sv_id].append(na)
                processed_vcf_ids.add(mate_vcf_id)
            processed_vcf_ids.add(sv_id)
        else:
            if sv_type == "UNK":
                continue
            base_sv_id = "_".join(sv_id_entries)
            if sv_type == "INV":
                chr11, chr12 = get_inv_na1_chrs_from_vcf_record(vcf_record=record)
                chr21, chr22 = get_inv_na2_chrs_from_vcf_record(vcf_record=record)
                coord11, coord12 = get_inv_na1_coords_from_vcf_record(vcf_record=record)
                coord21, coord22 = get_inv_na2_coords_from_vcf_record(vcf_record=record)
                strand11, strand12 = get_inv_na1_strands_from_vcf_record(vcf_record=record)
                strand21, strand22 = get_inv_na2_strands_from_vcf_record(vcf_record=record)
                if setup.get(STRIP_CHR, True):
                    chr11 = strip_chr(chr_string=chr11)
                    chr12 = strip_chr(chr_string=chr12)
                    chr21 = strip_chr(chr_string=chr21)
                    chr22 = strip_chr(chr_string=chr22)
                pos11 = Position(chromosome=chr11, coordinate=coord11, strand=strand11)
                pos12 = Position(chromosome=chr12, coordinate=coord12, strand=strand12)
                pos21 = Position(chromosome=chr21, coordinate=coord21, strand=strand21)
                pos22 = Position(chromosome=chr22, coordinate=coord22, strand=strand22)
                extra[EXTERNAL_NA_ID] = base_sv_id
                extra1 = deepcopy(extra)
                extra2 = deepcopy(extra)
                na1 = Adjacency(position1=pos11, position2=pos12, extra=extra1)
                na2 = Adjacency(position1=pos21, position2=pos22, extra=extra2)
                nas_by_ids[base_sv_id].append(na1)
                nas_by_ids[base_sv_id].append(na2)
                processed_vcf_ids.add(base_sv_id)
            elif sv_type in ["DEL", "DUP"]:
                chr1 = str(record.CHROM)
                if setup.get(STRIP_CHR, True):
                    chr1 = strip_chr(chr_string=chr1)
                chr2 = chr1
                coord1 = int(record.POS)
                coord2 = int(record.INFO["END"])
                if sv_type == "DUP":
                    strand1 = Strand.REVERSE
                    strand2 = Strand.FORWARD
                else:
                    strand1 = Strand.FORWARD
                    strand2 = Strand.REVERSE
                extra[EXTERNAL_NA_ID] = base_sv_id
                pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
                pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
                add_record_ref_alt_to_extra(sv_type, record, extra)
                na = Adjacency(position1=pos1, position2=pos2, extra=extra)
                nas_by_ids[base_sv_id].append(na)
                processed_vcf_ids.add(base_sv_id)
            else:
                raise Exception("Unknown SV type ({svtype}) for longranger".format(svtype=sv_type))
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    return list(nas_by_ids.values())


def add_record_ref_alt_to_extra(sv_type, record, extra):
    if any(stand_svtype.lower() in sv_type.lower() for stand_svtype in ["ins", "del"]) and hasattr(record.ALT[0], "sequence"):
        extra["ALT"] = [alt.sequence for alt in record.ALT]
        if len(record.REF) > 0:
            extra["REF"] = record.REF


def get_nas_from_manta_vcf_records(manta_vcf_records, setup=None):
    nas_by_ids = defaultdict(list)
    if setup is None:
        setup = {}
    records_by_ids = get_vcf_records_by_ids(vcf_records=manta_vcf_records)
    processed_ids = set()
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        svtype = record.INFO["SVTYPE"]
        record_id = str(record.ID)
        if record_id in processed_ids:
            continue
        if svtype == "BND":
            mate_vcf_id = record.INFO["MATEID"][0]
            if mate_vcf_id not in records_by_ids:
                continue
            mate_record = records_by_ids[mate_vcf_id]
            chr1 = str(record.CHROM)
            chr2 = str(mate_record.CHROM)
            coord1 = int(record.POS)
            coord2 = int(mate_record.POS)
            strand1 = Strand.from_pm_string(string=get_strand_string_from_alt_breakend(alt_breakend=mate_record.ALT))
            strand2 = Strand.from_pm_string(string=get_strand_string_from_alt_breakend(alt_breakend=record.ALT))
            if setup.get(STRIP_CHR, True):
                chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = str(record.ID)
            extra[VCF_FILTER2] = get_string_vcf_filter(record=mate_record)
            extra[EXTERNAL_NA_ID2] = mate_vcf_id
            extra = clear_duplicated_entries_extra(extra=extra)
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            nas_by_ids[na.extra[EXTERNAL_NA_ID]].append(na)
            processed_ids.add(record_id)
            processed_ids.add(mate_vcf_id)
        elif svtype in ["INV", "DUP", "DEL", "INS"]:
            chr1 = str(record.CHROM)
            if setup.get(STRIP_CHR, True):
                chr1 = strip_chr(chr_string=chr1)
            chr2 = chr1
            coord1 = int(record.POS)
            coord2 = int(record.INFO["END"])
            if svtype == "INV":
                if "INV3" in record.INFO:
                    strand1 = Strand.REVERSE
                    strand2 = Strand.REVERSE
                elif "INV5" in record.INFO:
                    strand1 = Strand.FORWARD
                    strand2 = Strand.FORWARD
                else:
                    raise Exception("Unknown INV type for Manta (i.e., not INV3, nor INV5)")
            elif svtype == "DUP":
                strand1 = Strand.REVERSE
                strand2 = Strand.FORWARD
            elif svtype == "DEL":
                strand1 = Strand.FORWARD
                strand2 = Strand.REVERSE
            else:
                # insertion
                strand1 = Strand.FORWARD
                strand2 = Strand.REVERSE
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = record_id
            add_record_ref_alt_to_extra(svtype, record, extra)
            if svtype == "INS":
                if "SVLEN" in record.INFO:
                    extra[SVLEN] = int(record.INFO["SVLEN"][0])
                    del extra["SVLEN"]
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            nas_by_ids[record_id].append(na)
            processed_ids.add(record_id)
        else:
            raise Exception("Unknown SV type {svtype}".format(svtype=svtype))
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return nas_by_ids.values()


def get_nas_from_sniffles_vcf_records(sniffles_vcf_records, setup=None):
    if setup is None:
        setup = {}
    nas_by_ids = defaultdict(list)
    records_by_ids = get_vcf_records_by_ids(vcf_records=sniffles_vcf_records)
    processed_ids = set()
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        for vcf_sample in record.samples:
            extra["OR_GT"] = "/".join(map(str, vcf_sample.gt_alleles)).replace("None", ".")
        record_id = str(record.ID)
        if record_id in processed_ids:
            continue
        svtype = record.INFO.get("SVTYPE", "")
        chr1 = record.CHROM
        coord1 = int(record.POS)
        if "CHR2" in record.INFO:
            chr2 = record.INFO["CHR2"]
            coord2 = int(record.INFO["END"]) if "END" in record.INFO else None
            if coord2 is None:
                if "SVLEN" in extra:
                    svlen = int(extra["SVLEN"])
                elif len(record.REF) > 0 and hasattr(record.ALT[0], "sequence"):
                    ref_len = len(record.REF)
                    alt_len = len(record.ALT[0].sequence)
                    if ref_len == 1:
                        svlen = alt_len
                    else:
                        svlen = ref_len
                else:
                    raise ValueError(f"Error while processing SV with id {record_id}. No indication of SVLEN nor LAT/REF sequences")
                coord2 = coord1 + abs(svlen)
        elif "ins" in svtype.lower() or "del" in svtype.lower():
            chr2 = chr1
            coord2 = int(record.INFO["END"]) if "END" in record.INFO else None
            if coord2 is None:
                if "SVLEN" in extra:
                    svlen = int(extra["SVLEN"])
                elif len(record.REF) > 0 and hasattr(record.ALT[0], "sequence"):
                    ref_len = len(record.REF)
                    alt_len = len(record.ALT[0].sequence)
                    if ref_len == 1:
                        svlen = alt_len
                    else:
                        svlen = ref_len
                else:
                    raise ValueError(f"Error while processing SV with id {record_id}. No indication of SVLEN nor LAT/REF sequences")
                coord2 = coord1 + abs(svlen)
        else:
            record_breakend = record.ALT[0]
            chr2 = record_breakend.chr
            coord2 = record_breakend.pos
        try:
            strands = record.INFO["STRANDS"]
            if isinstance(strands, (list, tuple)):
                strands = strands[0]
        except KeyError:
            if "ins" in svtype.lower() or "del" in svtype.lower():
                strands = "+-"
            else:
                raise ValueError("Unknown strands {missing STRANDS entry}")
        if "ins" in svtype.lower():
            coord2 = coord1 + 1
        add_record_ref_alt_to_extra(svtype, record, extra)
        strand1 = Strand.from_pm_string(string=strands[0])
        strand2 = Strand.from_pm_string(string=strands[1])
        if setup.get(STRIP_CHR, True):
            chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
        extra.update({
            EXTERNAL_NA_ID: record_id,
        })
        extra = clear_duplicated_entries_extra(extra=extra)
        pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
        pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
        na = Adjacency(position1=pos1, position2=pos2, extra=extra)
        nas_by_ids[record_id].append(na)
        processed_ids.add(record_id)
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return nas_by_ids.values()


def get_nas_from_survivor_vcf_records(survivor_vcf_records, setup=None, adjacencies_by_ids_by_sample_name=None, suffix_sample_extra=False, survivor_prefix=""):
    if setup is None:
        setup = {}
    if adjacencies_by_ids_by_sample_name is None:
        adjacencies_by_ids_by_sample_name = {}
    adjacencies_by_coordinates_by_sample_name = {}
    for sample, adj_by_ids in adjacencies_by_ids_by_sample_name.items():
        adjacencies_by_coordinates_by_sample_name[sample] = {}
        for adj_id, adj in adj_by_ids.items():
            adjacencies_by_coordinates_by_sample_name[sample][f'{adj.position1.chromosome}_{adj.position1.coordinate}-{adj.position2.chromosome}_{adj.position2.coordinate}'] = adj
    nas_by_ids = defaultdict(list)
    records_by_ids = get_vcf_records_by_ids(vcf_records=survivor_vcf_records)
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        svtype = record.INFO.get("SVTYPE", "")
        record_id = str(record.ID)
        chr1 = record.CHROM
        coord1 = int(record.POS)
        if "CHR2" in record.INFO:
            chr2 = record.INFO["CHR2"]
            coord2 = int(record.INFO["END"])
        else:
            record_breakend = record.ALT[0]
            chr2 = record_breakend.chr
            coord2 = record_breakend.pos
        try:
            strands = record.INFO["STRANDS"]
            if isinstance(strands, (list, tuple)):
                strands = strands[0]
        except KeyError:
            if svtype == "INS":
                strands = "+-"
            else:
                raise ValueError("Unknown strands {missing STRANDS entry}")
        strand1 = Strand.from_pm_string(string=strands[0])
        strand2 = Strand.from_pm_string(string=strands[1])
        if setup.get(STRIP_CHR, True):
            chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
        extra.update({EXTERNAL_NA_ID: record_id})
        extra = {key.lower(): value for key, value in extra.items()}
        pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
        pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
        supporting_source_ids = set()
        source_vector = []
        for vcf_sample in record.samples:
            source_vector.append(vcf_sample.sample)
            source_id = vcf_sample.data.ID.split(",")[0]
            if source_id == "NaN":
                continue
            source_ids = set()
            source_ids.add(source_id)
            coordinates = vcf_sample.data.CO.split(",")
            if len(coordinates) > 1:
                for coordinate_entry in coordinates:
                    if coordinate_entry in adjacencies_by_coordinates_by_sample_name[vcf_sample.sample]:
                        adj = adjacencies_by_coordinates_by_sample_name[vcf_sample.sample][coordinate_entry]
                        source_ids.add(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased))
            for source_id in sorted(source_ids):
                supporting_source_ids.add(source_id)
            if vcf_sample.sample in adjacencies_by_ids_by_sample_name:
                sample_source_adjacencies = adjacencies_by_ids_by_sample_name[vcf_sample.sample]
                source_adjacency = sample_source_adjacencies[source_id]
                for key, value in source_adjacency.extra.items():
                    if suffix_sample_extra:
                        key = str(vcf_sample.sample) + "_" + str(key)
                    if key not in extra:
                        extra[key] = value
        separator = "_" if len(survivor_prefix) > 0 else ""
        extra[survivor_prefix + separator + "supporting_source_ids"] = sorted(supporting_source_ids)
        extra[survivor_prefix + separator + "sources"] = source_vector
        extra[survivor_prefix + separator + "supporting_sources"] = [source for source, flag in zip(source_vector, extra.get("supp_vec", [0 for _ in source_vector])) if
                                                                     int(flag) != 0]
        na = Adjacency(position1=pos1, position2=pos2, extra=extra)
        nas_by_ids[record_id].append(na)
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return nas_by_ids.values()


def get_nas_from_svaba_vcf_records(svaba_vcf_records, source_type="sv", setup=None, samples=None, samples_all_any="any", samples_only=False):
    if setup is None:
        setup = {}
    nas_by_ids = defaultdict(list)
    records_by_ids = get_vcf_records_by_ids(vcf_records=svaba_vcf_records)
    processed_vcf_ids = set()
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        if samples is not None:
            present = {}
            for format_sample in record.samples:
                present[format_sample.sample] = "1" in format_sample.gt_alleles
            refute = False
            if samples_only:
                refute = any([value for key, value in present.items() if key not in samples])
            if samples_all_any == "any":
                present = any([present.get(sample_name, True) for sample_name in samples])
            else:
                present = all([present.get(sample_name, True) for sample_name in samples])
            present &= not refute
            if not present:
                continue
        if source_type == "indel":
            if 1 not in [len(record.REF), len(record.ALT)]:
                raise ValueError("Something wierd with record {id}".format(id=record.ID))
            if len(record.REF) > len(record.ALT):
                svtype = "DEL"
            else:
                svtype = "INS"
            chr1 = record.CHROM
            if setup.get(STRIP_CHR, True):
                chr1 = strip_chr(chr_string=chr1)
            chr2 = chr1
            coord1 = int(record.POS)
            coord2 = coord1 + 1
            if svtype == "DEL":
                coord2 += int(record.INFO["SPAN"])
            strand1 = Strand.FORWARD
            strand2 = Strand.REVERSE
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = str(record.ID)
            extra["svtype"] = svtype
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            add_record_ref_alt_to_extra(svtype, record, extra)
            nas_by_ids[str(record.ID)].append(na)
        elif source_type == "sv":
            if record.ID in processed_vcf_ids:
                continue
            processed_vcf_ids.add(record.ID)
            mate_vcf_record_id = record.INFO["MATEID"]
            if mate_vcf_record_id not in records_by_ids:
                continue
            processed_vcf_ids.add(mate_vcf_record_id)
            mate_vcf_record = records_by_ids[mate_vcf_record_id]
            if samples is not None:
                present = {}
                for format_sample in mate_vcf_record.samples:
                    present[format_sample.sample] = "1" in format_sample.gt_alleles
                refute = False
                if samples_only:
                    refute = any([value for key, value in present.items() if key not in samples])
                if samples_all_any == "any":
                    present = any([present.get(sample_name, True) for sample_name in samples])
                else:
                    present = all([present.get(sample_name, True) for sample_name in samples])
                present &= not refute
                if not present:
                    continue
            chr1 = record.CHROM
            chr2 = mate_vcf_record.CHROM
            if setup.get(STRIP_CHR, True):
                chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
            coord1 = int(record.POS)
            coord2 = int(mate_vcf_record.POS)
            strand1 = get_strand_from_alt_breakend(vcf_recrod=mate_vcf_record)
            strand2 = get_strand_from_alt_breakend(vcf_recrod=record)
            extra.update(mate_vcf_record.INFO)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = str(record.ID).split(":")[0]
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            nas_by_ids[extra[EXTERNAL_NA_ID]].append(na)
        else:
            raise ValueError("unknown SvABA input type {source_type}. Only 'sv' or 'indel' are allowed".format(source_type=source_type))
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return nas_by_ids.values()


def get_nas_from_grocsv_vcf_records(grocsv_vcf_records, setup=None, samples=None, samples_all_any="any", samples_only=False):
    if setup is None:
        setup = {}
    nas_by_ids = defaultdict(list)
    records_by_ids = get_vcf_records_by_ids(vcf_records=grocsv_vcf_records)
    processed_ids = set()
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        record_id = str(record.ID)
        if record_id in processed_ids:
            continue
        svtype = record.INFO["SVTYPE"]
        if svtype == "BND":
            mate_vcf_record_id = record.INFO["MATEID"]
            svid = record.INFO["EVENT"]
            mate_vcf_record = records_by_ids[mate_vcf_record_id]
            processed_ids.add(record_id)
            processed_ids.add(mate_vcf_record_id)
            if samples is not None:
                present = {}
                for format_sample in record.samples:
                    present[format_sample.sample] = str(format_sample.data.GT) in ["1", "."]
                refute = False
                if samples_only:
                    refute = any([value for key, value in present.items() if key not in samples])
                if samples_all_any == "any":
                    present = any([present.get(sample_name, True) for sample_name in samples])
                else:
                    present = all([present.get(sample_name, True) for sample_name in samples])
                present &= not refute
                if not present:
                    continue
            extra.update(mate_vcf_record.INFO)
            chr1 = record.CHROM
            chr2 = mate_vcf_record.CHROM
            coord1 = int(record.POS)
            coord2 = int(mate_vcf_record.POS)
            strand1 = get_strand_from_alt_breakend(vcf_recrod=mate_vcf_record)
            strand2 = get_strand_from_alt_breakend(vcf_recrod=record)
            if setup.get(STRIP_CHR, True):
                chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra[EXTERNAL_NA_ID] = svid
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            nas_by_ids[svid].append(na)
        else:
            raise Exception("Unknown SV type \"{svtype}\" for GROCSVS VCF entry".format(svtype=svtype))
    nas_by_ids = update_nas_ids(nas_by_ids_defaultdict=nas_by_ids, setup=setup)
    update_adjacencies_svtype(adjacencies=nas_by_ids.values())
    return nas_by_ids.values()


def get_nas_from_naibr_file(naibr_file_name, setup=None):
    if setup is None:
        setup = {}
    with open(naibr_file_name, "rt") as source:
        return get_nas_from_naibr_source(source=source, setup=setup)


def get_nas_from_naibr_source(source, setup=None):
    if setup is None:
        setup = {}
    result = defaultdict(list)
    csv_reader = csv.DictReader(source, delimiter="\t")
    for cnt, row in enumerate(csv_reader):
        chr1 = row["Chr1"]
        chr2 = row["Chr2"]
        if setup.get(STRIP_CHR, True):
            chr1 = strip_chr(chr_string=chr1)
            chr2 = strip_chr(chr_string=chr2)
        coord1 = int(row["Break1"])
        coord2 = int(row["Break2"])
        strands = row["Orientation"]
        strand1 = Strand.from_pm_string(string=strands[0])
        strand2 = Strand.from_pm_string(string=strands[1])
        pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
        pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
        na_id = "{cnt}".format(cnt=cnt)
        na = Adjacency(position1=pos1, position2=pos2, extra={EXTERNAL_NA_ID: na_id})
        result[na_id].append(na)
    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    update_adjacencies_svtype(adjacencies=result.values())
    return list(result.values())


def get_nas_from_delly_vcf_records(delly_vcf_records, setup=None):
    if setup is None:
        setup = {}
    result = defaultdict(list)
    records_by_ids = get_vcf_records_by_ids(vcf_records=delly_vcf_records)
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        svtype = record.INFO["SVTYPE"]
        svid = str(record.ID)
        extra[SVTYPE] = svtype
        extra[EXTERNAL_NA_ID] = svid
        if svtype in ["DEL", "DUP", "INV", "INS", "BND"]:
            chr1 = record.CHROM
            coord1 = int(record.POS)
            chr2 = record.INFO["CHR2"]
            coord2 = int(record.INFO["END"])
            ct_string = record.INFO["CT"]
            if setup.get(STRIP_CHR, True):
                chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
            if svtype in ["DEL", "INS"]:
                if "INSLEN" in record.INFO:
                    extra[SVLEN] = int(record.INFO["INSLEN"])
                strand1 = Strand.FORWARD
                strand2 = Strand.REVERSE
            elif svtype == "DUP":
                strand1 = Strand.REVERSE
                strand2 = Strand.FORWARD
            elif svtype in ["INV", "BND"]:
                strand1, strand2 = get_strands_from_CT_vcf_string(ct_string=ct_string)
            else:
                raise Exception("SVTYPE parsing issue")
            if svtype != "INS":
                assert (strand1, strand2) == get_strands_from_CT_vcf_string(ct_string=ct_string)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            add_record_ref_alt_to_extra(svtype, record, extra)
            extra = clear_duplicated_entries_extra(extra=extra)
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            result[svid].append(na)
        else:
            raise Exception("Unknown SV type \"{svtype}\" for Delly VCF input".format(svtype=svtype))
    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    return result.values()


def get_strands_from_CT_vcf_string(ct_string):
    strand1 = Strand.FORWARD if ct_string.startswith("3") else Strand.REVERSE
    strand2 = Strand.FORWARD if ct_string.endswith("3") else Strand.REVERSE
    return strand1, strand2


def delly_vcf_to_nas_stream(source, dest, setup=None, extra=None):
    if setup is None:
        setup = {}
    adjacencies = generate_nas_from_delly_vcf_records(source=source, setup=setup)
    write_adjacencies_to_destination(destination=dest, adjacencies=adjacencies, sort_adjacencies=False, extra=extra)


def generate_nas_from_delly_vcf_records(source, setup=None):
    vcf_reader = vcf.Reader(source)
    for record in vcf_reader:
        extra = deepcopy(record.INFO)
        svtype = record.INFO["SVTYPE"]
        svid = str(record.ID)
        extra[SVTYPE] = svtype
        extra[EXTERNAL_NA_ID] = svid
        if svtype in ["DEL", "DUP", "INV", "INS", "BND"]:
            chr1 = record.CHROM
            coord1 = int(record.POS)
            chr2 = record.INFO["CHR2"]
            coord2 = int(record.INFO["END"])
            ct_string = record.INFO["CT"]
            if setup.get(STRIP_CHR, True):
                chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
            if svtype in ["DEL", "INS"]:
                if "INSLEN" in record.INFO:
                    extra[SVLEN] = int(record.INFO["INSLEN"])
                strand1 = Strand.FORWARD
                strand2 = Strand.REVERSE
            elif svtype == "DUP":
                strand1 = Strand.REVERSE
                strand2 = Strand.FORWARD
            elif svtype in ["INV", "BND"]:
                strand1, strand2 = get_strands_from_CT_vcf_string(ct_string=ct_string)
            else:
                raise Exception("SVTYPE parsing issue")
            if svtype != "INS":
                assert (strand1, strand2) == get_strands_from_CT_vcf_string(ct_string=ct_string)
            pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
            pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
            extra = clear_duplicated_entries_extra(extra=extra)
            na = Adjacency(position1=pos1, position2=pos2, extra=extra)
            if setup.get(APPEND_ID, True):
                if not na.extra[EXTERNAL_NA_ID].endswith(setup.get(ID_SUFFIX, "")):
                    na.extra[EXTERNAL_NA_ID] += "_{suffix}".format(suffix=setup.get(ID_SUFFIX, ""))
            yield na
        else:
            raise Exception("Unknown SV type \"{svtype}\" for Delly VCF input".format(svtype=svtype))


def get_nas_from_pbsv_vcf_records(pbsv_vcf_records, setup=None, sample=None, silent_skip_unknown_sv=True):
    if setup is None:
        setup = {}
    result = defaultdict(list)
    processed_ids = set()
    records_by_ids = get_vcf_records_by_ids(vcf_records=pbsv_vcf_records)
    for record in records_by_ids.values():
        extra = deepcopy(record.INFO)
        svid = str(record.ID)
        if svid in processed_ids:
            continue
        svtype = record.INFO["SVTYPE"]
        extra[SVTYPE] = svtype
        extra[VCF_FILTER] = get_string_vcf_filter(record=record)
        extra[EXTERNAL_NA_ID] = svid
        if svtype in ["DEL", "INS", "INV", "BND", "DUP"]:
            for vcf_sample in record.samples:
                if sample is None or vcf_sample.sample.lower() == sample.lower():
                    extra["OR_GT"] = "/".join(map(str, vcf_sample.gt_alleles)).replace("None", ".")
                    support, total = 0, 0
                    l, r = vcf_sample.gt_alleles
                    if l in ['1', '.']:
                        support += vcf_sample.data.AD[0]
                    if r in ['1', '.']:
                        support += vcf_sample.data.AD[1]
                    total = vcf_sample.data.DP
                    extra["re"] = support
                    extra["rcnt"] = total
                    break
            if svtype in ["DEL", "INS", "DUP"]:
                chr1 = str(record.CHROM)
                coord1 = int(record.POS)
                if setup.get(STRIP_CHR, True):
                    chr1 = strip_chr(chr_string=chr1)
                chr2 = chr1
                coord2 = record.INFO["END"]
                if svtype in ["DEL", "INS"]:
                    strand1 = Strand.FORWARD
                    strand2 = Strand.REVERSE
                elif "ins" in svid.lower():
                    strand1 = Strand.FORWARD
                    strand2 = Strand.REVERSE
                    coord1 = coord2
                    coord2 = coord1 + 1
                    assert len(record.INFO["SVLEN"]) == 1
                    extra[SVLEN] = int(record.INFO["SVLEN"][0])
                    del extra["SVLEN"]
                else:
                    strand1 = Strand.REVERSE
                    strand2 = Strand.FORWARD
                if svtype == "INS" and coord1 == coord2:
                    coord2 += 1
                if svtype == "INS":
                    assert len(record.INFO["SVLEN"]) == 1
                    extra[SVLEN] = int(record.INFO["SVLEN"][0])
                pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
                pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
                extra = clear_duplicated_entries_extra(extra=extra)
                add_record_ref_alt_to_extra(svtype, record, extra)
                na = Adjacency(position1=pos1, position2=pos2, extra=extra)
                result[svid].append(na)
                processed_ids.add(svid)
            elif svtype == "INV":
                chr1 = str(record.CHROM)
                if setup.get(STRIP_CHR, True):
                    chr1 = strip_chr(chr_string=chr1)
                coord11, coord12 = get_inv_na1_coords_from_vcf_record(vcf_record=record)
                coord21, coord22 = get_inv_na2_coords_from_vcf_record(vcf_record=record)
                strand11, strand12 = get_inv_na1_strands_from_vcf_record(vcf_record=record)
                strand21, strand22 = get_inv_na2_strands_from_vcf_record(vcf_record=record)
                pos11 = Position(chromosome=chr1, coordinate=coord11, strand=strand11)
                pos12 = Position(chromosome=chr1, coordinate=coord12, strand=strand12)
                pos21 = Position(chromosome=chr1, coordinate=coord21, strand=strand21)
                pos22 = Position(chromosome=chr1, coordinate=coord22, strand=strand22)
                extra1 = deepcopy(extra)
                extra2 = deepcopy(extra)
                extra1[EXTERNAL_NA_ID] = extra1[EXTERNAL_NA_ID] + "_1"
                extra2[EXTERNAL_NA_ID] = extra2[EXTERNAL_NA_ID] + "_2"
                extra1 = clear_duplicated_entries_extra(extra=extra1)
                extra2 = clear_duplicated_entries_extra(extra=extra2)
                na1 = Adjacency(position1=pos11, position2=pos12, extra=extra1)
                na2 = Adjacency(position1=pos21, position2=pos22, extra=extra2)
                result[svid].append(na1)
                result[svid].append(na2)
                processed_ids.add(svid)
            else:
                assert len(record.INFO["MATEID"]) == 1
                mate_vcf_id = record.INFO["MATEID"][0]
                if mate_vcf_id not in records_by_ids:
                    continue
                mate_record = records_by_ids[mate_vcf_id]
                extra.update(mate_record.INFO)
                chr1 = str(record.CHROM)
                coord1 = int(record.POS)
                chr2 = str(mate_record.CHROM)
                coord2 = int(mate_record.POS)
                strand1 = get_strand_from_alt_breakend(vcf_recrod=mate_record)
                strand2 = get_strand_from_alt_breakend(vcf_recrod=record)
                if setup.get(STRIP_CHR, True):
                    chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
                pos1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
                pos2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
                extra[EXTERNAL_NA_ID2] = mate_vcf_id
                extra = clear_duplicated_entries_extra(extra=extra)
                na = Adjacency(position1=pos1, position2=pos2, extra=extra)
                result[svid].append(na)
                processed_ids.add(svid)
                processed_ids.add(mate_vcf_id)
        else:
            if silent_skip_unknown_sv:
                continue
            else:
                raise Exception("Unknonw SV type \"{svtype}\" for PBSV VCF input".format(svtype=svtype))
    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    update_adjacencies_svtype(adjacencies=result.values())
    return result.values()


GUNDEM_SAMPLE_NAME = "patient"
GUNDEM_CHR1 = "chr1"
GUNDEM_CHR2 = "chr2"
GUNDEM_COORDINATE1 = "pos1"
GUNDEM_COORDINATE2 = "pos2"
GUNDEM_STRAND1 = "str1"
GUNDEM_STRAND2 = "str2"
GUNDEM_SAMPLES_LIST = "samples"
GUNDEM_SAMPLES_COUNTS = "tumour-counts"
GUNDEM_PATIENT = "patient"

GUNDEM_PER_SAMPLE_SUPPORT = "per_sample_support"


def gundem2015_flip_strand(strand_string):
    if strand_string == "-":
        return "+"
    elif strand_string == "+":
        return "-"
    raise ValueError("{} is not a valid strand string".format(strand_string))


def get_nas_from_gundem2015_file(file_name, separator="\t", flip_second_strand=True):
    with open(file_name, "rt") as source:
        return get_nas_from_gundem2015_source(source=source, separator=separator, flip_second_strand=flip_second_strand)


def get_nas_from_gundem2015_source(source, setup=None, separator="\t", flip_second_strand=True):
    if setup is None:
        setup = {}
    result = defaultdict(list)
    reader = csv.DictReader(source, delimiter=separator)
    for cnt, row in enumerate(reader):
        patient_name = row[GUNDEM_PATIENT]
        sample_list = row[GUNDEM_SAMPLES_LIST].split("/")
        sample_list = [patient_name + "-" + sample_name for sample_name in sample_list]
        sample_cnts = list(map(int, row[GUNDEM_SAMPLES_COUNTS].split("/")))
        per_sample_counts = {sample: cnt for sample, cnt in zip(sample_list, sample_cnts)}
        chromosome1, chromosome2 = row[GUNDEM_CHR1], row[GUNDEM_CHR2]
        if setup.get(STRIP_CHR, True):
            chromosome1, chromosome2 = strip_chr(chr_string=chromosome1), strip_chr(chr_string=chromosome2)
        position1 = Position(chromosome=chromosome1, coordinate=int(row[GUNDEM_COORDINATE1]), strand=Strand.from_pm_string(row[GUNDEM_STRAND1]))
        strand2_string = gundem2015_flip_strand(row[GUNDEM_STRAND2]) if flip_second_strand else row[GUNDEM_STRAND2]
        position2 = Position(chromosome=chromosome2, coordinate=int(row[GUNDEM_COORDINATE2]), strand=Strand.from_pm_string(strand2_string))
        aid = str(cnt)

        extra = {
            EXTERNAL_NA_ID: aid,
            GUNDEM_PER_SAMPLE_SUPPORT: per_sample_counts,
        }
        novel_adjacency = Adjacency(position1=position1, position2=position2, adjacency_type=AdjacencyType.NOVEL, extra=extra)
        result[aid].append(novel_adjacency)

    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    return result.values()


REMIXT_PREDICTION_ID = "prediction_id"
REMIXT_CHROMOSOME_1 = "chromosome_1"
REMIXT_POSITION_1 = "position_1"
REMIXT_STRAND_1 = "strand_1"
REMIXT_CHROMOSOME_2 = "chromosome_2"
REMIXT_POSITION_2 = "position_2"
REMIXT_STRAND_2 = "strand_2"
REMIXT_CN_1 = "cn_1"
REMIXT_CN_2 = "cn_2"


def get_nas_from_remixt_file(file_name, setup=None, separator="\t", clone_ids=None, remixt_na_correction=True, skip_absent=False):
    with open(file_name, "rt") as source:
        return get_nas_from_remixt_source(source=source, setup=setup, separator=separator, clone_ids=clone_ids, remixt_na_correction=remixt_na_correction, skip_absent=skip_absent)


def get_nas_from_remixt_source(source, setup=None, separator="\t", clone_ids=None, remixt_na_correction=True, skip_absent=False):
    if clone_ids is None:
        clone_ids = ["1", "2"]
    if setup is None:
        setup = {}
    result = defaultdict(list)
    reader = csv.DictReader(source, delimiter=separator)
    for line in reader:
        naid = line[REMIXT_PREDICTION_ID]
        chr1_name = line[REMIXT_CHROMOSOME_1]
        chr1_coordinate = int(line[REMIXT_POSITION_1])
        chr1_strand = Strand.from_pm_string(line[REMIXT_STRAND_1])
        chr2_name = line[REMIXT_CHROMOSOME_2]
        chr2_coordinate = int(line[REMIXT_POSITION_2])
        chr2_strand = Strand.from_pm_string(line[REMIXT_STRAND_2])
        clone1_cn = int(line[REMIXT_CN_1])
        clone2_cn = int(line[REMIXT_CN_2])
        cn_data = {}
        skip = True
        for clone_cn, clone_name in zip([clone1_cn, clone2_cn], ["1", "2"]):
            if clone_name in clone_ids:
                cn_data[clone_name] = {str(Phasing.AA): clone2_cn}
                skip &= clone2_cn == 0
        skip = skip and skip_absent
        if skip:
            continue
        if remixt_na_correction and chr1_strand == Strand.FORWARD:
            chr1_coordinate -= 1
        if remixt_na_correction and chr2_strand == Strand.FORWARD:
            chr2_coordinate -= 1
        if setup.get(STRIP_CHR, True):
            chr1_name, chr2_name = strip_chr(chr_string=chr1_name), strip_chr(chr_string=chr2_name)
        position1 = Position(chromosome=chr1_name, coordinate=chr1_coordinate, strand=chr1_strand)
        position2 = Position(chromosome=chr2_name, coordinate=chr2_coordinate, strand=chr2_strand)
        novel_adjacency = Adjacency(position1=position1, position2=position2,
                                    adjacency_type=AdjacencyType.NOVEL, extra={EXTERNAL_NA_ID: naid, COPY_NUMBER: cn_data})
        aid = novel_adjacency.stable_id_non_phased
        result[aid].append(novel_adjacency)
    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    return result.values()


BREAKDANCER_SCORE = "score"
BREAKDANCER_NUM_READS = "num_reads"
BREAKDANCER_NUM_READS_PER_LIB = "num_reads_lib"
BREAKDANCER_INTR_CHROMOSOMAL_SV_TYPE = "CTX"


def get_strand_from_breakdancer_strand_string(strand_string):
    pcntstr, ncntstr, _ = re.split("[+-]", strand_string)
    pcnt, ncnt = int(pcntstr), int(ncntstr)
    return Strand.FORWARD if pcnt >= ncnt else Strand.REVERSE


def get_nas_from_breakdancer_file(file_name, setup=None):
    with open(file_name, "rt") as source:
        return get_nas_from_breakdancer_source(source=source, setup=setup)


def get_nas_from_breakdancer_source(source, setup=None):
    if setup is None:
        setup = {}
    result = defaultdict(list)
    for cnt, line in enumerate(source):
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        data = line.split("\t")
        chr1 = data[0]
        coord1 = int(data[1])
        strand1 = get_strand_from_breakdancer_strand_string(strand_string=data[2])
        chr2 = data[3]
        coord2 = int(data[4])
        strand2 = get_strand_from_breakdancer_strand_string(strand_string=data[5])
        if setup.get(STRIP_CHR, True):
            chr1, chr2 = strip_chr(chr_string=chr1), strip_chr(chr_string=chr2)
        sv_type = data[6]
        size = int(data[7])
        if sv_type == BREAKDANCER_INTR_CHROMOSOMAL_SV_TYPE:
            size = -1

        score = int(float(data[8]))
        read_cnt = int(float(data[9]))
        reads_per_lib_string = data[10]
        naid = str(cnt)
        position1 = Position(chromosome=chr1, coordinate=coord1, strand=strand1)
        position2 = Position(chromosome=chr2, coordinate=coord2, strand=strand2)
        novel_adjacency = Adjacency(position1=position1, position2=position2,
                                    adjacency_type=AdjacencyType.NOVEL, extra={EXTERNAL_NA_ID: naid, SVTYPE: sv_type,
                                                                               SVLEN: size, BREAKDANCER_SCORE: score,
                                                                               BREAKDANCER_NUM_READS: read_cnt, BREAKDANCER_NUM_READS_PER_LIB: reads_per_lib_string})
        aid = novel_adjacency.stable_id_non_phased
        result[aid].append(novel_adjacency)
    result = update_nas_ids(nas_by_ids_defaultdict=result, setup=setup)
    update_adjacencies_svtype(adjacencies=result.values())
    return result.values()


def build_setup(args):
    result = {
        STRIP_CHR: args.strip_chr,
        APPEND_ID: args.append_id_suffix,
        EXTRA_FIELDS: args.o_extra_fields.split(","),
        ID_SUFFIX: args.id_suffix,
    }
    return result


def get_chrs_regions_string_list_from_file(file_name):
    with open(file_name, "rt") as source:
        return get_chrs_regions_string_lists_from_source(source=source)


def get_chrs_regions_string_lists_from_source(source):
    result = []
    for line in source:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        result.append(line)
    return result


def parse_segment_chr_region(string):
    if ":" not in string.split()[0] and "\t" not in string:
        return Segment.from_chromosome_coordinates(chromosome=string.strip(), start=0, end=3000000000)
    if ":" in string.split()[0]:
        chromosome, start_end = string.split(":")
        start, end = start_end.split("-")
        start, end = int(start), int(end)
        return Segment.from_chromosome_coordinates(chromosome=chromosome.lower(), start=start, end=end)
    data = string.split("\t")
    chromosome, start, end = data[:3]
    start, end = int(start), int(end)
    segment = Segment.from_chromosome_coordinates(chromosome=chromosome.lower(), start=start, end=end)
    if len(data) > 3:
        segment.extra.update(parse_segment_extra(extra=data[3], extra_separator=";"))
    return segment
