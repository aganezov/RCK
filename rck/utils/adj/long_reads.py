import logging
from collections import defaultdict
import networkx as nx
import pysam

from rck.core.io import read_adjacencies_from_source

#####


def get_mode_str(format="bam", input=False):
    result = "r" if input else "w"
    if format == "bam":
        result += "b"
    elif format == "cram":
        result += "c"
    return result

#####


def extract_long_reads():
    pass


def filter_alignment():
    pass


def infer_labeling_constraints(rck_nas_source, alignment_file, i_alignment_format, lr_field, min_sv_cnt, logger=None):
    logger = logger or logging.getLogger('dummy')
    logger.info("RCK NAs file object {file_name}".format(file_name=str(rck_nas_source)))
    nas = read_adjacencies_from_source(source=rck_nas_source)
    reads_to_nas = defaultdict(list)
    for na in nas:
        reads_str = na.extra.get(lr_field, "")
        reads = reads_str.split(",")
        for read in reads:
            if len(read) == 0:
                continue
            reads_to_nas[read].append(na)
    logger.debug("{reads_cnt} -- number of reads".format(reads_cnt=len(reads_to_nas)))
    reads = {read for read in reads_to_nas if len(reads_to_nas[read]) >= min_sv_cnt}
    logger.debug("{reads_cnt} -- number of reads that each span {min_sv_cnt}+ NAs".format(reads_cnt=len(reads), min_sv_cnt=min_sv_cnt))
    location_graph = nx.Graph()
    mode = get_mode_str(format=i_alignment_format, input=True)
    current_read_name = None
    current_entries = []
    with pysam.AlignmentFile(alignment_file, mode) as i_stream:
        if "SO:queryname" not in i_stream.text:
            logger.critical("Input alignment file {alignment_file} is not sorted by read (i.e., query) name".format(alignment_file=alignment_file))
            raise Exception("Input bam file needs to be sorted by read (i.e., query) name")
        for entry in i_stream:
            if entry.qname != current_read_name:
                if len(current_entries) > 0:
                    reads_novel_adjacencies = reads_to_nas[current_read_name]
                    add_lr_labeling_constraints(location_graph=location_graph, alignment_entries=current_entries, nas=reads_novel_adjacencies)
                current_read_name = entry.qname
                current_entries = [entry]
            else:
                current_entries.append(entry)
        # last reads streak has to be processed as well
        if current_read_name in reads_to_nas:
            pass


def add_lr_labeling_constraints(location_graph, alignment_entries, nas):
    entries = sorted(alignment_entries, key=lambda e: (e.query_alignment_start, e.query_alignment_end))
    internal = len(nas) != len(entries) - 1

    pass


def label_constraints_combining():
    pass

