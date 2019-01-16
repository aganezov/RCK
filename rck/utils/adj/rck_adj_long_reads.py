import argparse
import logging
import sys
from collections import defaultdict
import pysam
import networkx as nx

from rck.core.io import read_adjacencies_from_source, get_logging_cli_parser, get_standard_logger_from_args
from rck.utils.adj.long_reads import infer_labeling_constraints


def get_mode_str(format="bam", input=False):
    result = "r" if input else "w"
    if format == "bam":
        result += "b"
    elif format == "cram":
        result += "c"
    return result


def get_reads_set_from_source(source):
    reads = set()
    for line in source:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        reads.add(line)
    return reads


def main():
    parser = argparse.ArgumentParser()
    logging_parser = get_logging_cli_parser()
    ########
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ########
    lr_extraction_parser = subparsers.add_parser("extract-lr", parents=[logging_parser])
    lr_extraction_parser.add_argument("rck_nas", type=argparse.FileType("rt"), default=sys.stdin)
    lr_extraction_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    lr_extraction_parser.add_argument("--min-sv-cnt", type=int, default=2)
    lr_extraction_parser.add_argument("--lr-field", default="support_read_names")
    #########
    lr_alignment_filter_parser = subparsers.add_parser("filter-alignment", parents=[logging_parser])
    lr_alignment_filter_parser.add_argument("alignment", nargs="?", type=str, default="-")
    lr_alignment_filter_parser.add_argument("--i-alignment-format", type=str, choices=["bam", "sam", "cram"], default="bam")
    lr_alignment_filter_parser.add_argument("-r", "--reads", type=argparse.FileType("rt"), required=True)
    lr_alignment_filter_parser.add_argument("--r-separator", default="\t")
    lr_alignment_filter_parser.add_argument("--s-separator", default="\t")
    lr_alignment_filter_parser.add_argument("-o", "--output", type=str, default="-")
    lr_alignment_filter_parser.add_argument("--o-alignment-format", type=str, choices=["bam", "sam", "cram"], default="bam")
    #########
    labeling_constraint_inference_parser = subparsers.add_parser("label-const-inf", parents=[logging_parser])
    labeling_constraint_inference_parser.add_argument("alignment", type=str, default="-")
    labeling_constraint_inference_parser.add_argument("--i-alignment-format", type=str, choices=["bam", "sam", "cram"], default="bam")
    labeling_constraint_inference_parser.add_argument("--rck-nas", type=argparse.FileType("rt"), required=True)
    labeling_constraint_inference_parser.add_argument("--min-sv-cnt", type=int, default=2)
    labeling_constraint_inference_parser.add_argument("--lr-field", default="support_read_names")
    labeling_constraint_inference_parser.add_argument("-o", "--output", type=argparse.FileType("rt"), default=sys.stdout)
    #########
    labeling_constraint_combine_parser = subparsers.add_parser("label-const-com", parents=[logging_parser])
    labeling_constraint_combine_parser.add_argument("label-constr", type=argparse.FileType("rt"), nargs="+")
    labeling_constraint_combine_parser.add_argument("-o", "--output", type=argparse.FileType("rt"), default=sys.stdout)
    #########
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-LR")
    if args.command == "extract-lr":
        nas = read_adjacencies_from_source(source=args.rck_nas)
        reads_to_nas = defaultdict(list)
        for na in nas:
            reads_str = na.extra.get(args.lr_field, "")
            reads = reads_str.split(",")
            for read in reads:
                if len(read) == 0:
                    continue
                reads_to_nas[read].append(na)
        extracted_read_names = {read for read in reads_to_nas if len(reads_to_nas[read]) >= args.min_sv_cnt}
        for read_name in extracted_read_names:
            print(read_name, file=args.output)
    elif args.command == "filter-alignment":
        reads = get_reads_set_from_source(source=args.reads)
        imode = get_mode_str(format=args.i_alignment_format, input=True)
        omode = get_mode_str(format=args.o_alignment_format, input=False)
        with pysam.AlignmentFile(args.alignment, imode) as i_stream:
            with pysam.AlignmentFile(args.output, omode, template=i_stream) as o_stream:
                for entry in i_stream:
                    if entry.qname in reads:
                        o_stream.write(entry)
    elif args.command == "label-const-inf":
        constraints = infer_labeling_constraints(rck_nas_source=args.rck_nas, alignment_file=args.alignment, i_alignment_format=args.i_alignment_format,
                                                 lr_field=args.lr_field, min_sv_cnt=args.min_sv_cnt, logger=logger)

    elif args.command == "label-constr-com":
        pass

if __name__ == "__main__":
    main()
