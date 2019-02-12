import argparse

import os
import sys
from collections import defaultdict

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import get_standard_logger_from_args, read_adjacencies_from_source, write_adjacency_groups_to_destination, get_logging_cli_parser
from rck.utils.adj.adjacency_group_inference import infer_sniffles_molecule_groups, infer_short_nas_labeling_groups, infer_alignment_labeling_groups, filter_alignment
from rck.utils.adj.adjacency_group_process import refined_labeling_groups


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADJ-GROUPS-infer")
    parser.add_argument("--version", action="version", version=rck.version)
    cli_logging_parser = get_logging_cli_parser()

    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ###
    sniffles_molecule_group_parser = subparsers.add_parser("sniffles-m", parents=[cli_logging_parser])
    sniffles_molecule_group_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    sniffles_molecule_group_parser.add_argument("--i-separator", default="\t")
    sniffles_molecule_group_parser.add_argument("--i-extra-separator", default=";")
    sniffles_molecule_group_parser.add_argument("--extra-rnames-field", default="rnames")
    sniffles_molecule_group_parser.add_argument("--fp", type=float, default=0.5)
    sniffles_molecule_group_parser.add_argument("--gid-suffix", dest="gid_suffix", default="sniffles-M")
    sniffles_molecule_group_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    sniffles_molecule_group_parser.add_argument("--o-separator", default="\t")
    sniffles_molecule_group_parser.add_argument("--o-aids-separator", default=",")
    sniffles_molecule_group_parser.add_argument("--o-extra-separator", default=";")
    ###
    short_nas_labeling_group_parser = subparsers.add_parser("short-l", parents=[cli_logging_parser])
    short_nas_labeling_group_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    short_nas_labeling_group_parser.add_argument("--i-separator", default="\t")
    short_nas_labeling_group_parser.add_argument("--i-extra-separator", default=";")
    short_nas_labeling_group_parser.add_argument("--max-size", type=int, default=1000)
    short_nas_labeling_group_parser.add_argument("--allow-intermediate-same", action="store_true", dest="allow_intermediate_same")
    short_nas_labeling_group_parser.add_argument("--allow-intermediate-tra", action="store_true", dest="allow_intermediate_tra")
    short_nas_labeling_group_parser.add_argument("--no-refine", action="store_false", dest="refine")
    short_nas_labeling_group_parser.add_argument("--fp", type=float, default=1)
    short_nas_labeling_group_parser.add_argument("--gid-suffix", dest="gid_suffix", default="short-nas-L")
    short_nas_labeling_group_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    short_nas_labeling_group_parser.add_argument("--o-separator", default="\t")
    short_nas_labeling_group_parser.add_argument("--o-aids-separator", default=",")
    short_nas_labeling_group_parser.add_argument("--o-extra-separator", default=";")
    ###
    sniffles_labeling_group_parser = subparsers.add_parser("sniffles-l", parents=[cli_logging_parser])
    sniffles_labeling_group_parser.add_argument("--rck-adj", type=argparse.FileType("rt"), required=True)
    sniffles_labeling_group_parser.add_argument("--i-separator", default="\t")
    sniffles_labeling_group_parser.add_argument("--i-extra-separator", default=";")
    sniffles_labeling_group_parser.add_argument("--alignment", required=True)
    sniffles_labeling_group_parser.add_argument("--alignment-format", choices=["sam", "bam", "cram"], default="bam")
    sniffles_labeling_group_parser.add_argument("--extra-rnames-field", default="rnames")
    sniffles_labeling_group_parser.add_argument("--no-refine", action="store_false", dest="refine")
    sniffles_labeling_group_parser.add_argument("--fp", type=float, default=1)
    sniffles_labeling_group_parser.add_argument("--gid-suffix", default="sniffles-L")
    sniffles_labeling_group_parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    sniffles_labeling_group_parser.add_argument("--o-separator", default="\t")
    sniffles_labeling_group_parser.add_argument("--o-aids-separator", default=",")
    sniffles_labeling_group_parser.add_argument("--o-extra-separator", default=";")
    ###
    filter_alignment_parser = subparsers.add_parser("filter-alignment", parents=[cli_logging_parser])
    filter_alignment_parser.add_argument("--rck-adj", type=argparse.FileType("rt"), required=True)
    filter_alignment_parser.add_argument("--i-separator", default="\t")
    filter_alignment_parser.add_argument("--i-extra-separator", default=";")
    filter_alignment_parser.add_argument("--extra-rnames-field", default="rnames")
    filter_alignment_parser.add_argument("--alignment", required=True)
    filter_alignment_parser.add_argument("--alignment-format", choices=["sam", "bam", "cram"], default="bam")
    filter_alignment_parser.add_argument("-o", "--output", required=True)
    filter_alignment_parser.add_argument("--output-format", choices=["sam", "bam", "cram"], default="bam")
    ###
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-ADJ-GROUPS-infer")
    if args.command == "sniffles-m":
        logger.info("Inferring molecule adjacency groups from adjacencies with Sniffles RNAMES support extra info.")
        logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
        adjacencies = read_adjacencies_from_source(source=args.rck_adj, separator=args.i_separator, extra_separator=args.i_extra_separator)
        logger.info("Inferring molecule adjacency groups from read adjacencies")
        adj_groups = infer_sniffles_molecule_groups(adjacencies=adjacencies, extra_rnames_field=args.extra_rnames_field, gid_suffix=args.gid_suffix)
        logger.info("Inferred {cnt} molecule adjacency groups".format(cnt=len(adj_groups)))
        logger.info("Writing inferred molecule adjacency groups to {file}".format(file=args.output))
        write_adjacency_groups_to_destination(destination=args.output, adjacency_groups=adj_groups,
                                              separator=args.o_separator, extra_separator=args.o_extra_separator, aids_separator=args.o_aids_separator,
                                              extra_fill="")
    elif args.command == "short-l":
        logger.info("Inferring labeling adjacency groups from adjacencies from adjacencies.")
        logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
        adjacencies = read_adjacencies_from_source(source=args.rck_adj, separator=args.i_separator, extra_separator=args.i_extra_separator)
        logger.info("Inferring labeling adjacency groups from read adjacencies")
        adj_groups = infer_short_nas_labeling_groups(adjacencies=adjacencies, gid_suffix=args.gid_suffix, max_size=args.max_size,
                                                     allow_intermediate_same=args.allow_intermediate_same,
                                                     allow_intermediate_tra=args.allow_intermediate_tra)
        logger.info("Inferred {cnt} labeling adjacency groups".format(cnt=len(adj_groups)))
        if args.refine:
            logger.info("Refining inferred labeling adjacency groups")
            adj_groups = refined_labeling_groups(adj_groups=adj_groups, gid_suffix=args.gid_suffix)
        logger.info("A total of {cnt} refined labeling adjacency groups remain".format(cnt=len(adj_groups)))
        logger.info("Writing inferred labeling adjacency group s to {file}".format(file=args.output))
        write_adjacency_groups_to_destination(destination=args.output, adjacency_groups=adj_groups,
                                              separator=args.o_separator, aids_separator=args.o_aids_separator, extra_separator=args.o_extra_separator,
                                              extra_fill="")
    elif args.command == "sniffles-l":
        logger.info("Inferring labeling adjacency groups from adjacencies, and their reads-of-origin alignments")
        logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
        adjacencies = read_adjacencies_from_source(source=args.rck_adj, extra_separator=args.i_extra_separator, separator=args.i_separator)
        logger.info("Inferring labeling adjacency groups from read adjacencies and their reads-of-origin alignments")
        adj_groups = infer_alignment_labeling_groups(adjacencies=adjacencies, alignment_file_name=args.alignment, alignment_format=args.alignment_format,
                                                     extra_rnames_field=args.extra_rnames_field, gid_suffix=args.gid_suffix)
        logger.info("Inferred {cnt} labeling adjacency groups. There can be many duplicates, refinement shall take care of it.".format(cnt=len(adj_groups)))
        if args.refine:
            logger.info("Refining inferred labeling adjacency groups")
            adj_groups = refined_labeling_groups(adj_groups=adj_groups, gid_suffix=args.gid_suffix)
        logger.info("A total of {cnt} refined labeling adjacency groups remain".format(cnt=len(adj_groups)))
        logger.info("Writing inferred labeling adjacency group s to {file}".format(file=args.output))
        write_adjacency_groups_to_destination(destination=args.output, adjacency_groups=adj_groups,
                                              separator=args.o_separator, aids_separator=args.o_aids_separator, extra_separator=args.o_extra_separator,
                                              extra_fill="")
    elif args.command == "filter-alignment":
        logger.info("Filtering input read alignment to retain only reads mentioned as supporting adjacencies from the input")
        logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
        adjacencies = read_adjacencies_from_source(source=args.rck_adj, extra_separator=args.i_extra_separator, separator=args.i_separator)
        logger.info("Filtering input alignment form file {file} and writing result in {o_file}".format(file=args.alignment, o_file=args.output))
        filter_alignment(adjacencies=adjacencies, alignment_file_name=args.alignment, alignment_format=args.alignment_format, extra_rnames_field=args.extra_rnames_field,
                         output_alignment_file_name=args.output, output_alignment_format=args.output_format)
        exit(0)


if __name__ == "__main__":
    main()
