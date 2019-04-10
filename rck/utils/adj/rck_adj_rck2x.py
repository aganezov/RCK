import argparse
import os
import sys


current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, read_adjacencies_from_source, write_adjacencies_to_vcf_sniffles_destination, \
    write_adjacencies_to_circa_destination, read_chr_sizes_from_source, write_segments_to_circa_destination
from rck.core.structures import AdjacencyType
from rck.utils.adj.process import get_circa_adj_cnt


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-NAS-rck2x")
    parser.add_argument('--version', action='version', version=rck.version)
    cli_logging_parser = get_logging_cli_parser()
    ###
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ###
    vcf_parser = subparsers.add_parser("vcf-sniffles", parents=[cli_logging_parser], help="Convert RCK Adjacencies to the VCF (Sniffles) format")
    vcf_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    vcf_parser.add_argument("--separator", default="\t")
    vcf_parser.add_argument("--extra-separator", default=";")
    vcf_parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    vcf_parser.add_argument("--o-extra-fields", default="all")
    vcf_parser.add_argument("--o-no-include-ref", action="store_false", dest="include_ref")
    vcf_parser.add_argument("--clone-suffix", default="")
    vcf_parser.add_argument("--dummy-clone", default="dummy_clone")
    ###
    circa_parser = subparsers.add_parser("circa", parents=[cli_logging_parser], help="Convert RCK Adjacencies to the TSV format supported by Circa")
    circa_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    circa_parser.add_argument("--separator", default="\t")
    circa_parser.add_argument("--extra-separator", default=";")
    circa_parser.add_argument("--size-extra-field")
    circa_parser.add_argument("--size-extra-field-no-abs", action="store_false", dest="size_extra_field_abs")
    circa_parser.add_argument("--size-extra-seq-field")
    circa_parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    ###
    circa_density_parser = subparsers.add_parser("circa-dens", parents=[cli_logging_parser],
                                                 help="Convert RCK Adjacencies to the TSV format with adjacencies density cnt per window supported by Circa")
    circa_density_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    circa_density_parser.add_argument("--separator", default="\t")
    circa_density_parser.add_argument("--extra-separator", default=";")
    circa_density_parser.add_argument("--window-size", type=int, default=10000000)
    circa_density_parser.add_argument("--chr-sizes", type=argparse.FileType("rt"))
    circa_density_parser.add_argument("--element", choices=["breakend", "adj"], default="breakend")
    circa_density_parser.add_argument("--element-adj-cnt-full", action="store_true", dest="circa_element_adj_cnt_full")
    circa_density_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ###
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-NAS-rck2x")
    logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
    adjacencies = read_adjacencies_from_source(source=args.rck_adj, extra_separator=args.extra_separator, separator=args.separator)
    if args.command == "vcf-sniffles":
        if not args.include_ref:
            logger.debug("Reference adjacencies were excluded from the output.")
            adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, adjacencies))
        if args.o_extra_fields is None or len(args.o_extra_fields) == 0 or args.o_extra_fields == ",":
            extra = None
        elif args.o_extra_fields != "all":
            extra = args.o_extra_fields.split(",")
        else:
            extra = args.o_extra_fields
        logger.debug("Output extra fields are identified as {o_extra}".format(o_extra=",".join(extra) if extra is not None else ""))
        logger.info("Converting RCK formatted adjacencies to the VCF (Sniffles) format")
        logger.info("Writing adjacencies to {file}".format(file=args.output))
        write_adjacencies_to_vcf_sniffles_destination(destination=args.output, adjacencies=adjacencies, extra=extra,
                                                      dummy_clone=args.dummy_clone, clone_suffix=args.clone_suffix)
    elif args.command == "circa":
        logger.info("Converting input RCK formatted adjacencies into a Circa suitable format (extra column get transformed into a size column)")
        logger.info("Writing adjacencies info suitable for Circa to {file}".format(file=args.output))
        write_adjacencies_to_circa_destination(destination=args.output, adjacencies=adjacencies, size_extra_field=args.size_extra_field,
                                               size_extra_seq_field=args.size_extra_seq_field, size_abs=args.size_extra_field_abs)
    elif args.command == "circa-dens":
        logger.info("Computing cnt of input RCK formatted adjacencies per window into a CIRCA suitable format")
        chr_sizes = args.chr_sizes
        if args.chr_sizes is not None:
            chr_sizes = read_chr_sizes_from_source(source=args.chr_sizes)
        circa_adj_cnts = get_circa_adj_cnt(adjacencies=adjacencies, window_size=args.window_size, chr_sizes=chr_sizes, element=args.element,
                                           adj_full_cnt=args.circa_element_adj_cnt_full)
        segments = []
        for segment, cnt in circa_adj_cnts.items():
            segment.extra[args.element + "_cnt"] = cnt
            segments.append(segment)
        write_segments_to_circa_destination(destination=args.output, segments=segments, extra=[args.element + "_cnt"])
    logger.info("Success")


if __name__ == "__main__":
    main()
