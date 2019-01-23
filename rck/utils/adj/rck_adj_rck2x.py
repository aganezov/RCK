import argparse
import os
import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, read_adjacencies_from_source, write_adjacencies_to_vcf_sniffles_destination
from rck.core.structures import AdjacencyType


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-NAS-rck2x")
    parser.add_argument('--version', action='version', version=rck.version)
    cli_logging_parser = get_logging_cli_parser()
    ###
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ###
    vcf_parser = subparsers.add_parser("vcf-sniffles", parents=[cli_logging_parser], help="Convert RCK Adjacencies to the VCF (Sniffles) format")
    vcf_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin,)
    vcf_parser.add_argument("--separator", default="\t")
    vcf_parser.add_argument("--extra-separator", default=";")
    vcf_parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    vcf_parser.add_argument("--o-extra-fields", default="all")
    vcf_parser.add_argument("--o-no-include-ref", action="store_false", dest="include_ref")
    vcf_parser.add_argument("--clone-suffix", default="")
    vcf_parser.add_argument("--dummy-clone", default="dummy_clone")
    ###
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-NAS-rck2x")
    if args.o_extra_fields is None or len(args.o_extra_fields) == 0 or args.o_extra_fields == ",":
        extra = None
    elif args.o_extra_fields != "all":
        extra = args.o_extra_fields.split(",")
    else:
        extra = args.o_extra_fields
    logger.debug("Output extra fields are identified as {o_extra}".format(o_extra=",".join(extra) if extra is not None else ""))
    logger.info("Reading adjacencies from {file}".format(file=args.rck_adj))
    adjacencies = read_adjacencies_from_source(source=args.rck_adj, extra_separator=args.extra_separator, separator=args.separator)
    if not args.include_ref:
        logger.debug("Reference adjacencies were excluded from the output. ")
        adjacencies = list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, adjacencies))
    if args.command == "vcf-sniffles":
        logger.info("Converting RCK formatted adjacencies to the VCF (Sniffles) format")
        logger.info("Writing adjacencies to {file}".format(file=args.output))
        write_adjacencies_to_vcf_sniffles_destination(destination=args.output, adjacencies=adjacencies, extra=extra,
                                                      dummy_clone=args.dummy_clone, clone_suffix=args.clone_suffix)


if __name__ == "__main__":
    main()
