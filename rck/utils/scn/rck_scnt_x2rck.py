import argparse
import sys
import os

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, write_scnt_to_destination, get_full_path
from rck.utils.scn.convert import get_scnt_from_battenberg_source, get_scnt_from_hatchet_source, hatchet_get_clone_ids_from_file, \
    get_scnt_from_remixt_source, titan_get_clone_ids_from_file, get_scnt_from_titan_source
from rck.utils.adj.process import get_chromosome_strip_parser


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SCNT-x2rck")
    cli_logging_parser = get_logging_cli_parser()
    chr_strip_parser = get_chromosome_strip_parser()
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    ####
    titan_parser = subparsers.add_parser("titan", parents=[cli_logging_parser, chr_strip_parser])
    titan_parser.add_argument("titan_ichor_seg")
    titan_parser.add_argument("--sample-name", required=True)
    titan_parser.add_argument("--clone-ids", default=None)
    titan_parser.add_argument("--separator", default="\t")
    titan_parser.add_argument("--corrected-cn-fix", choices=["None", "equal", "relative-dist"], default="None")
    titan_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ####
    battenberg_parser = subparsers.add_parser("battenberg", parents=[cli_logging_parser, chr_strip_parser])
    battenberg_parser.add_argument("battenberg", type=argparse.FileType("rt"), default=sys.stdin)
    battenberg_parser.add_argument("--separator", default="\t")
    battenberg_parser.add_argument("--sample-name", required=True)
    battenberg_parser.add_argument("--clone-ids", choices=["1", "2", "1,2"], default="1,2")
    battenberg_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ####
    hatchet_parser = subparsers.add_parser("hatchet", parents=[cli_logging_parser, chr_strip_parser])
    hatchet_parser.add_argument("hatchet", type=str)
    hatchet_parser.add_argument("--separator", default="\t")
    hatchet_parser.add_argument("--min-usage", type=float, default=0.01)
    hatchet_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    group = hatchet_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample-name", default=None)
    group.add_argument("--clone-ids", default=None)
    ####
    remixt_parser = subparsers.add_parser("remixt", parents=[cli_logging_parser, chr_strip_parser])
    remixt_parser.add_argument("remixt", type=argparse.FileType("rt"), default=sys.stdin)
    remixt_parser.add_argument("--separator", default="\t")
    remixt_parser.add_argument("--clone-ids", choices=["1", "2", "1,2"], default="1,2")
    remixt_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ####
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-SCNT")

    if args.command == "titan":
        logger.info("Converting allele-specific segment copy values form TitanCNA format to RCK")
        titan_full_path = get_full_path(path=args.titan_ichor_seg)
        if args.clone_ids is None:
            logger.debug("Clone ids were not provided, extracting all clone ids from {file}".format(file=titan_full_path))
            clone_ids = titan_get_clone_ids_from_file(file_name=titan_full_path, sample_name=args.sample_name, separator=args.separator)
        else:
            clone_ids = sorted(set(args.clone_ids.split(",")))
        logger.debug("Clone ids are identified as {clone_ids}".format(clone_ids=",".join(clone_ids)))
        with open(args.titan_ichor_seg, "rt") as source:
            logger.info("Reading allele-specific segment copy number values from {file}".format(file=titan_full_path))
            segments, scnt = get_scnt_from_titan_source(source=source, sample_name=args.sample_name, clone_ids=clone_ids, separator=args.separator,
                                                        corrected_cn_fix=args.corrected_cn_fix, chr_strip=args.strip_chr)
            logger.info("Writing allele-specific segment copy number values in RCK format to {file}".format(file=args.output))
            write_scnt_to_destination(destination=args.output, segments=segments, scnt=scnt, clone_ids=clone_ids, separator=args.separator)
    elif args.command == "battenberg":
        logger.info("Converting allele-specific segment copy values form Battenberg format to RCK")
        clone_ids = args.clone_ids.split(",")
        logger.debug("Clone ids are identified as {clone_ids}".format(clone_ids=",".join(clone_ids)))
        logger.info("Reading allele-specific segment copy number values form {file}".format(file=args.battenberg))
        segments, scnt = get_scnt_from_battenberg_source(source=args.battenberg, sample_name=args.sample_name, separator=args.separator, chr_strip=args.strip_chr)
        logger.info("Writing allele-specific segment copy number values in RCK format to {file}".format(file=args.output))
        write_scnt_to_destination(destination=args.output, segments=segments, scnt=scnt, separator=args.separator, clone_ids=clone_ids)
    elif args.command == "hatchet":
        hatchet_full_path = get_full_path(path=args.hatchet)
        logger.info("Converting allele-specific segment copy values form HATCHet format to RCK")
        if args.clone_ids is None:
            logger.debug("Clone ids were not provided, extracting all clone ids from {file}".format(file=hatchet_parser))
            clone_ids = hatchet_get_clone_ids_from_file(file_name=hatchet_full_path, sample_name=args.sample_name, separator=args.separator, min_usage=args.min_usage)
        else:
            clone_ids = sorted(set(args.clone_ids.split(",")))
        logger.debug("Clone ids were identified as {clone_ids}".format(clone_ids=",".join(clone_ids)))
        with open(hatchet_full_path) as source:
            logger.info("Reading allele-specific segment copy number values from {file}".format(file=hatchet_full_path))
            segments, scnt = get_scnt_from_hatchet_source(source=source, clone_ids=clone_ids, separator=args.separator, chr_strip=args.strip_chr)
            logger.info("Writing allele-specific segment copy number values in RCK format to {file}".format(file=args.output))
            write_scnt_to_destination(destination=args.output, segments=segments, scnt=scnt, clone_ids=clone_ids, separator=args.separator)
    elif args.command == "remixt":
        logger.info("Converting allele-specific segment copy values form ReMixT format to RCK")
        clone_ids = args.clone_ids.split(",")
        logger.debug("Clone ids were identified as {clone_ids}".format(clone_ids=",".join(clone_ids)))
        logger.info("Reading allele-specific segment copy number values from {file}".format(file=args.remixt))
        segments, scnt = get_scnt_from_remixt_source(source=args.remixt, separator=args.separator, chr_strip=args.strip_chr)
        logger.info("Writing allele-specific segment copy number values in RCK format to {file}".format(file=args.output))
        write_scnt_to_destination(destination=args.output, segments=segments, scnt=scnt, separator=args.separator, clone_ids=clone_ids)
    logger.info("Success!")


if __name__ == "__main__":
    main()
