import argparse
import sys
import os

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from rck.utils.scn.process import get_haploid_scnt, get_circa_segments_cna_fractions
from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, read_scnt_from_source, write_scnt_to_shatterseek_destination, read_chr_sizes_from_source, \
    write_segments_to_circa_destination
from rck.utils.adj.process import get_chromosome_strip_parser


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SCNT-rck2x")
    cli_logging_parser = get_logging_cli_parser()
    chr_strip_parser = get_chromosome_strip_parser()
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    ####
    shatterseek_parser = subparsers.add_parser("shatterseek", parents=[cli_logging_parser, chr_strip_parser])
    shatterseek_parser.add_argument("rck_scnt", type=argparse.FileType("rt"), default=sys.stdin)
    shatterseek_parser.add_argument("--clone-id", required=True)
    shatterseek_parser.add_argument("--separator", default="\t")
    shatterseek_parser.add_argument("--extra-separator", default=";")
    shatterseek_parser.add_argument("--default-cn", type=int, default=0)
    shatterseek_parser.add_argument("--output-header", action="store_true", dest="output_header")
    shatterseek_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ####
    circa_dens_parser = subparsers.add_parser("circa-dens", parents=[cli_logging_parser, chr_strip_parser])
    circa_dens_parser.add_argument("rck_scnt", type=argparse.FileType("rt"), default=sys.stdin)
    circa_dens_parser.add_argument("--clone-id", required=True)
    circa_dens_parser.add_argument("--separator", default="\t")
    circa_dens_parser.add_argument("--extra-separator", default=";")
    circa_dens_parser.add_argument("--cna-type", choices=["ampl", "del"], default="ampl")
    circa_dens_parser.add_argument("--haploid", action="store_true", dest="haploid")
    circa_dens_parser.add_argument("--inverse", action="store_true", dest="inverse")
    circa_dens_parser.add_argument("--window-size", type=int, default=10000000)
    circa_dens_parser.add_argument("--chr-sizes", type=argparse.FileType("rt"))
    circa_dens_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ####
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-SCNT")

    if args.command == "shatterseek":
        logger.info("Starting converting RCK Segment Copy Number Tensor data to ShatterSeek")
        logger.debug("Specified clone is {clone_id}".format(clone_id=args.clone_id))
        logger.info("Reading RCK formatted data from {file}".format(file=args.rck_scnt))
        segments, scnt = read_scnt_from_source(source=args.rck_scnt, separator=args.separator, extra_separator=args.extra_separator)
        logger.info("Read CN data is translated into a haploid (!!!) version of itself.")
        haploid_scnt = get_haploid_scnt(segments=segments, scnt=scnt)
        logger.info("Writing data for clone {clone_id} in a ShatterSeek suitable format to {file}".format(clone_id=args.clone_id, file=args.output))
        write_scnt_to_shatterseek_destination(destination=args.output, segments=segments, scnt=haploid_scnt, clone_id=args.clone_id,
                                              default=args.default_cn, output_header=args.output_header)
    elif args.command == "circa-dens":
        logger.info("Starting computing ampl/del statistics from RKC Segment Copy Number Tensor Format")
        logger.debug("Specified clone is {clone_id}".format(clone_id=args.clone_id))
        logger.info("Reading RCK formatted data from {file}".format(file=args.rck_scnt))
        segments, scnt = read_scnt_from_source(source=args.rck_scnt, separator=args.separator, extra_separator=args.extra_separator)
        chr_sizes = args.chr_sizes
        if args.chr_sizes is not None:
            chr_sizes = read_chr_sizes_from_source(source=args.chr_sizes)
        circa_segments_cna_fractions = get_circa_segments_cna_fractions(segments=segments, scnt=scnt, clone_id=args.clone_id,
                                                                        window_size=args.window_size, chr_sizes=chr_sizes, cna_type=args.cna_type,
                                                                        haploid=args.haploid)
        segments = []
        total_average = 0
        total_length = 0
        for segment, cna_fraction in circa_segments_cna_fractions.items():
            value = cna_fraction * segment.length / args.window_size
            if args.inverse:
                value = 1 - value
            segment.extra[args.cna_type + "_fraction"] = value
            total_length += segment.length
            total_average += cna_fraction * segment.length
            segments.append(segment)
        print("total average cna fraction", total_average / total_length)
        write_segments_to_circa_destination(destination=args.output, segments=segments, extra=[args.cna_type + "_fraction"])
    logger.info("Success!")


if __name__ == "__main__":
    main()
