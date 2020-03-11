import argparse
import sys

from rck.core.io import get_logging_cli_parser, read_scnt_from_source
from rck.utils.scn.stats import cn_distance


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SCNT-stats")
    cli_logging_parser = get_logging_cli_parser()
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    #####
    distance_parser = subparsers.add_parser("distance", parents=[cli_logging_parser])
    distance_parser.add_argument('--scnt1', type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt1-separator", default="\t")
    distance_parser.add_argument("--scnt1-extra-separator", default=";")
    distance_parser.add_argument("--scnt1-clone-ids", default=None)
    distance_parser.add_argument('--scnt2', type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt2-separator", default="\t")
    distance_parser.add_argument("--scnt2-extra-separator", default=";")
    distance_parser.add_argument("--scnt2-clone-ids", default=None)
    distance_parser.add_argument("--topn", type=int, default=3)
    distance_parser.add_argument("--verbose", action="store_true", dest="verbose")
    distance_parser.add_argument("--both-haplotype-specific", action="store_true", dest="both_haplotype_specific")
    distance_parser.add_argument('-o', '--output', type=argparse.FileType("wt"), default=sys.stdout)
    #####
    args = parser.parse_args()
    if args.command == "distance":
        scnt1_clone_ids = args.scnt1_clone_ids if args.scnt1_clone_ids is None else args.scnt1_clone_ids.split(",")
        segments1, scnt1 = read_scnt_from_source(source=args.scnt1, separator=args.scnt1_separator, extra_separator=args.scnt1_extra_separator, clone_ids=scnt1_clone_ids)
        scnt2_clone_ids = args.scnt2_clone_ids if args.scnt2_clone_ids is None else args.scnt2_clone_ids.split(",")
        segments2, scnt2 = read_scnt_from_source(source=args.scnt2, separator=args.scnt2_separator, extra_separator=args.scnt2_extra_separator, clone_ids=scnt2_clone_ids)
        result = cn_distance(segments1=segments1, scnt1=scnt1, segments2=segments2, scnt2=scnt2, both_haplotype_specific=args.both_haplotype_specific)
        sorted_result = sorted([(key, value) for key, value in result.items()], key=lambda entry: sum(entry[1].values()))
        output_result = sorted_result[:args.topn]
        if args.verbose:
            print(f'Length-weighted segment copy number distance for tensors in {args.scnt1.name} and {args.scnt2.name}', file=args.output)
        for cnt, (case, clone_specific_distance) in enumerate(output_result, start=1):
            print(f'{cnt}. Best distance (total) of {sum(clone_specific_distance.values()):,} with clone-specific ones {clone_specific_distance}, for case {case}', file=args.output)


if __name__ == "__main__":
    main()
