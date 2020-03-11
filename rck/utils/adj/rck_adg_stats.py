import argparse
import sys

import rck
from rck.core.io import get_logging_cli_parser, read_adjacency_groups_from_source
from rck.core.structures import AdjacencyGroupType
from rck.utils.adj.adjacency_group_stats import groups_size_tally


def main():
    cli_logging_parser = get_logging_cli_parser()
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADG-STATS")
    parser.add_argument('--version', action='version', version=rck.version)
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    #######
    labeling_group_size_parser = subparsers.add_parser("size-l", parents=[cli_logging_parser], help="Group size for RCK AdjGROUP in input file")
    labeling_group_size_parser.add_argument("rck_adg", type=argparse.FileType("rt"), nargs="?", default=sys.stdin)
    labeling_group_size_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    labeling_group_size_parser.add_argument("--no-allow-zero-values", action="store_false", dest="allow_zero_values")
    labeling_group_size_parser.add_argument("--min", type=int, default=-1)
    labeling_group_size_parser.add_argument("--max", type=int, default=-1)
    #######
    args = parser.parse_args()
    if args.command == "size-l":
        adj_groups = read_adjacency_groups_from_source(source=args.rck_adg)
        labeling_adg = [ag for ag in adj_groups if ag.group_type == AdjacencyGroupType.LABELING]
        tally = groups_size_tally(adjacency_groups=labeling_adg)
        min_key, max_key = min(tally.keys()), max(tally.keys())
        if args.max != -1:
            max_key = args.max
        if args.min != -1:
            min_key = args.min
        min_value = 0
        for key in tally:
            if key < min_key:
                min_value += tally[key]
        print("<{min_key}".format(min_key=min_key), min_value, sep=",", file=args.output)
        for value in range(min_key, max_key):
            if value not in tally and not args.allow_zero_values:
                continue
            print(value, tally.get(value, 0), sep=",", file=args.output)
        max_value = 0
        for key in tally:
            if key >= max_key:
                max_value += tally[key]
        print(">={max_key}".format(max_key=max_key), max_value, sep=",", file=args.output)



if __name__ == "__main__":
    main()
