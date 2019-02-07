import argparse
import itertools

import os
import sys

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, stream_adjacency_groups_from_source, write_adjacency_groups_to_destination, \
    read_adjacency_groups_from_source
from rck.core.structures import AdjacencyGroupType
from rck.utils.adj.adjacency_group_process import refined_labeling_groups


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADJ-GROUPS-process")
    parser.add_argument("--version", action="version", version=rck.version)
    cli_logging_parser = get_logging_cli_parser()

    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    ###
    cat_parser = subparsers.add_parser("cat", parents=[cli_logging_parser])
    cat_parser.add_argument("rck_adg", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    cat_parser.add_argument("--i-separator", default="\t")
    cat_parser.add_argument("--i-extra-separator", default=";")
    cat_parser.add_argument("--i-aids-separator", default=",")
    cat_parser.add_argument("--enforce-unique-ids", action="store_true", dest="enforce_unique_ids")
    cat_parser.add_argument("--id-collision-strategy", choices=["skip", "error"], default="error")
    cat_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    cat_parser.add_argument("--o-separator", default="\t")
    cat_parser.add_argument("--o-aids-separator", default=",")
    cat_parser.add_argument("--o-extra-separator", default=";")
    ###
    refine_parser = subparsers.add_parser("refine", parents=[cli_logging_parser])
    refine_parser.add_argument("rck_adg", type=argparse.FileType("rt"), default=sys.stdin)
    refine_parser.add_argument("--i-separator", default="\t")
    refine_parser.add_argument("--i-extra-separator", default=";")
    refine_parser.add_argument("--i-aids-separator", default=",")
    # refine_parser.add_argument("--no-refine-m", action="store_false", dest="refine_m")
    # refine_parser.add_argument("--no-refine-l", action="store_false", dest="refine_l")
    # refine_parser.add_argument("--no-refine-n", action="store_false", dest="refine_n")
    refine_parser.add_argument("--gid-suffix", default="refined")
    refine_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    refine_parser.add_argument("--o-separator", default="\t")
    refine_parser.add_argument("--o-aids-separator", default=",")
    refine_parser.add_argument("--o-extra-separator", default=";")
    ###
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-ADJ-GROUPS-process")
    if args.command == "cat":
        adj_groups = itertools.chain(*(stream_adjacency_groups_from_source(source=adj_group_source, separator=args.i_separator,
                                                                           aids_separator=args.i_aids_separator, extra_separator=args.i_extra_separator)
                                       for adj_group_source in args.rck_adg))
        if args.enforce_unique_ids:
            pass
        write_adjacency_groups_to_destination(destination=args.output, adjacency_groups=adj_groups, separator=args.o_separator,
                                              aids_separator=args.o_aids_separator, extra_separator=args.o_extra_separator)
    elif args.command == "refine":
        logger.info("Refining input adjacency groups")
        logger.info("Reading adjacency groups from {file}".format(file=args.rck_adg))
        adg_groups = read_adjacency_groups_from_source(source=args.rck_adg, separator=args.i_separator,
                                                       extra_separator=args.i_extra_separator, aids_separator=args.i_aids_separator)
        logger.info("A total of {cnt} adjacency groups has been read".format(cnt=len(adg_groups)))
        molecule_groups = [ag for ag in adg_groups if ag.group_type == AdjacencyGroupType.MOLECULE]
        logger.info("A total of {cnt} molecule adjacency groups has been read".format(cnt=len(molecule_groups)))
        labeling_groups = [ag for ag in adg_groups if ag.group_type == AdjacencyGroupType.LABELING]
        logger.info("A total of {cnt} labeling adjacency groups has been read".format(cnt=len(labeling_groups)))
        general_groups = [ag for ag in adg_groups if ag.group_type == AdjacencyGroupType.GENERAL]
        logger.info("A total of {cnt} general adjacency groups has been read".format(cnt=len(general_groups)))
        logger.info("Refining molecule adjacency groups")
        refined_molecule_groups = molecule_groups
        logger.info("A total of {cnt} refined molecule adjacency groups remains".format(cnt=len(refined_molecule_groups)))
        logger.info("Refining labeling adjacency groups")
        refined_labeling_groups = refined_labeling_groups(adj_groups=labeling_groups, gid_suffix="" if len(args.gid_suffix) == 0 else args.gid_suffix + "-L",
                                                          retain_source_gids=True)
        logger.info("A total of {cnt} refined labeling adjacency groups remains".format(cnt=len(refined_labeling_groups)))
        logger.info("Refining general adjacency groups")
        refined_general_groups = general_groups
        logger.info("A total of {cnt} refined labeling general adjacency groups remains".format(cnt=len(refined_general_groups)))
        adj_groups = itertools.chain(refined_molecule_groups, refined_labeling_groups, refined_general_groups)
        logger.info("Writing refined adjacency groups to {file}".format(file=args.output))
        write_adjacency_groups_to_destination(destination=args.output, adjacency_groups=adj_groups, separator=args.o_separator, aids_separator=args.o_aids_separator)


if __name__ == "__main__":
    main()
