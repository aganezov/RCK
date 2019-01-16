import argparse
import itertools
import os
import sys
from collections import defaultdict

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import read_adjacencies_from_source, write_adjacencies_to_destination, EXTERNAL_NA_ID, stream_adjacencies_from_source, get_logging_cli_parser, \
    get_standard_logger_from_args
from rck.utils.adj.process import get_shared_nas_parser, Merger, get_chromosome_regions_dict, iter_over_string_entries_from_source, get_extra_field_regexes, filter_adjacencies, \
    KEEP, REMOVE, refined_adjacencies_reciprocal


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADJ-process")
    parser.add_argument('--version', action='version', version=rck.version)
    ####
    shared_parser = get_shared_nas_parser()
    cli_logging_parser = get_logging_cli_parser()
    shared_parser.add_argument("--output", "-o", dest="rck_adj_file", type=argparse.FileType("wt"), default=sys.stdout)
    shared_parser.add_argument("--no-sort", action="store_false", dest="sort")
    ####
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ####
    filter_parser = subparsers.add_parser("filter", parents=[shared_parser, cli_logging_parser], usage=argparse.SUPPRESS)
    filter_parser.add_argument("rck_adj", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    filter_parser.add_argument("--keep-chrs-regions", nargs="+", default=None)
    filter_parser.add_argument("--keep-chrs-regions-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--remove-chrs-regions", nargs="+", default=None)
    filter_parser.add_argument("--remove-chrs-regions-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--keep-extra-field-regex", nargs="+", default=None)
    filter_parser.add_argument("--keep-extra-field-regex-file", type=argparse.FileType, default=None)
    filter_parser.add_argument("--keep-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--remove-extra-field-regex", nargs="+", default=None)
    filter_parser.add_argument("--remove-extra-field-regex-file", type=argparse.FileType, default=None)
    filter_parser.add_argument("--remove-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--merged-field", default="origin_ids")
    filter_parser.add_argument("--merged-distinct-origin-min-cnt", type=int, default=-1)
    filter_parser.add_argument("--merged-distinct-origin-min-cnt-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--merged-origin-regex", default=".*_(?P<source>.*)")
    ####
    cat_parser = subparsers.add_parser("cat", parents=[shared_parser, cli_logging_parser], help="Concatenate Adjacencies in input files (NOTE: different from \"merge\")")
    cat_parser.add_argument("rck_adj", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    cat_parser.add_argument("--enforce-unique-ids", action="store_true", dest="enforce_unique_ids")
    cat_parser.add_argument("--id-collision-strategy", choices=['skip', 'error'], default='error')
    ####
    reciprocal_parser = subparsers.add_parser("reciprocal", parents=[shared_parser, cli_logging_parser], help="ensure that reciprocal novel adjacencies are treated as such")
    reciprocal_parser.add_argument("rck_adj", type=argparse.FileType("rt"), default=sys.stdin)
    reciprocal_parser.add_argument("--max-distance", type=int, default=50)
    ####
    merge = subparsers.add_parser("merge", parents=[shared_parser, cli_logging_parser], help="EXPERIMENTAL! Merge Adjacencies in input files (NOTE: different from \"cat\").")
    merge.add_argument("rck_adj", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    merge.add_argument("--max-distance", type=int, default=500)
    merge.add_argument("--enforce-sv-types", action="store_true", dest="enforce_sv_types")
    merge.add_argument("--merged-field", default="origin_ids")
    merge.add_argument("--merged-id-template", default="{cnt}_merged")
    ####
    test_merge = subparsers.add_parser("test-merge", parents=[shared_parser, cli_logging_parser], usage=argparse.SUPPRESS)
    test_merge.add_argument("sources", type=argparse.FileType("rt"), nargs="+")
    test_merge.add_argument("--merged", type=argparse.FileType("rt"), required=True)
    test_merge.add_argument("--max-distance", default=500)
    test_merge.add_argument("--merged-field", default="origin_ids")
    test_merge.add_argument("--origin-sep", default=",")
    test_merge.add_argument("--origin-regex", default=".*_(?P<source>.*)")
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-ADK-process")
    processed_adjacencies = []
    if args.o_extra_fields is None or len(args.o_extra_fields) == 0 or args.o_extra_fields == ",":
        extra = None
    elif args.o_extra_fields != "all":
        extra = args.o_extra_fields.split(",")
    else:
        extra = args.o_extra_fields
    if args.command == "merge":
        all_adjacencies = []
        for rck_adj_source in args.rck_adj:
            all_adjacencies.extend(read_adjacencies_from_source(source=rck_adj_source))
        merger = Merger(origin_ids_field=args.merged_field)
        for na in all_adjacencies:
            merger.add_adjacency(adjacency=na, max_distance=args.max_distance)
        processed_adjacencies = merger.get_merged_adjacencies(merged_template=args.merged_id_template)
    elif args.command == "test-merge":
        source_adjacencies = []
        for rck_adj_source in args.sources:
            source_adjacencies.extend(read_adjacencies_from_source(source=rck_adj_source))
        merged_nas = read_adjacencies_from_source(source=args.merged)
        source_nas_by_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.idx): adj for adj in source_adjacencies}
        logger.info("Total number of merged nas :: {mlen}".format(mlen=len(merged_nas)))
        logger.info("Total number of source nas :: {slen}".format(slen=len(source_adjacencies)))
        source_nas_groups = defaultdict(list)
        cnt = 0
        for merged_na in merged_nas:
            if args.merged_field not in merged_na.extra:
                logger.info("no source ids field {{{field}}} in merged_na {naid}".format(field=args.merged_field, naid=merged_na.extra.get(EXTERNAL_NA_ID, merged_na.idx)))
                continue
            source_ids = merged_na.extra[args.merged_field].split(args.origin_sep)
            for source_id in source_ids:
                source_nas_groups[source_id].append(merged_na.extra.get(EXTERNAL_NA_ID, merged_na.idx))
            for ona1_id, ona2_id in itertools.combinations(source_ids, 2):
                ona1 = source_nas_by_ids[ona1_id]
                ona2 = source_nas_by_ids[ona2_id]
                if not Merger.adjacencies_mergeable(adjacency1=ona1, adjacency2=ona2, max_distance=args.max_distance):
                    cnt += 1
        logger.info("Non mergable pairs cnt :: {cnt}".format(cnt=cnt))
        for source_id, groups in source_nas_groups.items():
            if len(groups) > 1:
                logger.info("source na id {source_id} is present in several merged adjacencies {merged_ids}".format(source_id=source_id, merged_ids=str(groups)))
    elif args.command == "cat":
        adjacencies = itertools.chain(*(stream_adjacencies_from_source(source=rck_adj_source) for rck_adj_source in args.rck_adj))
        if args.enforce_unique_ids:
            processed_ids = set()
            result = []
            for adj in adjacencies:
                aid = adj.extra.get(EXTERNAL_NA_ID, adj.idx)
                if aid in processed_ids:
                    logger.debug("Adjacency id {aid} has been encountered more than once".format(aid=aid))
                    if args.id_collision_strategy == "skip":
                        continue
                    elif args.id_collision_strategy == "error":
                        raise ValueError("More than one adjacency with id {aid}".format(aid=aid))
                result.append(adj)
                processed_ids.add(aid)
            adjacencies = result
        processed_adjacencies = list(adjacencies)
    elif args.command == "filter":
        adjacencies = itertools.chain(*(stream_adjacencies_from_source(source=rck_adj_source) for rck_adj_source in args.rck_adj))
        keep_chromosome_regions_entries = args.keep_chrs_regions if args.keep_chrs_regions is not None else []
        if args.keep_chrs_regions_file is not None:
            keep_chromosome_regions_entries.extend(list(iter_over_string_entries_from_source(source=args.keep_chrs_regions_file)))
        remove_chr_regions_entries = args.remove_chrs_regions if args.remove_chrs_regions is not None else []
        if args.remove_chrs_regions_file is not None:
            remove_chr_regions_entries.extend(list(iter_over_string_entries_from_source(source=args.remove_chrs_regions_file)))
        keep_regions = get_chromosome_regions_dict(string_entries=keep_chromosome_regions_entries)
        remove_regions = get_chromosome_regions_dict(string_entries=remove_chr_regions_entries)
        keep_extra_field_entries = args.keep_extra_field_regex if args.keep_extra_field_regex is not None else []
        if args.keep_extra_field_regex_file is not None:
            keep_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.keep_extra_field_regex_file)))
        remove_extra_field_entries = args.remove_extra_field_regex if args.remove_extra_field_regex is not None else []
        if args.remove_extra_field_regex_file is not None:
            remove_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.remove_extra_field_regex_file)))
        keep_extra_field = get_extra_field_regexes(string_entries=keep_extra_field_entries)
        remove_extra_field = get_extra_field_regexes(string_entries=remove_extra_field_entries)
        result = filter_adjacencies(adjacencies=adjacencies, keep_regions=keep_regions, remove_regions=remove_regions,
                                    keep_extra_field=keep_extra_field, keep_extra_field_missing_strategy=args.keep_extra_field_missing_strategy,
                                    remove_extra_field=remove_extra_field, remove_extra_field_missing_strategy=args.remove_extra_field_missing_strategy,
                                    merged_field=args.merged_field,
                                    merged_distinct_origin_cnt=args.merged_distinct_origin_min_cnt,
                                    merged_origin_regex_string=args.merged_origin_regex,
                                    merged_distinct_origin_min_cnt_missing_strategy=args.merged_distinct_origin_min_cnt_missing_strategy)
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=result, sort_adjacencies=False, extra=extra)
        exit(0)
    elif args.command == "reciprocal":
        adjacencies = read_adjacencies_from_source(source=args.rck_adj)
        processed_adjacencies = refined_adjacencies_reciprocal(novel_adjacencies=adjacencies, max_distance=args.max_distance, inplace=True)

    if len(processed_adjacencies) > 0:
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=processed_adjacencies, extra=extra, sort_adjacencies=args.sort)


if __name__ == "__main__":
    main()
