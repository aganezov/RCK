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
from rck.utils.adj.process import get_shared_nas_parser, Merger, iter_over_string_entries_from_source, get_extra_field_regexes, \
    filter_adjacencies_by_extra, \
    KEEP, REMOVE, refined_adjacencies_reciprocal, update_adjacencies
from rck.utils.adj.convert import get_chrs_regions_string_lists_from_source, get_chrs_regions_string_list_from_file, parse_segment_chr_region
from rck.utils.adj.process import filter_adjacencies_by_chromosomal_regions, filter_adjacencies_by_size, iter_haploid_adjacencies


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
    filter_parser = subparsers.add_parser("filter", parents=[shared_parser, cli_logging_parser])
    filter_parser.add_argument("rck_adj", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    filter_parser.add_argument("--keep-extra-field-regex", action="append", default=None)
    filter_parser.add_argument("--keep-extra-field-regex-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--keep-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--keep-annotate", action="store_true", dest="annotate_retained")
    filter_parser.add_argument("--keep-annotate-s-extra-field", default=None, dest="annotate_seg_extra_field")
    filter_parser.add_argument("--keep-annotate-short-circ", action="store_true", dest="annotate_shirt_circ")
    filter_parser.add_argument("--keep-annotate-extra-prefix", dest="annotate_extra_prefix")
    filter_parser.add_argument("--remove-extra-field-regex", action="append", default=None)
    filter_parser.add_argument("--remove-extra-field-regex-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--remove-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--min-size", type=int, default=0)
    filter_parser.add_argument("--max-size", type=int, default=1000000000)
    filter_parser.add_argument("--no-allow-inter-chr", action="store_false", dest="allow_inter_chr")
    filter_parser.add_argument("--no-allow-intra-chr", action="store_false", dest="allow_intra_chr")
    filter_parser.add_argument("--size-extra-field", default="svlen")
    filter_parser.add_argument("--size-extra-field-no-abs", action="store_false", dest="size_extra_field_abs")
    filter_parser.add_argument("--size-extra-seq-field")
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
    haploid_parser = subparsers.add_parser("haploid", parents=[shared_parser, cli_logging_parser], help="collapse any info that is allele/haplotype-specific into a haploid mode")
    haploid_parser.add_argument("rck_adj", type=argparse.FileType("rt"), nargs="+", default=[sys.stdin])
    ####
    update_parser = subparsers.add_parser("update", parents=[shared_parser, cli_logging_parser],
                                          help="Updates adjacencies in the 'adj' with the info from --source based on aid matches. Outputs updated --target entries")
    update_parser.add_argument("rck_adj", type=argparse.FileType("rt"))
    update_parser.add_argument("--source", type=argparse.FileType("rt"), required=True)
    update_parser.add_argument("--exclude-extra-fields", default="")
    update_parser.add_argument("--include-extra-fields", default="")
    update_parser.add_argument("--no-include-missing", action="store_false", dest="include_missing")
    update_parser.add_argument("--no-coords-update", action="store_false", dest="coord_update")
    update_parser.add_argument("--no-coord1-update", action="store_false", dest="coord1_update")
    update_parser.add_argument("--no-coord2-update", action="store_false", dest="coord2_update")
    update_parser.add_argument("--no-strands-update", action="store_false", dest="strands_update")
    update_parser.add_argument("--no-strand1-update", action="store_false", dest="strand1_update")
    update_parser.add_argument("--no-strand2-update", action="store_false", dest="strand2_update")
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-ADK-process")
    processed_adjacencies = []
    if args.o_extra_fields is None or len(args.o_extra_fields) == 0 or args.o_extra_fields == ",":
        extra = None
    elif args.o_extra_fields != "all":
        extra = args.o_extra_fields.split(",")
    else:
        extra = args.o_extra_fields
    if args.command == "cat":
        adjacencies = itertools.chain(*(stream_adjacencies_from_source(source=rck_adj_source) for rck_adj_source in args.rck_adj))
        if args.enforce_unique_ids:
            processed_ids = set()
            adjacencies = []
            for adj in adjacencies:
                aid = adj.extra.get(EXTERNAL_NA_ID, adj.idx)
                if aid in processed_ids:
                    logger.debug("Adjacency id {aid} has been encountered more than once".format(aid=aid))
                    if args.id_collision_strategy == "skip":
                        continue
                    elif args.id_collision_strategy == "error":
                        raise ValueError("More than one adjacency with id {aid}".format(aid=aid))
                adjacencies.append(adj)
                processed_ids.add(aid)
            adjacencies = adjacencies
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=adjacencies, extra=extra, sort_adjacencies=args.sort)
        exit(0)
    elif args.command == "filter":
        logger.info("Filtering input adjacencies from following sources {sources}".format(sources=",".join(map(str, args.rck_adj))))
        adjacencies = itertools.chain(*(stream_adjacencies_from_source(source=rck_adj_source) for rck_adj_source in args.rck_adj))
        include_chrs_regions_strings = []
        exclude_chrs_regions_strings = []
        if args.chrs_include is not None:
            for chrs_lists in args.chrs_include:
                for chrs_list in chrs_lists:
                    for chr_name in chrs_list.split(","):
                        include_chrs_regions_strings.append(chr_name)
        if args.chrs_include_file is not None:
            for chr_name in get_chrs_regions_string_lists_from_source(source=args.chrs_include_file):
                include_chrs_regions_strings.append(chr_name)
        if args.chrs_exclude is not None:
            for chrs_lists in args.chrs_exclude:
                for chrs_list in chrs_lists:
                    for chr_name in chrs_list.split(","):
                        exclude_chrs_regions_strings.append(chr_name)
        if args.chrs_exclude_file is not None:
            for chr_name in get_chrs_regions_string_list_from_file(file_name=args.chrs_exclude_file):
                exclude_chrs_regions_strings.append(chr_name)
        include_regions = [parse_segment_chr_region(string) for string in include_chrs_regions_strings]
        exclude_regions = [parse_segment_chr_region(string) for string in exclude_chrs_regions_strings]
        adjacencies = filter_adjacencies_by_chromosomal_regions(adjacencies=adjacencies, include=include_regions, exclude=exclude_regions,
                                                                include_both=args.include_both, exclude_both=args.exclude_both,
                                                                include_spanning=args.include_spanning, exclude_spanning=args.exclude_spanning,
                                                                annotate_retained=args.annotate_retained, annotate_retained_extra_field_prefix=args.annotate_extra_prefix,
                                                                annotated_retained_segments_extra_field=args.annotate_seg_extra_field, annotate_short_circ=args.annotate_shirt_circ)
        keep_extra_field_entries = args.keep_extra_field_regex if args.keep_extra_field_regex is not None else []
        if args.keep_extra_field_regex_file is not None:
            keep_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.keep_extra_field_regex_file)))
        remove_extra_field_entries = args.remove_extra_field_regex if args.remove_extra_field_regex is not None else []
        if args.remove_extra_field_regex_file is not None:
            remove_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.remove_extra_field_regex_file)))
        keep_extra_field = get_extra_field_regexes(string_entries=keep_extra_field_entries)
        remove_extra_field = get_extra_field_regexes(string_entries=remove_extra_field_entries)
        adjacencies = filter_adjacencies_by_extra(adjacencies=adjacencies,
                                                  keep_extra_field=keep_extra_field, keep_extra_field_missing_strategy=args.keep_extra_field_missing_strategy,
                                                  remove_extra_field=remove_extra_field, remove_extra_field_missing_strategy=args.remove_extra_field_missing_strategy)
        adjacencies = filter_adjacencies_by_size(adjacencies=adjacencies, min_size=args.min_size, max_size=args.max_size, size_extra_field=args.size_extra_field,
                                                 size_extra_seq_field=args.size_extra_seq_field, allow_inter_chr=args.allow_inter_chr,
                                                 size_extra_field_abs=args.size_extra_field_abs, allow_intra_chr=args.allow_intra_chr,)
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=adjacencies, sort_adjacencies=False, extra=extra)
        exit(0)
    elif args.command == "reciprocal":
        adjacencies = read_adjacencies_from_source(source=args.rck_adj)
        processed_adjacencies = refined_adjacencies_reciprocal(novel_adjacencies=adjacencies, max_distance=args.max_distance, inplace=True)
    elif args.command == "haploid":
        adjacencies = itertools.chain(*(stream_adjacencies_from_source(source=rck_adj_source) for rck_adj_source in args.rck_adj))
        haploid_adjacencies = iter_haploid_adjacencies(adjacencies=adjacencies, copy=False)
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=haploid_adjacencies, sort_adjacencies=False, extra=extra)
        exit(0)
    elif args.command == "update":
        adjacencies = read_adjacencies_from_source(source=args.rck_adj)
        source_adjacencies = read_adjacencies_from_source(source=args.source)
        extra_include = {v for v in args.include_extra_fields.split(",") if len(v) > 0}
        extra_exclude = {v for v in args.exclude_extra_fields.split(",") if len(v) > 0}
        processed_adjacencies = update_adjacencies(target_adjacencies=adjacencies, source_adjacencies=source_adjacencies,
                                                   update_coords=args.update_coords, update_coord1=args.update_coord1, update_coord2=args.update_coord2,
                                                   update_strands=args.update_strands, update_strand1=args.update_strand1, update_strand2=args.update_strand2,
                                                   extra_exclude=extra_exclude, extra_include=extra_include, include_missing=args.include_missing)
    if len(processed_adjacencies) > 0:
        write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=processed_adjacencies, extra=extra, sort_adjacencies=args.sort)


if __name__ == "__main__":
    main()
