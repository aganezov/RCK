import argparse
import itertools
import sys
import os

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, get_full_path, read_scnt_from_file, write_scnt_to_file, \
    write_scnt_to_destination, read_scnt_from_source, stream_segments_from_source, write_segments_to_destination
from rck.core.structures import aligned_scnts, refined_scnt, cn_distance_inter_scnt
from rck.utils.adj.process import KEEP, REMOVE, iter_over_string_entries_from_source, get_extra_field_regexes
from rck.utils.scn.process import iter_haploid_segments, filter_segments_by_chromosomal_regions, filter_segments_by_extra, filter_segments_by_size
from rck.utils.adj.convert import get_chrs_regions_string_list_from_file, parse_segment_chr_region, get_chrs_regions_string_lists_from_source


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SCNT-process")
    cli_logging_parser = get_logging_cli_parser()
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    ###
    refine_parser = subparsers.add_parser("refine", parents=[cli_logging_parser])
    refine_parser.add_argument('scnt', type=argparse.FileType("rt"), default=sys.stdin)
    refine_parser.add_argument("--separator", default="\t")
    refine_parser.add_argument("--no-allow-missing-clones", action="store_false", dest="allow_missing_clones")
    refine_parser.add_argument("--clone-ids", default=None)
    refine_parser.add_argument("--no-merge-fragments", action="store_false", dest="merge_fragments")
    refine_parser.add_argument("--max-merge-gap", type=int, default=1000000)
    refine_parser.add_argument("--no-fill-gaps", action="store_false", dest="fill_gaps")
    refine_parser.add_argument("--max-fill-gap", type=int, default=1000000)
    refine_parser.add_argument('--output', type=argparse.FileType("wt"), default=sys.stdout)
    ###
    align_parser = subparsers.add_parser("align", parents=[cli_logging_parser])
    align_parser.add_argument("scnt", nargs="+")
    align_parser.add_argument("--separator", default="\t")
    align_parser.add_argument("--output-suffix", default="aligned")
    align_parser.add_argument("--no-allow-unit-segments", action="store_false", dest="allow_unit_segments")
    align_parser.add_argument("--output-dir", default="")
    ###
    distance_parser = subparsers.add_parser("distance", parents=[cli_logging_parser])
    distance_parser.add_argument("--scnt1", type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt1-separator", default="\t")
    distance_parser.add_argument("--scnt1-extra-separator", default=";")
    distance_parser.add_argument("--scnt2", type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt2-separator", default="\t")
    distance_parser.add_argument("--scnt2-extra-separator", default=";")
    distance_parser.add_argument("--clone-ids", default=None)
    distance_parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    ###
    filter_parser = subparsers.add_parser("filter", parents=[cli_logging_parser])
    filter_parser.add_argument("scnt", type=argparse.FileType("rt"), default=sys.stdin)
    filter_parser.add_argument("--separator", default="\t")
    filter_parser.add_argument("--extra-separator", default=";")
    filter_parser.add_argument("--o-extra-fields", default="all")
    filter_parser.add_argument("--chrs-include", action="append", nargs=1)
    filter_parser.add_argument("--chrs-include-file", type=argparse.FileType("rt"))
    filter_parser.add_argument("--chrs-include-no-full", action="store_false", dest="include_full")
    filter_parser.add_argument("--chrs-exclude", action="append", nargs=1)
    filter_parser.add_argument("--chrs-exclude-file", type=argparse.FileType("rt"))
    filter_parser.add_argument("--chrs-exclude-full", action="store_true", dest="exclude_full")
    filter_parser.add_argument("--keep-extra-field-regex", nargs="+", default=None)
    filter_parser.add_argument("--keep-extra-field-regex-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--keep-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--remove-extra-field-regex", nargs="+", default=None)
    filter_parser.add_argument("--remove-extra-field-regex-file", type=argparse.FileType("rt"), default=None)
    filter_parser.add_argument("--remove-extra-field-missing-strategy", choices=[KEEP, REMOVE], default=KEEP)
    filter_parser.add_argument("--min-size", type=int, default=0)
    filter_parser.add_argument("--max-size", type=int, default=1000000000)
    filter_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ###
    haploid_parser = subparsers.add_parser("haploid", parents=[cli_logging_parser])
    haploid_parser.add_argument("scnt", type=argparse.FileType("rt"), default=sys.stdin)
    haploid_parser.add_argument("--separator", default="\t")
    haploid_parser.add_argument("--extra-separator", default=";")
    haploid_parser.add_argument("--output", "-o", type=argparse.FileType("wt"), default=sys.stdout)
    ###
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-SCNT-process")

    if args.command == "refine":
        clone_ids = args.clone_ids.split(",") if args.clone_ids is not None else None
        logger.debug("Clone ids identified as {clone_ids}. If None -- all clone ids will be processed.".format(clone_ids=",".join(clone_ids)))
        logger.info("Reading Segment Copy Number Tensor form {file}".format(file=args.scnt))
        segments, scnt = read_scnt_from_source(source=args.scnt, clone_ids=clone_ids, separator=args.separator)
        logger.info("Refining Segment Copy Number Tensor from {file}".format(file=args.scnt))
        segments, scnt, _ = refined_scnt(segments=segments, scnt=scnt,
                                         merge_fragments=args.merge_fragments, max_merge_gap=args.max_merge_gap,
                                         fill_gaps=args.fill_gaps, max_fill_gap=args.max_fill_gap)
        logger.info("Writing refined Segment Copy Number Tensor to {file}".format(file=args.output))
        write_scnt_to_destination(destination=args.output, scnt=scnt, segments=segments, clone_ids=clone_ids, separator=args.separator)
    elif args.command == "align":
        scnt_files = {}
        for path in args.scnt:
            full_path = get_full_path(path=path)
            name = os.path.splitext(os.path.basename(full_path))[0]
            if name.endswith(".scnt"):
                name = name[:-5]
            if name.endswith("."):
                name = name[:-1]
            scnt_files[name] = full_path
        logger.debug("Input Segment Copy Number Tensors (SCNT) identified as {input_scnts}".format(input_scnts=" , ".join(scnt_files.values())))
        scnts_by_name = {}
        segments_by_name = {}
        clone_ids_by_scnt = {}
        logger.info("Reading input SCNTs")
        for name, path in scnt_files.items():
            logger.debug("Reading SCNT from {file}".format(file=scnt_files[name]))
            segments, scnt = read_scnt_from_file(file_name=scnt_files[name], separator=args.separator)
            clone_ids_by_scnt[name] = sorted(scnt.keys())
            scnts_by_name[name] = scnt
            segments_by_name[name] = segments
        if len(scnts_by_name.values()) == 1:
            logger.warning("Only one input SCNT identified. Doing nothing with it, outputting as is.")
            aligned_segments_by_name, aligned_scnts_by_name = segments_by_name, scnts_by_name
        else:
            logger.info("Aligning input SCNTs.")
            aligned_segments_by_name, aligned_scnts_by_name = aligned_scnts(segments_by_sample_names=segments_by_name, scnts_by_sample_names=scnts_by_name)
        result_base_names = {}
        cnt = 0
        for name in sorted(scnt_files.keys()):
            new_name = name
            if name in result_base_names:
                new_name = name + str(cnt)
                cnt += 1
            new_name = new_name + "." + args.output_suffix
            result_base_names[name] = new_name
        output_dir = args.output_dir if args.output_dir != "" else os.getcwd()
        output_dir = get_full_path(path=output_dir)
        logger.info("Writing aligned SCNTs")
        for name, new_name in result_base_names.items():
            scnt = aligned_scnts_by_name[name]
            segments = aligned_segments_by_name[name]
            scnt_path = os.path.join(output_dir, new_name + "rck.scnt.tsv")
            logger.debug("Writing aligned SCNT {scnt_name} to {file}".format(scnt_name=name, file=scnt_path))
            write_scnt_to_file(file_name=scnt_path, segments=segments, scnt=scnt, separator=args.separator)
    elif args.command == "filter":
        logger.info("Filtering input segments from following sources {sources}".format(sources=args.scnt))
        segments = stream_segments_from_source(source=args.scnt, separator=args.separator, extra_separator=args.extra_separator)
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
        segments = filter_segments_by_chromosomal_regions(segments=segments, include=include_regions, exclude=exclude_regions,
                                                          include_full=args.include_full, exclude_full=args.exclude_full)
        keep_extra_field_entries = args.keep_extra_field_regex if args.keep_extra_field_regex is not None else []
        if args.keep_extra_field_regex_file is not None:
            keep_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.keep_extra_field_regex_file)))
        remove_extra_field_entries = args.remove_extra_field_regex if args.remove_extra_field_regex is not None else []
        if args.remove_extra_field_regex_file is not None:
            remove_extra_field_entries.extend(list(iter_over_string_entries_from_source(source=args.remove_extra_field_regex_file)))
        keep_extra_field = get_extra_field_regexes(string_entries=keep_extra_field_entries)
        remove_extra_field = get_extra_field_regexes(string_entries=remove_extra_field_entries)
        segments = filter_segments_by_extra(segments=segments, keep_extra_field=keep_extra_field, keep_extra_field_missing_strategy=args.keep_extra_field_missing_strategy,
                                            remove_extra_field=remove_extra_field, remove_extra_field_missing_strategy=args.remove_extra_field_missing_strategy)
        segments = filter_segments_by_size(segments=segments, min_size=args.min_size, max_size=args.max_size)
        write_segments_to_destination(destination=args.output, segments=segments)

    elif args.command == "haploid":
        segments = stream_segments_from_source(source=args.scnt, separator=args.separator, extra_separator=args.extra_separator)
        haploid_segments = iter_haploid_segments(segments=segments, copy=False)
        write_segments_to_destination(destination=args.output, segments=haploid_segments)
    elif args.command == "distance":
        clone_ids = args.clone_ids
        if args.clone_ids is not None:
            clone_ids = args.clone_ids.split(",")
        segments1, scnt1 = read_scnt_from_source(source=args.scnt1, clone_ids=clone_ids, separator=args.scnt1_separator,
                                                 extra_separator=args.scnt1_extra_separator, remove_cn_data_from_segs=True)
        segments2, scnt2 = read_scnt_from_source(source=args.scnt2, clone_ids=clone_ids, separator=args.scnt2_separator,
                                                 extra_separator=args.scnt2_extra_separator, remove_cn_data_from_segs=True)
        segments_by_sample_names = {"1": segments1, "2": segments2}
        scnts_by_sample_names = {"1": scnt1, "2": scnt2}
        segments_by_sample_names, scnts_by_sample_names = aligned_scnts(segments_by_sample_names=segments_by_sample_names,
                                                                        scnts_by_sample_names=scnts_by_sample_names)
        segments = segments_by_sample_names["1"]
        scnt1, scnt2 = scnts_by_sample_names["1"], scnts_by_sample_names["2"]
        distance = cn_distance_inter_scnt(tensor1=scnt1, tensor2=scnt2, segments=segments, check_clone_ids_match=True)
        print("distance = ", distance)

    logger.info("Success!")


if __name__ == "__main__":
    main()
