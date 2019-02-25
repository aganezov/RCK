import argparse
import datetime
import logging
import os
import re
import shutil
import sys
from copy import deepcopy

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import rck
from rck.utils.adj.process import refined_adjacencies_reciprocal
from rck.core.io import read_adjacencies_from_file, write_acnt_to_file, get_logging_cli_parser, get_standard_logger_from_args, read_adjacency_groups_from_file, \
    read_segments_from_file, \
    extract_scnt_from_segments, extract_scnb_from_segments, read_scnb_from_file, get_full_path, write_scnt_to_file, EXTERNAL_NA_ID, FALSE_POSITIVE, remove_cn_data_from_segments, \
    write_scnb_to_file, write_adjacencies_to_file, write_segments_to_file, remove_cn_data_from_adjacencies, read_positions_from_file, \
    write_positions_to_file, remove_cnb_data_from_segments
from rck.core.process import positions_aligned, adj_groups_concur
from rck.core.structures import get_segments_for_fragments_ids_dict, get_ref_telomeres_from_segments, get_ref_adjacencies_from_segments, SegmentCopyNumberBoundaries, refined_scnt, \
    refined_scnb, refined_scnt_with_adjacencies_and_telomeres, extract_spanned_extremities, SCNBoundariesStrategies, LengthSpreadRelationships, AdjacencyType, AdjacencyGroupType, \
    Phasing, Haplotype, CNBoundaries, AdjacencyCopyNumberProfile, SegmentCopyNumberProfile
from rck.core.graph import construct_hiag_inflate_from_haploid_data, IntervalAdjacencyGraph


def main():
    parser = argparse.ArgumentParser(parents=[get_logging_cli_parser()])
    parser.add_argument('--version', action='version', version=rck.version)

    parser.add_argument("--scnt", required=True)
    parser.add_argument("--scnt-separator", default="\t")
    parser.add_argument("--scnt-extra-separator", default=";")

    parser.add_argument("--adjacencies", required=True)
    parser.add_argument("--adjacencies-separator", default="\t")
    parser.add_argument("--adjacencies-extra-separator", default=";")

    parser.add_argument("--workdir", default=None)

    parser.add_argument("--clone-ids", default=None)

    parser.add_argument("--scnb", default=None)
    parser.add_argument("--scnb-separator", default="\t")
    parser.add_argument("--scnb-extra-separator", default=";")

    parser.add_argument("--telomere-segments", default=None)
    parser.add_argument("--telomere-segments-separator", default="\t")
    parser.add_argument("--telomere-segments-extra-separator", default=";")

    parser.add_argument("--telomere-positions", default=None)
    parser.add_argument("--telomere-positions-separator", default="\t")
    parser.add_argument("--telomere-positions-no-reciprocal", action="store_false", dest="telomere_reciprocal")

    parser.add_argument("--adjacency-groups", default=None)
    parser.add_argument("--adjacency-group-separator", default="\t")
    parser.add_argument("--adjacency-group-extra-separator", default=";")

    parser.add_argument("--fragments", default=None)
    parser.add_argument("--fragments-separator", default="\t")
    parser.add_argument("--fragments-extra-separator", default=";")

    # parser.add_argument("--trees", default=None)
    parser.add_argument("--no-allow-unit-segments", action="store_false", dest="allow_unit_segments")
    ###
    pre_group = parser.add_argument_group()
    pre_group.add_argument("--no-pre", action="store_false", dest="do_preprocess")
    ##
    # Arguments for Segment Copy Number Tensor preprocessing
    ##
    pre_group.add_argument("--pre-no-scnt", action="store_false", dest="do_pre_scnt")
    pre_group.add_argument("--pre-scnt-no-merge-fragments", action="store_false", dest="pre_scnt_merge_fragments")
    pre_group.add_argument("--pre-scnt-max-merge-gap", type=int, default=1000000000)
    pre_group.add_argument("--pre-scnt-no-fill-gaps", action="store_false", dest="pre_scnt_fill_gaps")
    pre_group.add_argument("--pre-scnt-max-fill-gap", type=int, default=1000000000)
    ##
    # Arguments for Segment Copy Number Boundaries
    ##
    pre_group.add_argument("--pre-scnb", action="store_true", dest="do_pre_bnd")
    pre_group.add_argument("--pre-scnb-strategy", choices=[strategy.value for strategy in SCNBoundariesStrategies], type=SCNBoundariesStrategies.from_string,
                           default=SCNBoundariesStrategies.UNIFORM_MIN_MAX.value)
    pre_group.add_argument("--pre-scnb-uniform-spread-size", type=int, default=1)
    pre_group.add_argument("--pre-scnb-length-spread-relation", choices=[rel.value for rel in LengthSpreadRelationships], type=LengthSpreadRelationships.from_string,
                           default=LengthSpreadRelationships.DUMMY.value)
    pre_group.add_argument("--pre-scnb-uniform-min", type=int, default=0)
    pre_group.add_argument("--pre-scnb-uniform-max", type=int, default=10)
    pre_group.add_argument("--pre-scnb-min-allow-zero-for-positive", type=int, default=-1)
    pre_group.add_argument("--pre-scnb-max-allow-zero-for-positive", type=int, default=1000000000)
    pre_group.add_argument("--pre-scnb-min-allow-positive-for-zero", type=int, default=-1)
    pre_group.add_argument("--pre-scnb-max-allow-positive-for-zero", type=int, default=1000000000)
    pre_group.add_argument("--pre-scnb-is-male", action="store_false", dest="pre_scnb_is_female")
    ##
    # Arguments for adjacencies
    ##
    pre_group.add_argument("--pre-no-adj", action="store_false", dest="do_pre_adj")
    pre_group.add_argument("--pre-adj-no-reciprocal", action="store_false", dest="pre_adj_reciprocal")
    pre_group.add_argument("--pre-adj-reciprocal-include-ref", action="store_true", dest="pre_adj_reciprocal_include_ref")
    pre_group.add_argument("--pre-adj-reciprocal-max-distance", type=int, default=50)
    ###
    run_group = parser.add_argument_group()
    run_group.add_argument("--no-run", action="store_false", dest="do_run")
    run_group.add_argument("--run-g-mip-gap", type=float, default=0.015)
    run_group.add_argument("--run-g-time-limit", type=int, default=28800)
    run_group.add_argument("--run-g-threads", type=int, default=4)
    run_group.add_argument("--run-g-mip-focus", type=int, choices=[0, 1, 2, 3], default=0)
    run_group.add_argument("--run-nas-fp", type=float, default=0.1)
    run_group.add_argument("--run-group-m-default-fp", type=float, default=0.1)
    run_group.add_argument("--run-group-n-default-fp", type=float, default=0.1)
    run_group.add_argument("--run-segment-length-attr", choices=["length", "length_10", "length_100", "length_1000"], default="length_10")
    run_group.add_argument("--run-g-allow-interrupted", action="store_true")
    ###
    output_group = parser.add_argument_group()
    output_group.add_argument("--o-prefix-name", type=str, dest="out_prefix_name", default="")
    output_group.add_argument("--o-scnt-separator", default="\t")
    output_group.add_argument("--o-scnt-extra-separator", default=";")
    output_group.add_argument("--o-acnt-separator", default="\t")
    output_group.add_argument("--o-acnt-extra-separator", default=";")
    output_group.add_argument("--o-acnt-mix-novel-and-reference", action="store_true",)
    ###
    post_group = parser.add_argument_group()
    post_group.add_argument("--post-check-all", action="store_true", dest="post_check_all")
    post_group.add_argument("--post-check-scnb", action="store_true", dest="post_check_scnb")
    post_group.add_argument("--post-check-labeling", action="store_true", dest="post_check_labeling")
    post_group.add_argument("--post-check-balancing", action="store_true", dest="post_check_balancing")
    post_group.add_argument("--post-check-adj-groups", action="store_true", dest="post_check_adj_groups")
    post_group.add_argument("--post-check-adj-groups-m", action="store_true", dest="post_check_adj_groups_m")
    post_group.add_argument("--post-check-adj-groups-n", action="store_true", dest="post_check_adj_groups_n")
    post_group.add_argument("--post-check-nas-fp", action="store_true", dest="post_check_nas_fp")
    ###
    args = parser.parse_args()

    from rck.core.ilp_gurobi import OptModelMultiClone, DEFAULT_GROUP_M_FP, SEGMENT_LENGTH_ATTRIBUTE
    import gurobi as g
    logger = get_standard_logger_from_args(args=args, program_name="RCK")

    workdir = args.workdir if args.workdir is not None else "RCK-{date}".format(date=str(datetime.date.today()))
    workdir_path = get_full_path(path=workdir)
    logger.debug("Working directory is set as {workdir_path}".format(workdir_path=workdir_path))

    logger.debug("Creating working directory if does not exist")
    os.makedirs(workdir_path, exist_ok=True)
    internal_debug_log_file = os.path.join(workdir_path, "debug.log")
    formatter = logging.Formatter(args.log_format)
    fh = logging.FileHandler(internal_debug_log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    raw_input_dir_path = get_full_path(path=os.path.join(workdir_path, "raw_input"))
    logger.debug("Raw input directory is {raw_dir_path}".format(raw_dir_path=raw_input_dir_path))

    logger.debug("Creating raw input directory if does not exist")
    os.makedirs(raw_input_dir_path, exist_ok=True)

    scnt_file_path = get_full_path(path=args.scnt)
    logger.debug("Allele-specific segment copy number tensor file path {scnt_path}".format(scnt_path=scnt_file_path))
    logger.debug("Copying SCNT file {scnt_file} to raw input direcotry {raw_input_dir}".format(scnt_file=scnt_file_path, raw_input_dir=raw_input_dir_path))
    shutil.copy2(src=scnt_file_path, dst=raw_input_dir_path)
    logger.info("Reading allele-specific segment copy number tensor from {scnt_path}".format(scnt_path=scnt_file_path))
    segments = read_segments_from_file(file_name=scnt_file_path, separator=args.scnt_separator, extra_separator=args.scnt_extra_separator)
    try:
        logger.info("Trying to extract allele-specific SCNT from segment copy number tensor from {scnt_path}".format(scnt_path=scnt_file_path))
        scnt = extract_scnt_from_segments(segments=segments)
    except ValueError as e:
        logger.error("Unable to get the allele-specific segment copy number information from {scnt_path}. Terminating".format(scnt_path=scnt_file_path))
        raise e

    clone_ids = sorted(set(args.clone_ids.split(","))) if args.clone_ids is not None else sorted(scnt.keys())
    logger.info("Clones are set as {clone_ids}".format(clone_ids=",".join(clone_ids)))
    logger.debug("Clone ids in allele-specific segment copy number tensor are {clone_ids}.".format(clone_ids=sorted(scnt.keys())))
    if set(clone_ids) < set(scnt.keys()):
        logger.info("Segment copy number tensor has information about clones {scnt_clones}, while only clones {clone_ids} were specified for processing."
                    "".format(scnt_clones=sorted(scnt.keys()), clone_ids=clone_ids))
    elif set(clone_ids) > set(scnt.keys()):
        logger.error("Segment copy number tensor has information about clones {scnt_clones}, but only clones {clone_ids} were specified for processing. Error."
                     "".format(scnt_clones=sorted(scnt.keys()), clone_ids=clone_ids))
    new_scnt = {}
    for clone_id in clone_ids:
        if clone_id not in scnt:
            raise ValueError("Clone {clone_id} is not present for allele-specific segment copy number values".format(clone_id=clone_id))
        new_scnt[clone_id] = scnt[clone_id]
    scnt = new_scnt

    try:
        logger.info("Trying to extract segment copy number boundaries from {scnt_path}".format(scnt_path=scnt_file_path))
        scnb = extract_scnb_from_segments(segments=segments, clone_ids=clone_ids, allow_missing=True)
        logger.info("Successfully extracted (at least some) segment copy number boundaries from {scnt_path}".format(scnt_path=scnt_file_path))
    except ValueError:
        scnb = {clone_id: SegmentCopyNumberBoundaries() for clone_id in clone_ids}
        logger.info("Unsuccessfully tried to extract segment copy number boundaries from {scnt_path}. Not to worry just yet.".format(scnt_path=scnt_file_path))

    scnb_file_path = None if args.scnb is None else get_full_path(args.scnb)
    try:
        if scnb_file_path is not None:
            logger.debug("Copying Segment Copy Number Boundaries file {cnb_file} to raw input directory {raw_input_dir}".format(cnb_file=scnb_file_path,
                                                                                                                                raw_input_dir=raw_input_dir_path))
            shutil.copy2(src=scnb_file_path, dst=raw_input_dir_path)
            logger.info("Trying to read segment copy number boundaries data from file {scnb_path}.".format(scnb_path=scnb_file_path))
            segments_from_scnb_file, scnb_from_file = read_scnb_from_file(file_name=scnt_file_path, clone_ids=clone_ids, allow_missing=True,
                                                                          separator=args.scnb_separator, extra_separator=args.scnb_extra_separator)
            logger.info("Successfully extracted (at least some) segment copy numbers from {scnb_path}. Updating previous segment copy number boundaries data (if exists)."
                        "".format(scnb_path=scnb_from_file))
            for clone_id in clone_ids:
                scnb[clone_id].update(scnb_from_file[clone_id])
    except ValueError:
        logger.info("Unsuccessfully tried to extract segment copy number boundaries from {scnb_path}. Not to worry just yet.".format(scnb_path=scnb_file_path))

    if args.do_preprocess and args.do_pre_scnt:
        new_segments, scnt, segments_ids_mapping = refined_scnt(segments=segments, scnt=scnt,
                                                                merge_fragments=args.pre_scnt_merge_fragments, max_merge_gap=args.pre_scnt_max_merge_gap,
                                                                fill_gaps=args.pre_scnt_fill_gaps, max_fill_gap=args.pre_scnt_max_fill_gap)
        old_segment, segments = segments, new_segments
        scnb = refined_scnb(scnb=scnb, new_segments=segments, segments_ids_mapping=segments_ids_mapping)

    adjacencies_file_path = get_full_path(path=args.adjacencies)
    logger.debug("Adjacencies file path {file_path}".format(file_path=adjacencies_file_path))

    logger.debug("Copying NAs file {nas_file} to raw input directory {raw_input_dir}".format(nas_file=adjacencies_file_path, raw_input_dir=raw_input_dir_path))
    shutil.copy2(src=adjacencies_file_path, dst=raw_input_dir_path)

    logger.info("Reading Adjacencies from {file_path}".format(file_path=adjacencies_file_path))
    input_adjacencies = read_adjacencies_from_file(file_name=adjacencies_file_path, separator=args.adjacencies_separator, extra_separator=args.adjacencies_extra_separator)
    logger.debug("A total of {cnt} unlabeled novel adjacencies were obtained".format(cnt=len([a for a in input_adjacencies if a.adjacency_type == AdjacencyType.NOVEL])))
    logger.debug("A total of {cnt} unlabeled reference adjacencies were obtained".format(cnt=len([a for a in input_adjacencies if a.adjacency_type == AdjacencyType.REFERENCE])))
    if args.do_preprocess and args.do_pre_adj:
        logger.info("Preprocessing input adjacencies.")
        if args.pre_adj_reciprocal:
            input_ref = [adj for adj in input_adjacencies if adj.adjacency_type == AdjacencyType.REFERENCE]
            input_nov = [adj for adj in input_adjacencies if adj.adjacency_type == AdjacencyType.NOVEL]
            adjacencies = input_nov
            if args.pre_adj_reciprocal_include_ref:
                logger.debug("Both novel and reference input adjacencies are preprocessed w.r.t. reciprocality")
                adjacencies = adjacencies + input_ref
            logger.info("Preprocessing input adjacencies w.r.t. reciprocality.")
            adjacencies = refined_adjacencies_reciprocal(novel_adjacencies=adjacencies, max_distance=args.pre_adj_reciprocal_max_distance)
            if args.pre_adj_reciprocal_include_ref:
                input_adjacencies = adjacencies
            else:
                input_adjacencies = adjacencies + input_ref

    telomere_segments_file_path = None if args.telomere_segments is None else get_full_path(path=args.telomere_segments)
    logger.debug("Additional telomeres file path {tel_path}".format(tel_path=telomere_segments_file_path))
    if telomere_segments_file_path is not None:
        logger.debug("Copying telomere segments file {telomeres_file} to raw input directory {raw_input_dir}".format(telomeres_file=telomere_segments_file_path,
                                                                                                                     raw_input_dir=raw_input_dir_path))
        shutil.copy2(src=telomere_segments_file_path, dst=raw_input_dir_path)
        logger.info("Reading additional telomere locations (via segments) from {tel_path}".format(tel_path=telomere_segments_file_path))
        telomeres_segments = read_segments_from_file(file_name=telomere_segments_file_path, separator=args.telomere_segments_separator,
                                                     extra_separator=args.telomere_segments_extra_separator)
    else:
        telomeres_segments = []
        logger.info("No additional telomeres positions (via segments) were provided, only reference locations are allowed.")

    telomere_positions_file_path = None if args.telomere_positions is None else get_full_path(path=args.telomere_positions)
    if telomere_positions_file_path is not None:
        logger.debug("Copying telomere positions file {telomeres_file} to raw input directory {raw_input_dir}".format(telomeres_file=telomere_positions_file_path,
                                                                                                                      raw_input_dir=raw_input_dir_path))
        shutil.copy2(src=telomere_positions_file_path, dst=raw_input_dir_path)
        logger.info("Reading additional telomere locations (via positions) from {tel_path}".format(tel_path=telomere_positions_file_path))
        input_telomere_positions = read_positions_from_file(file_name=telomere_positions_file_path, separator=args.telomere_positions_separator)
    else:
        logger.info("No additional telomeres positions (via positions) were provided, only reference locations are allowed.")
        input_telomere_positions = []

    adjacency_group_file_path = None if args.adjacency_groups in ("", None) else get_full_path(path=args.adjacency_groups)
    logger.debug("Adjacency groups file path {nas_groups_file}".format(nas_groups_file=adjacency_group_file_path))
    if adjacency_group_file_path is not None:
        logger.info("Reading adjacencies groups from {nas_groups_file}".format(nas_groups_file=adjacency_group_file_path))
        adjacency_groups = read_adjacency_groups_from_file(file_name=adjacency_group_file_path)
    else:
        adjacency_groups = []
        logger.info("No information about adjacency groups was provided")

    fragments_file_path = None if args.fragments in (None, "") else get_full_path(path=args.fragments)
    if fragments_file_path is not None:
        logger.debug("Copying fragments file {fragments_file} to raw input directory {raw_input_dir}".format(fragments_file=fragments_file_path,
                                                                                                             raw_input_dir=raw_input_dir_path))
        shutil.copy2(src=fragments_file_path, dst=raw_input_dir_path)
        fragments = read_segments_from_file(file_name=fragments_file_path)
    else:
        fragments = None
        logger.debug("Fragments file path {fragments_file_path}".format(fragments_file_path=fragments_file_path))

    overall_nas_fp = args.run_nas_fp
    logger.info("Global False Positive value for novel adjacencies is set at {fp}".format(fp=overall_nas_fp))

    preprocessed_input_dir_path = get_full_path(os.path.join(workdir_path, "input"))
    os.makedirs(preprocessed_input_dir_path, exist_ok=True)
    if fragments is None:
        fragments = deepcopy(segments)

    if args.do_preprocess:
        logger.info("Refining input adjacencies and segments/fragments, so that extremities are aligned.")
        fragments = deepcopy(segments)
        segments, scnt, segments_ids_mapping = refined_scnt_with_adjacencies_and_telomeres(segments=segments, scnt=scnt, adjacencies=input_adjacencies,
                                                                                           telomere_positions=input_telomere_positions)
        scnb = refined_scnb(scnb=scnb, new_segments=segments, segments_ids_mapping=segments_ids_mapping)

    logger.info("Creating missing segment copy number boundaries will be automatically generated.")
    for clone_id in clone_ids:
        missing_only = not args.do_pre_bnd
        scnb[clone_id].fill(segments=segments, scnp=scnt[clone_id], missing_only=missing_only, strategy=args.pre_scnb_strategy,
                            min_allow_zero_for_positive=args.pre_scnb_min_allow_zero_for_positive,
                            max_allow_zero_for_positive=args.pre_scnb_max_allow_zero_for_positive,
                            min_allow_positive_for_zero=args.pre_scnb_min_allow_positive_for_zero,
                            max_allow_positive_for_zero=args.pre_scnb_max_allow_positive_for_zero,
                            uniform_spread_size=args.pre_scnb_uniform_spread_size,
                            length_spread_relation=args.pre_scnb_length_spread_relation,
                            uniform_min=args.pre_scnb_uniform_min,
                            uniform_max=args.pre_scnb_uniform_max,
                            is_female=args.pre_scnb_is_female)

    if telomeres_segments is not None:
        logger.info("Extracting additional permitted telomere locations")
        additional_telomeres = extract_spanned_extremities(source=segments, boundaries=telomeres_segments)
        logger.debug("Extracted {cnt} additional telomere locations".format(cnt=len(additional_telomeres)))
        input_telomere_positions.extend(additional_telomeres)

    logger.info("A total of input telomere locations is {cnt}".format(cnt=len(input_telomere_positions)))

    logger.info("Finished preprocessing. At this point every data part must be aligned with each other.")
    logger.info("Checking data alignment/concordance")
    if fragments is None:
        logger.error("FRAGMENTS must exist; their absence at this point indicates that preprocessing was turned off, while --fragments were not part of the original input")
        exit(1)
    if scnb is None:
        logger.error("Lower/Upper boundaries for Segment Copy Numbers must exist; "
                     "their absence at this point indicated that preprocessing was turned off, while --scn-boundaries were not part of the original input")
        exit(1)
    segments_to_fragments = None
    try:
        segments_to_fragments = get_segments_for_fragments_ids_dict(segments=segments, fragments=fragments, allow_non_covered=False)
    except Exception:
        logger.error("FRAGMENTS boundaries do not align with SEGMENTS boundaries; "
                     "having this problem at this point indicates that the preprocessing was turned off and the data in original input SCNT and --fragments is not aligned")
        exit(1)
    segments_positions = []
    for s in segments:
        segments_positions.append(s.start_position)
        segments_positions.append(s.end_position)
    nas_positions = []
    for na in input_adjacencies:
        nas_positions.append(na.position1)
        nas_positions.append(na.position2)
    logger.debug("Checking novel adjacency positions concordance")
    if not positions_aligned(segments_positions=segments_positions, other_positions=nas_positions):
        logger.error("NOVEL ADJACENCIES positions are not aligned with extremities of SEGMENTS; "
                     "having this problem at this point indicates that the preprocessing was turned off and the data in original input SCNT and NAS is not aligned")
        exit(1)
    telomeres = sorted(get_ref_telomeres_from_segments(segments=segments) + input_telomere_positions, key=lambda p: p.coordinate)
    logger.debug("Checking telomere positions concordance")
    if not positions_aligned(segments_positions=segments_positions, other_positions=telomeres):
        logger.error("TELOMERES positions are not aligned with extremities of SEGMENTS; "
                     "having this problem at this problem indicates that preprocessing was turned off and the data in the original SCNT and --telomeres is not aligned")
        exit(1)
    logger.debug("Checking adjacency groups concordance")
    if len(adjacency_groups) > 0 and not adj_groups_concur(adj_groups=adjacency_groups, adjacencies=input_adjacencies):
        logger.error("Some novel adjacencies group in NAS-GROUPS refers to an absent ADJACENCY; "
                     "having this problem at this time indicates that either the original --nas-groups has a group that is referring to the missing novel adjacency in NAS, "
                     "or the referred novel adjacency from NAS was removed during the preprocessing.")
        exit(1)
    logger.debug("Populating adjacency groups")
    adjacencies_by_external_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in input_adjacencies}
    for ag in adjacency_groups:
        ag.populate_adjacencies_via_ids(source=input_adjacencies, source_by_ids=adjacencies_by_external_ids)

    ########
    #
    # Creating files with data, that can be used without running preprocessing in the future
    #
    ########

    preprocessed_scnt_file = os.path.join(preprocessed_input_dir_path, "rck.scnt.tsv")
    remove_cn_data_from_segments(segments=segments)
    remove_cnb_data_from_segments(segments=segments)
    logger.info("Writing (preprocessed) allele-specific segment copy number data to {file}".format(file=preprocessed_scnt_file))
    write_scnt_to_file(file_name=preprocessed_scnt_file, segments=segments, scnt=scnt, clone_ids=clone_ids, inplace=False)

    preprocessed_scnb_file = os.path.join(preprocessed_input_dir_path, "rck.scnb.tsv")
    logger.info("Writing (preprocessed) segment copy number boundaries information to {file}".format(file=preprocessed_scnb_file))
    write_scnb_to_file(file_name=preprocessed_scnb_file, segments=segments, scnb=scnb, clone_ids=clone_ids, inplace=False)

    preprocessed_input_adjacencies_file = os.path.join(preprocessed_input_dir_path, "rck.adj.tsv")
    logger.info("Writing (preprocessed) input adjacencies data to {file}".format(file=preprocessed_input_adjacencies_file))
    write_adjacencies_to_file(file_name=preprocessed_input_adjacencies_file, adjacencies=input_adjacencies)
    remove_cn_data_from_adjacencies(input_adjacencies)

    preprocessed_fragments_file = os.path.join(preprocessed_input_dir_path, "rck.frag.tsv")
    remove_cn_data_from_segments(segments=fragments)
    logger.info("Writing (preprocessed) fragments data to {file}".format(file=preprocessed_fragments_file))
    write_segments_to_file(file_name=preprocessed_fragments_file, segments=fragments, extra=None)

    ########

    logger.debug("Extracting reference adjacencies from input segments")
    ref_adjacencies = get_ref_adjacencies_from_segments(segments=segments, assign_external_ids=True)
    logger.debug("A total of {cnt} reference unlabeled adjacencies are considered (labeled will be x2, except for X and Y chromosomes)".format(cnt=len(ref_adjacencies)))
    adjacencies = input_adjacencies + ref_adjacencies
    output_dir = os.path.join(workdir_path, "output")
    logger.debug("Output directory is {output}".format(output=output_dir))
    os.makedirs(output_dir, exist_ok=True)

    if args.telomere_reciprocal:
        logger.info("Reciprocal telomere positions are specified via --telomere-positions-reciprocal. Extracting.")
        iag = IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)
        iag.build_graph()
        for u, v in iag.ref_adjacency_edges(data=False):
            if u in telomeres and v not in telomeres:
                telomeres.append(v)
            if v in telomeres and u not in telomeres:
                telomeres.append(u)
    preprocessed_telomeres_file = os.path.join(preprocessed_input_dir_path, "rck.pos.tsv")
    logger.info("Writing (preprocessed) telomeres (both input and reference) to {file}".format(file=preprocessed_telomeres_file))
    write_positions_to_file(file_name=preprocessed_telomeres_file, positions=telomeres)

    logger.debug("Output directory is {output}".format(output=output_dir))
    os.makedirs(output_dir, exist_ok=True)

    ##########
    #
    # Gurobi model
    #
    ##########
    if not args.do_run:
        logger.info("A --no-run flag was set. Not performing the inference.")
        exit(0)

    gurobi_log_path = os.path.join(output_dir, "gurobi.log")
    logger.info("Gurobi log will be stored in {log_path}".format(log_path=gurobi_log_path))
    extra = {DEFAULT_GROUP_M_FP: args.run_group_m_default_fp,
             SEGMENT_LENGTH_ATTRIBUTE: args.run_segment_length_attr}
    logger.debug("Gurobi ILP extra is {extra}".format(extra=str(extra)))
    logger.info("Setting up Gurobi ILP model (includes construction of the IAG)")
    ilp_model = OptModelMultiClone(hapl_segments=segments,
                                   hapl_adjacencies=adjacencies,
                                   scnt=scnt,
                                   hapl_telomeres=telomeres,
                                   hapl_segments_to_fragments=segments_to_fragments,
                                   hapl_adjacencies_groups=adjacency_groups,
                                   scnb=scnb,
                                   hapl_nov_adjacencies_fp=overall_nas_fp,
                                   extra=extra)
    logger.debug("Building variables and constraints")
    ilp_model.build_gurobi_model()
    logger.debug("Setting gurobi parameters")
    ilp_model.gm.setParam("MIPGap", args.run_g_mip_gap)
    ilp_model.gm.setParam("MIPGapAbs", args.run_g_mip_gap)
    ilp_model.gm.setParam("MIPFocus", args.run_g_mip_focus)
    ilp_model.gm.setParam("LogFile", gurobi_log_path)
    ilp_model.gm.setParam("TimeLimit", args.run_g_time_limit)
    ilp_model.gm.setParam("Threads", args.run_g_threads)
    logger.info("Starting Gurobi to solve the optimization problem")
    ilp_model.solve_model()

    logger.info("Gurobi model solving has ended")
    status = ilp_model.gm.status
    solution_cnt = ilp_model.gm.solcount
    if status == g.GRB.Status.INFEASIBLE:
        logger.error("Constructed model was infeasible to solve. This usually happens because of the tight segment copy number boundaries.")
        logger.error("\t Try more flexible segment copy number boundaries or try increasing allowed Novel adjacencies False Positive value.")
        logger.info("Gurobi computes the IIS")
        ilp_model.gm.computeIIS()
        ilp_path = os.path.join(output_dir, "model.ilp")
        logger.info("Writing down Gurobi IIS to {file}".format(file=ilp_path))
        ilp_model.gm.write(ilp_path)
        logger.error("Inference was unsuccessful.")
        exit(1)
    allowed_statuses = [g.GRB.Status.OPTIMAL, g.GRB.Status.TIME_LIMIT]
    if args.run_g_allow_interrupted:
        allowed_statuses.append(g.GRB.Status.INTERRUPTED)
    if status not in allowed_statuses:
        logger.error("Gurobi finished with status {status}".format(status=status))
        logger.error("Inference was unsuccessful")
        exit(1)
    if solution_cnt == 0:
        logger.error("Inference was unsuccessful")
        exit(1)

    logger.info("Gurobi finished execution with status {status}".format(status=status))
    logger.info("Extracting inferred diploid segment copy number data")
    scnt = ilp_model.get_scnt_from_model()
    logger.info("Extracting inferred diploid adjacency copy number data")
    acnt = ilp_model.get_acnt_from_model()
    scnt_base_name = args.out_prefix_name + "rck.scnt.tsv"
    if scnt_base_name.startswith("."):
        scnt_base_name = scnt_base_name[1:]
    scnt_base_name = re.sub("\.+", ".", scnt_base_name)
    scnt_file_path = os.path.join(output_dir, scnt_base_name)
    logger.info("Writing inferred diploid segment copy number data to {file}".format(file=scnt_file_path))
    write_scnt_to_file(file_name=scnt_file_path, segments=segments, scnt=scnt)
    acnt_base_name = args.out_prefix_name + "rck.acnt.tsv"
    if acnt_base_name.startswith("."):
        acnt_base_name = acnt_base_name[1:]
    acnt_base_name = re.sub("\.+", ".", acnt_base_name)
    acnt_file_path = os.path.join(output_dir, acnt_base_name)
    logger.info("Writing inferred diploid adjacency copy number data to {file}".format(file=acnt_file_path))
    write_acnt_to_file(file_name=acnt_file_path, acnt=acnt, adjacencies=adjacencies, output_reference=True, mix_reference_and_novel=args.o_acnt_mix_novel_and_reference)

    #########
    #
    # POST PROCESSING
    #
    #########

    if args.post_check_all or args.post_check_scnb:
        logger.info("Performing post-inference check that inferred segment copy number values are within input and/or preprocessed bounds")
        for clone_id in clone_ids:
            logger.info("Working with clone {clone_id}".format(clone_id=clone_id))
            scnp = scnt[clone_id]
            scnbp = scnb[clone_id]
            problems = False
            for segment in segments:
                sid = segment.stable_id_non_hap
                sync_indicator = ilp_model.alleles_sync_result(segment=segment)
                cna, cnb = scnp.get_cn(sid=sid, haplotype=Haplotype.A), scnp.get_cn(sid=sid, haplotype=Haplotype.B)
                lower_a = scnbp.get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.LOWER)
                lower_b = scnbp.get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.LOWER)
                upper_a = scnbp.get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.UPPER)
                upper_b = scnbp.get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.UPPER)
                if sync_indicator == 1:
                    if (not lower_a <= cna <= upper_a) or (not lower_b <= cnb <= upper_b):
                        problems = True
                else:
                    if (not lower_b <= cna <= upper_b) or (not lower_a <= cnb <= upper_a):
                        problems = True
                if problems:
                    logger.error("Something went WRONG! For segment {sid} with allele-sync flag {flag} inferred copy numbers are A={cna}, B={cnb},"
                                 "the input copy number boundaries were A={lower_a}-{upper_a}, B={lower_b}-{upper_b}."
                                 "".format(sid=sid, flag=sync_indicator, cna=cna, cnb=cnb, lower_a=lower_a, upper_a=upper_a, lower_b=lower_b, upper_b=upper_b))
            if not problems:
                logger.info("Everything is OK!")

    if args.post_check_all or args.post_check_labeling:
        logger.info("Performing post-inference check on labels on adjacencies with positive copy numbers")
        logger.info("Checking that every unlabeled novel adjacency has at most one labeled 'realization' across all clones")
        novel_adjacencies = [a for a in adjacencies if a.adjacency_type == AdjacencyType.NOVEL]
        violations = False
        for adj in novel_adjacencies:
            labeling_realizations = set()
            for clone_id in clone_ids:
                acnp = acnt[clone_id]
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    if acnp.get_cn(aid=adj.stable_id_non_phased, phasing=ph) > 0:
                        labeling_realizations.add(ph)
            if len(labeling_realizations) > 1:
                violations = True
                logger.error("Something went WRONG! For unlabeled novel adjacency (id = {aid}) {stable_id} more than 1 labeled 'realization' is present across all clones"
                             "".format(aid=adj.extra.get(EXTERNAL_NA_ID, adj.idx), stable_id=adj.stable_id_non_phased))
        if not violations:
            logger.info("Everything is OK! No violations of multiple 'realizations' were identified.")

        logger.info("Checking that every reference location, that has reciprocal novel adjacencies, concurs with the Infinite Sites constraints")
        combined_scnp = SegmentCopyNumberProfile.combined(*[scnt[clone_id] for clone_id in clone_ids])
        combined_acnp = AdjacencyCopyNumberProfile.combined(*[acnt[clone_id] for clone_id in clone_ids])
        hiag = construct_hiag_inflate_from_haploid_data(hapl_segments=segments, hapl_adjacencies=adjacencies)
        hiag.assign_copy_numbers_from_scn_profile(scn_profile=combined_scnp)
        hiag.assign_copy_numbers_from_acn_profile(acn_profile=combined_acnp)
        hiag.remove_edges_with_zero_cn()

        logger.info("Checking for extremity-exclusivity constraint violations")
        extremity_exclusivity_violations = hiag.violations_of_extremity_exclusivity()
        if len(extremity_exclusivity_violations) > 0:
            logger.error("Something went WRONG! The following segments extremities violate extremity-exclusivity constraint of the Infinite Sites model.")
            logger.error(",".join(map(str, extremity_exclusivity_violations)))
        else:
            logger.info("Everything is OK! No extremity-exclusivity violations were found!")
        logger.info("Checking for homologous-extremity-exclusivity violations")
        homologous_extremity_exclusivity_violations = hiag.violations_of_homologous_extremity_exclusivity()
        if len(homologous_extremity_exclusivity_violations) > 0:
            logger.error("Something went WRONG! The following pairs of segments extremities violate homologous-extremity-exclusivity constraint of the Infinite Sites model.")
            string_data = []
            for p1, p2 in homologous_extremity_exclusivity_violations:
                string_data.append("{{{p1},{p2}}}".format(p1=str(p1), p2=str(p2[0])))
            logger.error(" , ".join(string_data))
        else:
            logger.info("Everything is OK! No homologous-extremity-exclusivity violations were found!"
                        "".format(cnt=len(ilp_model.reciprocal_locations)))
        logger.info("Checking for reciprocal-homologous-extremity-exclusivity violations")
        homologous_reciprocal_extremity_exclusivity_violations = hiag.violations_of_homologous_reciprocal_extremity_exclusivity()
        if len(homologous_reciprocal_extremity_exclusivity_violations) > 0:
            logger.error("Something went WRONG! "
                         "The following pairs of segments extremities violate homologous-reciprocal-extremity-exclusivity constraint of the Infinite Sites model.")
            string_data = []
            for p1, p2 in homologous_reciprocal_extremity_exclusivity_violations:
                string_data.append("{{{p1},{p2}}}".format(p1=str(p1), p2=str(p2)))
            logger.error(" , ".join(string_data))
        else:
            logger.info("Everything is OK! Ho homologous-reciprocal-extremity-exclusivity were found across {cnt} reciprocal locations!"
                        "".format(cnt=len(ilp_model.reciprocal_locations)))

    if args.post_check_all or args.post_check_balancing:
        logger.info("Performing post-inference check on balances/excesses on segments' extremities")
        for clone_id in clone_ids:
            logger.info("Processing clone {clone_id}".format(clone_id=clone_id))
            hiag = construct_hiag_inflate_from_haploid_data(hapl_segments=segments, hapl_adjacencies=adjacencies)
            scnp = scnt[clone_id]
            acnp = acnt[clone_id]
            hiag.assign_copy_numbers_from_scn_profile(scn_profile=scnp)
            hiag.assign_copy_numbers_from_acn_profile(acn_profile=acnp)
            hiag.remove_edges_with_zero_cn()
            logger.info("Checking that every vertex has a copy number excess >= 0.")
            for node in hiag.nodes(data=False):
                if hiag.node_imbalance(node=node) < 0:
                    logger.error("Something went WRONG! On segment extremity {node} there is a negative copy number excess...")
                    exit(1)
            logger.info("Getting inferred telomeres.")
            diploid_telomeres = hiag.get_telomeres()
            inferred_hapl_telomeres_ids = {p.stable_id_non_hap for p in diploid_telomeres}
            input_hapl_telomers_ids = {p.stable_id_non_hap for p in telomeres}
            if inferred_hapl_telomeres_ids > input_hapl_telomers_ids:
                logger.error("Something went WRONG! Following segments extremities, while not specified specified as possible telomere sites were inferred as such.")
                logger.error(",".join(map(str, sorted(input_hapl_telomers_ids - input_hapl_telomers_ids))))
            else:
                logger.info("Everything is OK! in clone {clone_id} all extremities have non-negative copy number excess, and inferred telomere sites concur with the input"
                            "".format(clone_id=clone_id))

    if args.post_check_all or args.post_check_adj_groups or args.post_check_adj_groups_m:
        logger.info("Performing post-inference check on adjacencies groups (molecule type)")
        molecule_groups = [ag for ag in adjacency_groups if ag.group_type == AdjacencyGroupType.MOLECULE]
        molecule_groups_fine_cnt = 0
        logger.info("There were {cnt} molecule adjacency groups in the input.".format(cnt=len(molecule_groups)))
        for group in molecule_groups:
            group_fp = group.extra.get(FALSE_POSITIVE, extra[DEFAULT_GROUP_M_FP])
            group_is_good = False
            good_inferred_fp = -1
            for clone_id in clone_ids:
                acnp = acnt[clone_id]
                clone_present = acnp.haploid_adjacencies_present(adjacencies=group.adjacencies)
                inferred_fp = 1 - (len(clone_present) * 1.0 / len(group.adjacencies))
                if inferred_fp <= group_fp:
                    good_inferred_fp = inferred_fp
                group_is_good |= inferred_fp <= group_fp
            if not group_is_good:
                logger.error("Something went WRONG! In Adjacency Group {gid} not a single clone has <= than the input FP of {input_fp}."
                             "".format(gid=group.gid, input_fp=group_fp))
            else:
                logger.debug("Group {gid} is fine (there is a clone with group FP of {fp:0.4f} which is >= than the input FP of {ifp:0.4f})."
                             "".format(gid=group.gid, fp=good_inferred_fp, ifp=group.extra.get(FALSE_POSITIVE, extra[DEFAULT_GROUP_M_FP])))
                molecule_groups_fine_cnt += 1
        if molecule_groups_fine_cnt == len(molecule_groups):
            logger.info("Everything is OK for all {good}/{all} molecule adjacency groups.".format(good=molecule_groups_fine_cnt, all=len(molecule_groups)))
        else:
            logger.error("Something went WRONG! In some molecule adjacency groups (see above) not a single clone has <= than the input FP of.")

    if args.post_check_all or args.post_check_nas_fp:
        logger.info("Performing post-inference check on overall novel adjacencies false positive parameter")
        novel_adjacencies = [a for a in adjacencies if a.adjacency_type == AdjacencyType.NOVEL]
        present = set()
        for clone_id in clone_ids:
            acnp = acnt[clone_id]
            clone_present = acnp.haploid_adjacencies_present(adjacencies=novel_adjacencies)
            clone_present_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.idx) for adj in clone_present}
            present |= clone_present_ids
        inferred_fp = 1 - (len(present) * 1.0 / len(novel_adjacencies))
        if inferred_fp <= overall_nas_fp:
            logger.info("Everything is OK! Inferred FP = {inf_fp:0.4f}, while the input FP was {input_fp:0.4f}.".format(inf_fp=inferred_fp, input_fp=overall_nas_fp))
        else:
            logger.error("Something went WRONG! Inferred FP = {inf_fp:0.4f}, while the input FP was {input_fp:0.4f}.".format(inf_fp=inferred_fp, input_fp=overall_nas_fp))
        logger.info("Finished everything. Hope you've enjoyed using RCK!")


if __name__ == "__main__":
    main()
