import argparse
import os
import re
import sys
import shutil
import gurobi as g



sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from dassp.core.io import read_scn_tensor, get_all_clone_ids_from_dasck_scnt_file, read_novel_adjacencies, read_positions, read_segments, read_scn_boundaries, write_scn_tensor, write_acnt
from dassp.core.preprocess import positions_aligned, nas_groups_aligned
from dassp.core.structures import get_segments_for_fragments_ids_dict, get_ref_telomeres_from_segments, get_ref_adjacencies_from_segments, construct_nas_groups
from dassp.algo.ilp import OptModelMultiClone


def get_full_path(path):
    return os.path.abspath(os.path.expanduser(path))


def create_dir_is_doesnt_exist(dirpath):
    if not os.path.exists(path=dirpath):
        os.makedirs(dirpath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scnt")
    parser.add_argument("nas")
    parser.add_argument("workdir")
    parser.add_argument("--telomeres", default="")
    parser.add_argument("--nas-groups", default="")
    parser.add_argument("--fragments", default="")
    parser.add_argument("--scn-boundaries", default="")
    parser.add_argument("--trees", default="")
    parser.add_argument("--no-allow-unit-segments", action="store_false", dest="allow_unit_segments")
    ###
    pre_group = parser.add_argument_group()
    pre_group.add_argument("--no-pre", action="store_false", dest="do_preprocess")
    ##
    # Arguments for Novel adjacencies preprocessing
    ##
    pre_group.add_argument("--pre-no-nas-refine", action="store_false", dest="do_pre_nas_refine")
    pre_group.add_argument("--pre-nas-no-remove-short", action="store_false", dest="pre_nas_remove_short")
    pre_group.add_argument("--pre-nas-short-min-size", type=int, default=500)
    pre_group.add_argument("--pre-nas-no-reciprocal-merging", action="store_false", dest="pre_nas_reciprocal_merging")
    pre_group.add_argument("--pre-nas-max-reciprocal-window", type=int, default=50)
    pre_group.add_argument("--pre-nas-no-check-is", action="store_false", dest="pre_nas_check_is")
    ##
    # Arguments for Segment Copy Number Tensor preprocessing
    ##
    pre_group.add_argument("--pre-no-scnt-refine", action="store_false", dest="do_pre_scnt_refine")
    pre_group.add_argument("--pre-scnt-no-merge-fragments", action="store_false", dest="pre_scnt_merge_fragments")
    pre_group.add_argument("--pre-scnt-max-merge-gap", type=int, default=1000000)
    pre_group.add_argument("--pre-scnt-no-fill-gaps", action="store_false", dest="pre_scnt_fill_gaps")
    pre_group.add_argument("--pre-scnt-max-fill-gap", type=int, default=1000000)
    ##
    #
    ##
    pre_group.add_argument("--pre-no-bnd", action="store_false", dest="do_pre_bnd")
    pre_group.add_argument("--pre-bnd-strategy", choices=["fixed", "uniform-spread", "length-spread", "uniform-min-max"], default="length-spread")
    pre_group.add_argument("--pre-bnd-uniform-spread-size", type=int, default=1)
    pre_group.add_argument("--pre-bnd-length-spread-relation", choices=["dummy"], default="dummy")
    pre_group.add_argument("--pre-bnd-uniform-min", type=int, default=0)
    pre_group.add_argument("--pre-bnd-uniform-max", type=int, default=10)
    pre_group.add_argument("--pre-bnd-min-allow-zero-for-positive", type=int, default=-1)
    pre_group.add_argument("--pre-bnd-max-allow-zero-for-positive", type=int, default=1000000000)
    pre_group.add_argument("--pre-bnd-min-allow-positive-for-zero", type=int, default=-1)
    pre_group.add_argument("--pre-bnd-max-allow-positive-for-zero", type=int, default=1000000000)
    pre_group.add_argument("--pre-bnd-is-male", action="store_false", dest="pre_bnd_is_female")
    ##
    #
    ##
    pre_group.add_argument("--pre-tel-no-hum-centromeres", action="store_false", dest="pre_tel_use_centromeres")
    ###
    run_group = parser.add_argument_group()
    run_group.add_argument("--no-run", action="store_false", dest="do_run")
    run_group.add_argument("--run-no-trees", action="store_false", dest="run_trees")
    run_group.add_argument("--run-g-mip-gap", type=float, default=0.015)
    run_group.add_argument("--run-g-time-limit", type=int, default=28800)
    run_group.add_argument("--run-nas-fps", type=str, default="0.0,0.25,0.5")
    ###
    output_group = parser.add_argument_group()
    output_group.add_argument("--o-prefix-name", type=str, dest="out_prefix_name", default="")
    output_group.add_argument("--o-tree-dir-template", type=str, dest="out_tree_dir_basename_template", default="tree_{tree_id}")
    output_group.add_argument("--o-scnt-basename-template", type=str, dest="out_scnt_basename_template", default="{prefix}.{fp:.2f}.scnt.tsv")
    output_group.add_argument("--o-acnt-basename-template", type=str, dest="out_acnt_basename_template", default="{prefix}.{fp:.2f}.acnt.tsv")
    args = parser.parse_args()
    scnt_file_path = get_full_path(path=args.scnt)
    clone_ids = get_all_clone_ids_from_dasck_scnt_file(file_name=scnt_file_path, allow_unit_segments=args.allow_unit_segments)
    segments, scnt = read_scn_tensor(file_name=scnt_file_path, clone_ids=clone_ids, allow_unit_segments=args.allow_unit_segments)
    nas_file_path = get_full_path(path=args.nas)
    nas = read_novel_adjacencies(file_name=nas_file_path)
    workdir_path = get_full_path(path=args.workdir)
    create_dir_is_doesnt_exist(dirpath=workdir_path)
    raw_input_dir_path = get_full_path(path=os.path.join(workdir_path, "raw_input"))
    create_dir_is_doesnt_exist(dirpath=raw_input_dir_path)
    shutil.copy2(src=scnt_file_path, dst=raw_input_dir_path)
    shutil.copy2(src=nas_file_path, dst=raw_input_dir_path)
    telomeres_file_path = None if args.telomeres == "" else get_full_path(path=args.telomeres)
    telomeres = None if telomeres_file_path is None else read_positions(file_name=telomeres_file_path)
    nas_groups_file_path = None if args.nas_groups == "" else get_full_path(path=args.nas_groups)
    # TODO: complete
    nas_groups_ids = None
    fragments_file_path = None if args.fragments == "" else get_full_path(path=args.fragments)
    fragments = None if fragments_file_path is None else read_segments(file_name=fragments_file_path)
    scn_boundaries_file_path = None if args.scn_boundaries == "" else get_full_path(path=args.scn_boundaries)
    scn_boundaries = None if scn_boundaries_file_path is None else read_scn_boundaries(file_name=scn_boundaries_file_path, segments=segments, clone_ids=clone_ids,
                                                                                       allow_unit_segments=args.allow_unit_segments)
    trees_file_path = None if args.trees == "" else get_full_path(path=args.trees)
    nas_fps = list(map(lambda entry: float(entry), args.run_nas_fps.split(",")))
    # TODO: complete
    trees = None
    for file_path in [telomeres_file_path, nas_groups_file_path, fragments_file_path, scn_boundaries_file_path, trees_file_path]:
        if file_path is not None:
            shutil.copy2(src=file_path, dst=raw_input_dir_path)
    tmp_dir_path = get_full_path(os.path.join(workdir_path, "tmp"))
    create_dir_is_doesnt_exist(dirpath=tmp_dir_path)
    preprocessed_input_dir_path = get_full_path(os.path.join(workdir_path, "input"))
    create_dir_is_doesnt_exist(dirpath=preprocessed_input_dir_path)
    if args.do_preprocess:
        if args.do_pre_nas_refine:
            pass
        if args.do_pre_scnt_refine:
            pass
        if args.do_pre_bnd:
            pass
    if args.do_run:
        if fragments is None:
            print("FRAGMENTS must exist; their absence at this point indicates that preprocessing was turned off, while --fragments were not part of the original input")
            exit(1)
        if scn_boundaries is None:
            print("Lower/Upper boundaries for Segment Copy Numbers must exist; "
                  "their absence at this point indicated that preprocessing was turned off, while --scn-boundaries were not part of the original input")
            exit(1)
        segments_to_fragments = None
        try:
            segments_to_fragments = get_segments_for_fragments_ids_dict(segments=segments, fragments=fragments, allow_non_covered=False)
        except Exception:
            print("FRAGMENTS boundaries do not align with SEGMENTS boundaries; "
                  "having this problem at this point indicates that the preprocessing was turned off and the data in original input SCNT and --fragments is not aligned")
            exit(1)
        segments_positions = []
        for s in segments:
            segments_positions.append(s.start_position)
            segments_positions.append(s.end_position)
        nas_positions = []
        for na in nas:
            nas_positions.append(na.position1)
            nas_positions.append(na.position2)
        if not positions_aligned(segments_positions=segments_positions, other_positions=nas_positions):
            print("NOVEL ADJACENCIES positions are not aligned with extremities of SEGMENTS; "
                  "having this problem at this point indicates that the preprocessing was turned off and the data in original input SCNT and NAS is not aligned")
            exit(1)
        if telomeres is None:
            telomeres = []
        telomeres.extend(get_ref_telomeres_from_segments(segments=segments))
        if not positions_aligned(segments_positions=segments_positions, other_positions=telomeres):
            print("TELOMERES positions are not aligned with extremities of SEGMENTS; "
                  "having this problem at this problem indicates that preprocessing was turned off and the data in the original SCNT and --telomeres is not aligned")
            exit(1)
        if nas_groups_ids is not None and not nas_groups_aligned(nas_groups_ids=nas_groups_ids, nas=nas):
            print("Some novel adjacencies group in NAS-GROUPS refers to an absent NOVEL ADJACENCY; "
                  "having this problem at this time indicates that either the original --nas-groups has a group that is referring to the missing novel adjacency in NAS, "
                  "or the referred novel adjacency from NAS was removed during the preprocessing.")
            exit(1)
        if nas_groups_ids is not None:
            nas_groups = construct_nas_groups(nas_groups_ids=nas_groups_ids, nas=nas)
            nas_groups_by_ids = {nas_group.stable_id_non_phased: nas_group for nas_group in nas_groups}
        else:
            nas_groups_by_ids = None
        ref_adjacencies = get_ref_adjacencies_from_segments(segments=segments, assign_external_ids=True)
        adjacencies = nas + ref_adjacencies
        if trees is None:
            trees = [None]
        else:
            trees = [None] + trees
        for tree_cnt, tree in enumerate(trees):
            tree_dir_basename = args.out_tree_dir_basename_template.format(tree_id=tree_cnt)
            tree_dir_path = os.path.join(workdir_path, tree_dir_basename)
            create_dir_is_doesnt_exist(dirpath=tree_dir_path)
            for nas_fp in nas_fps:
                instance_dir_path = os.path.join(tree_dir_path, str(nas_fp))
                create_dir_is_doesnt_exist(dirpath=instance_dir_path)
                gurobi_log_path = os.path.join(instance_dir_path, "gurobi.log")
                ilp_model = OptModelMultiClone(hapl_segments=segments,
                                               hapl_adjacencies=adjacencies,
                                               clone_specific_scnp=scnt,
                                               hapl_telomeres=telomeres,
                                               hapl_segments_to_fragments=segments_to_fragments,
                                               hapl_adjacencies_groups=nas_groups_by_ids,
                                               scn_boundaries=scn_boundaries,
                                               hapl_nov_adjacencies_fp=nas_fp)
                ilp_model.build_gurobi_model()
                ilp_model.gm.setParam("MIPGap", args.run_g_mip_gap)
                ilp_model.gm.setParam("MIPGapAbs", args.run_g_mip_gap)
                ilp_model.gm.setParam("LogFile", gurobi_log_path)
                ilp_model.gm.setParam("TimeLimit", args.run_g_time_limit)
                ilp_model.solve_model()

                status = ilp_model.gm.status
                if status == g.GRB.Status.OPTIMAL or status == g.GRB.Status.TIME_LIMIT:
                    inferred_scnt = ilp_model.get_clone_specific_scnp_from_model()
                    inferred_acnt = ilp_model.get_clone_specific_acnp_from_model()
                    scnt_base_name = args.out_scnt_basename_template.format(prefix=args.out_prefix_name, fp=nas_fp)
                    if scnt_base_name.startswith("."):
                        scnt_base_name = scnt_base_name[1:]
                    scnt_base_name = re.sub("\.+", ".", scnt_base_name)
                    scnt_file_path = os.path.join(instance_dir_path, scnt_base_name)
                    write_scn_tensor(file_name=scnt_file_path, scnt=inferred_scnt, segments=segments)
                    acnt_base_name = args.out_acnt_basename_template.format(prefix=args.out_prefix_name, fp=nas_fp)
                    if acnt_base_name.startswith("."):
                        acnt_base_name = acnt_base_name[1:]
                    acnt_base_name = re.sub("\.+", ".", acnt_base_name)
                    acnt_file_path = os.path.join(instance_dir_path, acnt_base_name)
                    write_acnt(file_name=acnt_file_path, acnt=inferred_acnt, adjacencies=adjacencies, output_reference=True)
                elif status == g.GRB.Status.INFEASIBLE:
                    ilp_model.gm.computeIIS()
                    ilp_path = os.path.join(instance_dir_path, "model.ilp")
                    ilp_model.gm.write(ilp_path)
