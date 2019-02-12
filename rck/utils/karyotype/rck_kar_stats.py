import argparse
from collections import defaultdict

from rck.core.graph import construct_hiag_inflate_from_haploid_data
from rck.core.io import read_scnt_from_source, read_acnt_from_source, read_scnb_from_source, read_adjacency_groups_from_source, read_positions_from_source, get_logging_cli_parser, \
    get_standard_logger_from_args, EXTERNAL_NA_ID
from rck.core.structures import get_ref_telomeres_from_segments, AdjacencyType, AdjacencyGroupType
from utils.karyotype.analysis import adjacency_groups_molecule_violations, adjacency_groups_labeling_violations, adjacency_groups_general_violations


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-KAR-stats", parents=[get_logging_cli_parser()])
    parser.add_argument("--verbose", choices=[0, 1, 2, 3, 4, 5], type=int, default=5)
    parser.add_argument("--acnt", required=True, type=argparse.FileType("rt"))
    parser.add_argument("--acnt-separator", default="\t")
    parser.add_argument("--acnt-extra-separator", default=";")
    parser.add_argument("--scnt", required=True, type=argparse.FileType("rt"))
    parser.add_argument("--scnt-separator", default="\t")
    parser.add_argument("--scnt-extra-separator", default=";")
    parser.add_argument("--scnb", type=argparse.FileType("rt"))
    parser.add_argument("--scnb-separator", default="\t")
    parser.add_argument("--scnb-extra-separator", default=";")
    parser.add_argument("--nas-fp", type=float, default=-1.0)
    parser.add_argument("--adjacency-groups", type=argparse.FileType("rt"))
    parser.add_argument("--adg-separator", default="\t")
    parser.add_argument("--adg-aids-separator", default=",")
    parser.add_argument("--adg-extra-separator", default=";")
    parser.add_argument("--telomere-positions", type=argparse.FileType("rt"))
    parser.add_argument("--telomere-positions-separator", default="\t")
    parser.add_argument("--telomere-positions-extra-separator", default=";")
    args = parser.parse_args()
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-KAR-stats")
    logger.info("Reading segment copy number tensor from {file}".format(file=args.scnt))
    segments, scnt = read_scnt_from_source(source=args.scnt, separator=args.scnt_separator, extra_separator=args.scnt_extra_separator, remove_cn_data_from_segs=True)
    logger.info("Reading adjacency copy number tensor from {file}".format(file=args.acnt))
    adjacencies, acnt = read_acnt_from_source(source=args.acnt, separator=args.acnt_separator, extra_separator=args.acnt_extra_separator, remove_cn_data_from_adj=True)
    if args.scnb is not None:
        logger.info("Reading segment copy number boundaries tensor from {file}".format(file=args.scnb))
        _, scnb = read_scnb_from_source(source=args.scnb, separator=args.scnb_separator, extra_separator=args.scnb_extra_separator, remove_cnb_data_from_segs=True)
    else:
        logger.info("No segment copy number boundaries tensor is provided via --scnb flag")
        scnb = None
    if args.adjacency_groups is not None:
        logger.info("Reading adjacency groups information from {file}".format(file=args.adjacency_groups))
        groups = read_adjacency_groups_from_source(source=args.adjacency_groups, separator=args.adg_separator,
                                                   extra_separator=args.adg_extra_separator, aids_separator=args.adg_aids_separator)
    else:
        logger.info("No adjacency groups information is provided via --adjacency-groups flag")
        groups = []
    if args.telomere_positions is not None:
        logger.info("Reading telomere positions from {file}".format(file=args.telomere_positions))
        telomeres = read_positions_from_source(source=args.telomere_positions, separator=args.telomeres_positions_separator,
                                               extra_separator=args.telomere_positions_extra_separator)
    else:
        logger.info("No telomere positions are provided via --telomere-positions flag. Defaulting to reference telomere positions".format(file=args.telomere_positions))
        telomeres = get_ref_telomeres_from_segments(segments=segments)
    segments_by_chrs = defaultdict(list)
    for segment in segments:
        segments_by_chrs[segment.chromosome].append(segment)
    print("A total of {cnt} chromosomes are observed".format(cnt=len(segments_by_chrs)))
    total_segments_cnt = 0
    for chr_name, chr_segments in segments_by_chrs.items():
        total_segments_cnt += len(chr_segments)
        if args.verbose >= 3:
            print("Chromosome {chr_name} has {cnt} segments".format(chr_name=chr_name, cnt=len(chr_segments)))
    print("A total of {cnt} segments are observed".format(cnt=total_segments_cnt))
    novel_adjacencies = [adj for adj in adjacencies if adj.adjacency_type == AdjacencyType.NOVEL]
    reference_adjacencies = [adj for adj in adjacencies if adj.adjacency_type == AdjacencyType.REFERENCE]
    print("A total of {cnt} adjacencies ({n_cnt} novel; {r_cnt} reference)".format(cnt=len(novel_adjacencies) + len(reference_adjacencies),
                                                                                   n_cnt=len(novel_adjacencies), r_cnt=len(reference_adjacencies)))

    adjacencies_by_external_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in adjacencies}
    if groups is not None:
        for ag in groups:
            ag.populate_adjacencies_via_ids(source=adjacencies, source_by_ids=adjacencies_by_external_ids)
        molecule_groups = [ag for ag in groups if ag.group_type == AdjacencyGroupType.MOLECULE]
        labeling_groups = [ag for ag in groups if ag.group_type == AdjacencyGroupType.LABELING]
        general_groups = [ag for ag in groups if ag.group_type == AdjacencyGroupType.GENERAL]
        if len(molecule_groups) > 0:
            logger.info("Checking compliance with {cnt} molecule groups".format(cnt=len(molecule_groups)))
            molecule_groups_violations = adjacency_groups_molecule_violations(groups=molecule_groups, acnt=acnt)
            if len(molecule_groups_violations):
                logger.error("A total of {cnt} molecule groups DO NOT agree with input karyotype. See molecule groups ids below".format(cnt=len(molecule_groups)))
                logger.error(", ".join([ag.gid for ag in molecule_groups_violations]))
            else:
                logger.info("All molecule groups agree with input karyotype")
        else:
            logger.info("No molecule groups were provided. Nothing to check.")
        if len(labeling_groups) > 0:
            logger.info("Checking compliance with {cnt} labeling groups".format(cnt=len(labeling_groups)))
            labeling_groups_violations = adjacency_groups_labeling_violations(groups=labeling_groups, acnt=acnt)
            if len(labeling_groups_violations):
                logger.error("A total of {cnt} labeling groups DO NOT agree with input karyotype. See labeling groups ids below".format(cnt=len(labeling_groups_violations)))
                logger.error(", ".join([ag.gid for ag in labeling_groups_violations]))
            else:
                logger.info("All labeling groups agree with input karyotype")
        else:
            logger.info("No labeling groups were provided. Nothing to check.")
        if len(general_groups) > 0:
            logger.info("Checking compliance with {cnt} general groups".format(cnt=len(general_groups)))
            general_groups_violations = adjacency_groups_general_violations(groups=general_groups, acnt=acnt)
            if len(general_groups_violations):
                logger.error("A total of {cnt} general groups DO NOT agree with input karyotype. See general groups ids below".format(cnt=len(general_groups_violations)))
                logger.error(", ".join([ag.gid for ag in general_groups_violations]))
            else:
                logger.info("All general groups agree with input karyotype")
    else:
        logger.info("No information about adjacency groups were provided. Nothing to check.")

    clone_ids = sorted(set(scnt.keys()) & set(acnt.keys()))
    for clone_id in clone_ids:
        logger.info("Checking balancing and telomeres for clone {clone_id}".format(clone_id=clone_id))
        hiag = construct_hiag_inflate_from_haploid_data(hapl_segments=segments, hapl_adjacencies=adjacencies)
        scnp = scnt[clone_id]
        acnp = acnt[clone_id]
        hiag.assign_copy_numbers_from_scn_profile(scn_profile=scnp)
        hiag.assign_copy_numbers_from_acn_profile(acn_profile=acnp)
        hiag.remove_edges_with_zero_cn()
        logger.info("Checking that every vertex has a copy number excess >= 0.")
        for node in hiag.nodes(data=False):
            if hiag.node_imbalance(node=node) < 0:
                logger.warning("Something went WRONG! On segment extremity {node} there is a negative copy number excess...".format(node=str(node)))
        logger.info("Getting inferred telomeres.")
        diploid_telomeres = hiag.get_telomeres()
        inferred_hapl_telomeres_ids = {p.stable_id_non_hap for p in diploid_telomeres}
        input_hapl_telomers_ids = {p.stable_id_non_hap for p in telomeres}
        if inferred_hapl_telomeres_ids > input_hapl_telomers_ids:
            logger.error("Something went WRONG! Following segments extremities, while not specified specified as possible telomere sites were inferred as such.")
            logger.error(",".join(map(str, sorted(inferred_hapl_telomeres_ids - input_hapl_telomers_ids))))
        else:
            logger.info("Everything is OK! in clone {clone_id} all extremities have non-negative copy number excess, and inferred telomere sites concur with the input"
                        "".format(clone_id=clone_id))


if __name__ == "__main__":
    main()
