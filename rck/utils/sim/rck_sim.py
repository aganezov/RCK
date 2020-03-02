import argparse
import os
import pathlib
from typing import Dict

import yaml

from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args, write_scnt_to_destination, write_acnt_to_destination
from rck.core.structures import SegmentCopyNumberProfile, AdjacencyCopyNumberProfile
from rck.utils.sim.manager import TREE, REFERENCE_FILE, SimManager
from rck.utils.sim.rck_sim_check import sim_check_as_strs


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SIMULATOR", parents=[get_logging_cli_parser()])
    parser.add_argument("setup", type=argparse.FileType("rt"))
    parser.add_argument("--tree", type=str, required=False)
    parser.add_argument("--ref", type=str, required=False)
    parser.add_argument("--store-intermediates", action="store_true", dest="store_inter_genomes")
    parser.add_argument("--no-post-check", action="store_false", dest="post_check")
    parser.add_argument("-o", "--output-dir", default="rck-sim-output")
    args = parser.parse_args()
    output_dir = args.output_dir
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if args.log_file is None:
        args.log_file = os.path.join(output_dir, "rck.sim.log")
    logger = get_standard_logger_from_args(args, "RCK-SIM")
    logger.info(f"Loading setting from file {args.setup.name}")
    settings = yaml.load(args.setup, Loader=yaml.FullLoader)
    if args.tree is not None:
        logger.debug(f"Phylo tree was not provided as command line argument, trying to get it from the setup config")
        settings[TREE] = args.tree
    if args.ref is not None:
        logger.debug(f"Reference was not provided as command line argument, trying to get it from the setup config")
        settings[REFERENCE_FILE] = args.ref
    logger.info(f"Creating simulation manager")
    sim_manager = SimManager.from_yaml_dict(settings, logger=logger)
    logger.info(f"Starting simulation")
    sim_manager.simulate(refine_after=True, store_intermediate_genomes=args.store_inter_genomes)
    logger.info(f"Success. Completed simulation")
    logger.info(f"Outputting simulation results")
    for genome_name, genome in sim_manager.genomes.items():
        with open(os.path.join(output_dir, f"raw_{genome_name}.rck.genome.txt"), "wt") as dest:
            logger.info(f"Outputting simulated (unrefined) genome {genome_name} to {dest.name}")
            print(str(genome), file=dest)
    scnt: Dict[str, SegmentCopyNumberProfile] = {}
    acnt: Dict[str, AdjacencyCopyNumberProfile] = {}
    all_adjacencies: set = set()
    for genome_name, genome in sim_manager.refined_genomes.items():
        genome.propagate_hap_labeles_to_positions()
        scnp: SegmentCopyNumberProfile = SegmentCopyNumberProfile.from_genome(genome)
        scnt[genome_name] = scnp
        acnp: AdjacencyCopyNumberProfile = AdjacencyCopyNumberProfile.from_genome(genome)
        acnt[genome_name] = acnp
        for adj in genome.iter_adjacencies():
            all_adjacencies.add(adj.get_non_phased_copy())
        with open(os.path.join(output_dir, f"{genome_name}.rck.genome.txt"), "wt") as dest:
            logger.info(f"Outputting simulated (refined) genome {genome_name} to {dest.name}")
            print(str(genome), file=dest)
    with open(os.path.join(output_dir, "true_rearr_history.yaml"), "wt") as dest:
        logger.info(f"Outputting exact rearrangement history to {dest.name}")
        yaml.dump(sim_manager.as_yaml(), stream=dest)
    with open(os.path.join(output_dir, "scnt.rck.tsv"), "wt") as dest:
        logger.info(f"Outputting true segment copy number tensor to {dest.name}")
        write_scnt_to_destination(destination=dest, segments={s.get_non_hap_copy() for s in sim_manager.refined_genomes[sim_manager.root_name].iter_segments()},
                                  scnt=scnt)
    with open(os.path.join(output_dir, "acnt.rck.tsv"), "wt") as dest:
        logger.info(f"Outputting true adjacency copy number tensor to {dest.name}")
        write_acnt_to_destination(destination=dest, acnt=acnt, adjacencies=all_adjacencies)

    if args.post_check:
        logger.info(f"Performing post simulation check to ensure results correctness")
        check_result_strs = sim_check_as_strs(ref_genome=sim_manager.refined_genomes[sim_manager.root_name],
                                              rearr_genomes={name: genome for name, genome in sim_manager.refined_genomes.items() if name != sim_manager.root_name},
                                              ref_genome_name=sim_manager.root_name)
        with open(os.path.join(output_dir, "post_check.txt"), "wt") as dest:
            logger.info(f"Outputting post simulation check results to {dest.name}")
            print("\n".join(check_result_strs), file=dest)


if __name__ == "__main__":
    main()
