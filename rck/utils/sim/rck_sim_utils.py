import argparse
import os
import pathlib

from rck.core.io import get_logging_cli_parser, read_acnt_from_file, read_scnt_from_file, write_adjacencies_to_file, write_scnt_to_file
from rck.utils.sim.rck_sim_run_prep import get_unlabeled_adjacencies, get_flipped_scnt, get_noisy_adjacencies, get_noisy_scnt, get_errored_scnt


def main():
    parser = argparse.ArgumentParser(prog="RCK-SIM-UTILS")
    subparsers = parser.add_subparsers(title="command", dest="command")
    cli_logging_parser = get_logging_cli_parser()
    subparsers.required = True
    #######
    run_prep_parser = subparsers.add_parser("prepare", parents=[cli_logging_parser])
    run_prep_parser.add_argument("-e", "--experiment", type=str, default=os.getcwd())
    run_prep_parser.add_argument("--clones", nargs="+")
    run_prep_parser.add_argument("--adj-noise", action="store_true", dest="adj_noise")
    run_prep_parser.add_argument("--adj-coord-n-prob", type=float, default=0.5)
    run_prep_parser.add_argument("--adj-coord-n-max", type=int, default=50)
    run_prep_parser.add_argument("--adj-fp", type=float, default=0.0)
    run_prep_parser.add_argument("--scnt-noise-average", action="store_true", dest="scnt_noise_average")
    run_prep_parser.add_argument("--scnt-noise-error", action="store_true", dest="scnt_noise_error")
    run_prep_parser.add_argument("--scnt-noise-error-p", type=float, default=0.03)
    run_prep_parser.add_argument("--scnt-noise-error-magnitude", type=int, default=1)
    run_prep_parser.add_argument("--scnt-chunk-size", type=int, default=50000)
    run_prep_parser.add_argument("-o", "--output", type=str, default="exp_prep")
    args = parser.parse_args()
    if args.command == "prepare":
        pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)
        clones = set(args.clones)
        adjacencies, acnt = read_acnt_from_file(file_name=os.path.join(args.experiment, "rck.acnt.tsv"))
        adjacencies = get_unlabeled_adjacencies(adjacencies=adjacencies, acnt=acnt, clones=clones)
        if args.adj_noise:
            adjacencies = get_noisy_adjacencies(adjacencies=adjacencies, coordinate_noise_prob=args.adj_coord_n_prob,
                                                coordinate_noise_max=args.adj_coord_n_max, fp_rate=args.adj_fp)
        segments, scnt = read_scnt_from_file(file_name=os.path.join(args.experiment, "rck.scnt.tsv"))
        scnt = {clone_id: scnp for clone_id, scnp in scnt.items() if clone_id in clones}
        if args.scnt_noise_average:
            segments, scnt = get_noisy_scnt(segments=segments, scnt=scnt, chunk_size=args.scnt_chunk_size)
        scnt = get_flipped_scnt(segments=segments, scnt=scnt)
        if args.scnt_noise_error:
            segments, scnt = get_errored_scnt(segments=segments, scnt=scnt, error_p=args.scnt_noise_error_p, error_magnitude=args.scnt_noise_error_magnitude)
        write_adjacencies_to_file(file_name=os.path.join(args.output, "rck.adj.tsv"), adjacencies=adjacencies)
        write_scnt_to_file(file_name=os.path.join(args.output, "rck.scnt.tsv"), segments=segments, scnt=scnt)


if __name__ == "__main__":
    main()
