import argparse
import sys

from rck.core.io import get_logging_cli_parser, read_scnt_from_source
from rck.utils.scn.stats import cn_distance, cn_change, changed_segments, T1, T2, changed_ampl_segments, changed_neutral_segments, changed_loss_segments, changed_haploid_segments


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-SCNT-stats")
    cli_logging_parser = get_logging_cli_parser()
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    #####
    distance_parser = subparsers.add_parser("distance", parents=[cli_logging_parser])
    distance_parser.add_argument('--scnt1', type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt1-separator", default="\t")
    distance_parser.add_argument("--scnt1-extra-separator", default=";")
    distance_parser.add_argument("--scnt1-clone-ids", default=None)
    distance_parser.add_argument('--scnt2', type=argparse.FileType("rt"), required=True)
    distance_parser.add_argument("--scnt2-separator", default="\t")
    distance_parser.add_argument("--scnt2-extra-separator", default=";")
    distance_parser.add_argument("--scnt2-clone-ids", default=None)
    distance_parser.add_argument("--topn", type=int, default=3)
    distance_parser.add_argument("--verbose", action="store_true", dest="verbose")
    distance_parser.add_argument("--both-haplotype-specific", action="store_true", dest="both_haplotype_specific")
    distance_parser.add_argument('-o', '--output', type=argparse.FileType("wt"), default=sys.stdout)
    #####
    cn_change_parser = subparsers.add_parser("cn-change", parents=[cli_logging_parser])
    cn_change_parser.add_argument("--scnt1", type=argparse.FileType("rt"), required=True)
    cn_change_parser.add_argument("--scnt1-separator", default="\t")
    cn_change_parser.add_argument("--scnt1-extra-separator", default=";")
    cn_change_parser.add_argument("--scnt1-clone-ids", default=None)
    cn_change_parser.add_argument("--scnt2", type=argparse.FileType("rt"), required=True)
    cn_change_parser.add_argument("--scnt2-separator", default="\t")
    cn_change_parser.add_argument("--scnt2-extra-separator", default=";")
    cn_change_parser.add_argument("--scnt2-clone-ids", default=None)
    cn_change_parser.add_argument("--topn", type=int, default=3)
    cn_change_parser.add_argument("--both-haplotype-specific", action="store_true", dest="both_haplotype_specific")
    cn_change_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    #####
    args = parser.parse_args()
    if args.command == "distance":
        scnt1_clone_ids = args.scnt1_clone_ids if args.scnt1_clone_ids is None else args.scnt1_clone_ids.split(",")
        segments1, scnt1 = read_scnt_from_source(source=args.scnt1, separator=args.scnt1_separator, extra_separator=args.scnt1_extra_separator, clone_ids=scnt1_clone_ids)
        scnt2_clone_ids = args.scnt2_clone_ids if args.scnt2_clone_ids is None else args.scnt2_clone_ids.split(",")
        segments2, scnt2 = read_scnt_from_source(source=args.scnt2, separator=args.scnt2_separator, extra_separator=args.scnt2_extra_separator, clone_ids=scnt2_clone_ids)
        result = cn_distance(segments1=segments1, scnt1=scnt1, segments2=segments2, scnt2=scnt2, both_haplotype_specific=args.both_haplotype_specific)
        sorted_result = sorted([(key, value) for key, value in result.items()], key=lambda entry: sum(entry[1].values()))
        output_result = sorted_result[:args.topn]
        if args.verbose:
            print(f'Length-weighted segment copy number distance for tensors in {args.scnt1.name} and {args.scnt2.name}', file=args.output)
        for cnt, (case, clone_specific_distance) in enumerate(output_result, start=1):
            print(f'{cnt}. Best distance (total) of {sum(clone_specific_distance.values()):,} with clone-specific ones {clone_specific_distance}, for case {case}', file=args.output)
    elif args.command == "cn-change":
        scnt1_clone_ids = args.scnt1_clone_ids if args.scnt1_clone_ids is None else args.scnt1_clone_ids.split(",")
        segments1, scnt1 = read_scnt_from_source(source=args.scnt1, separator=args.scnt1_separator, extra_separator=args.scnt1_extra_separator, clone_ids=scnt1_clone_ids)
        scnt2_clone_ids = args.scnt2_clone_ids if args.scnt2_clone_ids is None else args.scnt2_clone_ids.split(",")
        segments2, scnt2 = read_scnt_from_source(source=args.scnt2, separator=args.scnt2_separator, extra_separator=args.scnt2_extra_separator, clone_ids=scnt2_clone_ids)
        result = cn_change(segments1=segments1, scnt1=scnt1, segments2=segments2, scnt2=scnt2, both_haplotype_specific=args.both_haplotype_specific)
        sorted_result = sorted([(key, value) for key, value in result.items()], key=lambda entry: sum(entry[1][0].values()))
        output_result = sorted_result[:args.topn]
        for cnt, (case, (clone_specific_distance, clone_specific_cn_change)) in enumerate(output_result, start=1):
            print(f'{cnt}. Best distance (total) of {sum(clone_specific_distance.values()):,} with clone-specific ones {clone_specific_distance}, for case {case}', file=args.output)
            for clone_id, t1_t2_scnp in clone_specific_cn_change.items():
                cnt_total = len(segments1) * 2
                length_total = sum(s.length for s in segments1) * 2
                ch_haploid_segments = changed_haploid_segments(scnp1=t1_t2_scnp[T1], scnp2=t1_t2_scnp[T2], segments=segments1)
                cnt_ch_haploid = len(ch_haploid_segments)
                length_haploid_changed = sum(s.length for s in ch_haploid_segments)
                print(f'\tClone {clone_id}. '
                      f'Changed haploid segments cnt : {cnt_ch_haploid} / {cnt_total / 2} ({cnt_ch_haploid * 1.0 / (cnt_total / 2):0.2f}). '
                      f'length: {length_haploid_changed} / {length_total / 2} ({length_haploid_changed * 1.0 / (length_total / 2) :0.2f})')
                ch_segments = changed_segments(scnp1=t1_t2_scnp[T1], scnp2=t1_t2_scnp[T2], segments=segments1)
                cnt_changed = len(ch_segments)
                length_changed = sum(s.length for s in ch_segments)
                print(f'\tClone {clone_id}. '
                      f'Changed segments cnt : {cnt_changed} / {cnt_total} ({cnt_changed * 1.0 / cnt_total:0.2f}). '
                      f'length: {length_changed} / {length_total} ({length_changed * 1.0 / length_total:0.2f})')
                ch_ampl_segments = changed_ampl_segments(scnp1=t1_t2_scnp[T1], scnp2=t1_t2_scnp[T2], segments=segments1)
                cnt_ampl_changed = len(ch_ampl_segments)
                length_ampl_changed = sum(s.length for s in ch_ampl_segments)
                print(f'\t\tClone {clone_id}. '
                      f'Changed ampl segments cnt : {cnt_ampl_changed} / {cnt_total} ({cnt_ampl_changed * 1.0 / cnt_total:0.2f}). '
                      f'length: {length_ampl_changed} / {length_total} ({length_ampl_changed * 1.0 / length_total:0.2f})')
                ch_neutral_segments = changed_neutral_segments(scnp1=t1_t2_scnp[T1], scnp2=t1_t2_scnp[T2], segments=segments1)
                cnt_neutral_changed = len(ch_neutral_segments)
                length_neutral_changed = sum(s.length for s in ch_neutral_segments)
                print(f'\t\tClone {clone_id}. '
                      f'Changed neutral segments cnt : {cnt_neutral_changed} / {cnt_total} ({cnt_neutral_changed * 1.0 / cnt_total:0.2f}). '
                      f'length: {length_neutral_changed} / {length_total} ({length_neutral_changed * 1.0 / length_total:0.2f})')
                ch_del_segments = changed_loss_segments(scnp1=t1_t2_scnp[T1], scnp2=t1_t2_scnp[T2], segments=segments1)
                cnt_del_changed = len(ch_del_segments)
                length_dell_changed = sum(s.length for s in ch_del_segments)
                print(f'\t\tClone {clone_id}. '
                      f'Changed del segments cnt : {cnt_del_changed} / {cnt_total} ({cnt_del_changed * 1.0 / cnt_total:0.2f}). '
                      f'length: {length_dell_changed} / {length_total} ({length_dell_changed * 1.0 / length_total:0.2f})')


if __name__ == "__main__":
    main()
