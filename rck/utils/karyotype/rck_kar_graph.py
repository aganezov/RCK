import argparse
import sys

from rck.core.graph import construct_hiag_inflate_from_haploid_data
from rck.core.io import get_logging_cli_parser, read_scnt_from_source, read_acnt_from_source, write_graph_to_destination


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-KAR-graph")
    cli_logging_parser = get_logging_cli_parser()
    parser.add_argument("--acnt", required=True, type=argparse.FileType("rt"))
    parser.add_argument("--acnt-separator", default="\t")
    parser.add_argument("--acnt-extra-separator", default=";")
    parser.add_argument("--scnt", required=True, type=argparse.FileType("rt"))
    parser.add_argument("--scnt-separator", default="\t")
    parser.add_argument("--scnt-extra-separator", default=";")
    parser.add_argument("--clone")
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    writer_parser = subparsers.add_parser("write", parents=[cli_logging_parser])
    writer_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    writer_parser.add_argument("--style", choices=["edge-list"], default="edge-list")
    writer_parser.add_argument("--separator", default="\t")
    writer_parser.add_argument("--include-absent", action="store_true", dest="include_cn_0")
    args = parser.parse_args()
    segments, scnt = read_scnt_from_source(source=args.scnt, separator=args.scnt_separator, extra_separator=args.scnt_extra_separator, remove_cn_data_from_segs=True)
    adjacencies, acnt = read_acnt_from_source(source=args.acnt, separator=args.acnt_separator, extra_separator=args.acnt_extra_separator, remove_cn_data_from_adj=True)
    if args.command == "write":
        hiag = construct_hiag_inflate_from_haploid_data(hapl_segments=segments, hapl_adjacencies=adjacencies)
        if args.clone is None:
            common_clones = set(acnt.keys()) & set(scnt.keys())
            if len(common_clones) == 0:
                raise ValueError("No common clones in Adjacency and Segment Copy Number tensors")
            args.clone_id = sorted(common_clones)[0]
        acnp, scnp = acnt[args.clone_id], scnt[args.clone_id]
        hiag.assign_copy_numbers_from_scn_profile(scn_profile=scnp)
        hiag.assign_copy_numbers_from_acn_profile(acn_profile=acnp)
        if not args.include_cn_0:
            hiag.remove_edges_with_zero_cn()
        write_graph_to_destination(graph=hiag, destination=args.output, style=args.style)


if __name__ == "__main__":
    main()
