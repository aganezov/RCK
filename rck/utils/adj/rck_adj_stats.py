import argparse
import itertools
import random
import re

import sys
import os
import math
from collections import Counter
import pandas as pd

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import write_adjacencies_to_destination, read_adjacencies_from_source
from rck.utils.adj.convert import *
from rck.utils.adj.process import filter_adjacencies_by_chromosomal_regions, get_shared_nas_parser, ORIGIN_IDS


def get_str_size_label(int_size):
    if int_size < 1000:
        return str(int_size)
    elif int_size < 1000000:
        return str(math.floor(int_size / 1000)) + "KB"
    return str(math.floor(int_size / 1000000)) + "MB"


def add_bar_values(ax, bars):
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., 1.01 * height,
                '%d' % int(height),
                ha='center', va='bottom', fontsize=20)


def get_length(na, treat_ins_separately=True):
    if treat_ins_separately and na.extra.get(SVTYPE, "").lower() == "ins":
        result = int(na.extra.get(SVLEN, -1))
        if result > -1:
            return result
    return na.distance_non_hap


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADJ-STATS")
    parser.add_argument('--version', action='version', version=rck.version)
    ######
    shared_parser = get_shared_nas_parser()
    shared_parser.add_argument("--no-vis", action="store_false", dest="vis")
    shared_parser.add_argument("--vis-interactive", action="store_true")
    shared_parser.add_argument("--output-dir", "-o", dest="output_dir", default=os.getcwd())
    ######
    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True
    ######
    cnt_parser = subparsers.add_parser("cnt", parents=[shared_parser], help="Counting stats for RCK NAS in input file")
    cnt_parser.add_argument("rck_nas", type=argparse.FileType("rt"), default=sys.stdin)
    cnt_parser.add_argument("--ann", action="store_true")
    cnt_parser.add_argument("--ann-list", nargs=1)
    cnt_parser.add_argument("--ann-field", default="svtype", type=str)
    cnt_parser.add_argument("--ann-missing", default="unknown")
    cnt_parser.add_argument("--bins", nargs=1,
                            default=",".join(["-1", "0", "50", "100", "500", "1000", "5000", "10000", "50000",
                                              "100000", "500000", "1000000", "5000000", "10000000"]),
                            # default=",".join(["-1", "0", "50", "100", "200", "300", "400", "500"])
                            )
    cnt_parser.add_argument("--cnt-output-subdir", default="cnt")
    cnt_parser.add_argument("--title", default="")
    cnt_parser.add_argument("--no-bar-values", action="store_false", dest="bar_values")
    # cnt_parser.add_argument("--per-chr", action="store_true")
    ######
    lr_parser = subparsers.add_parser("lr", parents=[shared_parser], help="Counting stats w.r.t. long read information in the input RCK NAS file")
    lr_parser.add_argument("rck_nas", type=argparse.FileType("rt"), default=sys.stdin)
    lr_parser.add_argument("--lr-field", default="support_read_names")
    lr_parser.add_argument("--title", default="")
    lr_parser.add_argument("--no-bar-values", action="store_false", dest="bar_values")
    lr_parser.add_argument("--lr-output-subdir", default="lr")
    ######
    merged_parser = subparsers.add_parser("merged", parents=[shared_parser], help="Counting statistics over merged RCK NAS")
    merged_parser.add_argument("rck_nas", type=argparse.FileType("rt"), default=sys.stdin)
    merged_parser.add_argument("--origin-field", default=ORIGIN_IDS)
    merged_parser.add_argument("--origin-sep", default=",")
    merged_parser.add_argument("--origin-regex", default=".*_(?P<source>.*)")
    merged_parser.add_argument("--no-origin-field", choices=["skip", "self"])
    merged_parser.add_argument("---merged-output-subdir", default="merged")
    merged_parser.add_argument("--title", default="")
    #######
    support_parser = subparsers.add_parser("support", parents=[shared_parser], help="Counting support statistics over merged RCK NAS")
    support_parser.add_argument("rck_nas", type=argparse.FileType("rt"), default=sys.stdin)
    support_parser.add_argument("--sources", type=argparse.FileType("rt"), nargs="+")
    support_parser.add_argument("--title", default="")
    support_parser.add_argument("--no-bar-values", action="store_false", dest="bar_values")
    support_parser.add_argument("--support-output-subdir", default="support")
    #######
    args = parser.parse_args()
    if args.vis:
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.set(color_codes=True)
    if args.command == "cnt":
        nas = read_adjacencies_from_source(source=args.rck_nas)
        use_annotations = args.ann
        bins = set()
        for bin_str in args.bins.split(","):
            bin_value = int(bin_str)
            bins.add(bin_value)
        bins = sorted(bins)
        if bins[-1] < 500000000:
            bins.append(500000000)
        if use_annotations:
            annotations = set()
            if len(args.ann_list) == 1 and args.ann_list[0] == "all":
                allow_all_annotations = True
            else:
                allow_all_annotations = False
            for anns_str in args.ann_list:
                anns = [ann.lower() for ann in anns_str.split(",")]
                for ann in anns:
                    annotations.add(ann)
            annotations = sorted(annotations)
            nas_by_anns = defaultdict(list)
            for na in nas:
                na_ann = na.extra.get(args.ann_field, args.ann_missing).lower()
                if allow_all_annotations or na_ann in annotations:
                    nas_by_anns[na_ann].append(na)
        else:
            nas_by_anns = defaultdict(list)
            for na in nas:
                na_ann = str(na.position1.strand) + str(na.position2.strand)
                nas_by_anns[na_ann].append(na)
            # separate figures
        for na_ann, nas in nas_by_anns.items():
            lengths = [get_length(na=na) for na in nas]
            bin_cnts = defaultdict(int)
            for length in lengths:
                for l, r in zip(bins[:-1], bins[1:]):
                    if length < r:
                        bin_cnts[l] += 1
                        break

            if args.vis:
                values = [bin_cnts[i] for i in bins]
                x_axis_values = [i for i in range(len(bins))]
                plt.figure(figsize=(20, 10))
                bars = plt.bar(x_axis_values, values, label=na_ann, color="g")
                if args.bar_values:
                    add_bar_values(ax=plt.gca(), bars=bars)
                x_axis_values = [i for i in range(len(bins))]
                plt.xticks(x_axis_values, ["[{l}-\n{r})".format(l=get_str_size_label(int_size=l), r=get_str_size_label(int_size=r)) for l, r in zip(bins[:-1], bins[1:])])
                plt.tick_params(labelsize=20)
                plt.legend(prop={'size': 20})
                plt.xlabel("{ann} lengths".format(ann=na_ann), fontsize=20)
                plt.ylabel("# of {ann}".format(ann=na_ann), fontsize=20)
                plt.title(args.title, fontsize=20)
                # plt.xticks(bins, [str(b) for b in bins])
                plot_path = os.path.abspath(os.path.join(args.output_dir, args.cnt_output_subdir))
                if not os.path.exists(plot_path):
                    os.makedirs(plot_path)
                plot_file_name = os.path.join(plot_path, "{prefix}_cnt.png".format(prefix=na_ann.replace(os.path.sep, "_")))
                plt.savefig(plot_file_name)
                if args.vis_interactive:
                    plt.show()
                plt.clf()
        # one figure
        fig, ax = plt.subplots()
        fig.set_figheight(20)
        fig.set_figwidth(20)
        diff = .9 / len(list(nas_by_anns.keys()))
        if len(list(nas_by_anns.keys())) % 2 == 0:
            extra_diff = diff / 2
        else:
            extra_diff = 0
        x_axis_values = [i for i in range(len(bins))]
        colors = ["b", "g", "r", "c", "m", "y", "b", "w"]
        for cnt, (na_ann, nas) in enumerate(nas_by_anns.items()):
            lengths = [get_length(na=na) for na in nas]
            bin_cnts = defaultdict(int)
            for length in lengths:
                for l, r in zip(bins[:-1], bins[1:]):
                    if length < r:
                        bin_cnts[l] += 1
                        break
            if args.vis:
                values = [bin_cnts[i] for i in bins]
                x_axis_values_tmp = [x + diff * cnt for x in x_axis_values]
                color = colors[cnt]
                ax.bar(x_axis_values_tmp, values, width=diff, color=color, label="{ann}".format(ann=na_ann))
        plt.xticks(x_axis_values, ["[{l}-\n{r})".format(l=get_str_size_label(int_size=l), r=get_str_size_label(int_size=r)) for l, r in zip(bins[:-1], bins[1:])])
        plt.tick_params(labelsize=20)
        plt.legend(prop={'size': 20})
        plt.xlabel("lenghts", fontsize=20)
        plt.ylabel("# of SVs", fontsize=20)
        plt.title(args.title, fontsize=20)
        plot_path = os.path.abspath(os.path.join(args.output_dir, args.cnt_output_subdir))
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        plot_file_name = os.path.join(plot_path, "bins_cnt.png")
        plt.savefig(plot_file_name)
        if args.vis_interactive:
            plt.show()
        plt.gcf()

        fig, ax = plt.subplots()
        fig.set_figheight(20)
        fig.set_figwidth(20)
        x_axis_values = {ann: cnt for cnt, ann in enumerate(nas_by_anns.keys())}
        values = {ann: len(nas_by_anns[ann]) for ann in nas_by_anns.keys()}
        x_axis_values = [x_axis_values[ann] for ann in sorted(nas_by_anns.keys())]
        values = [values[ann] for ann in sorted(nas_by_anns.keys())]
        x_ticks = [ann for ann in sorted(nas_by_anns.keys())]
        bars = plt.bar(x_axis_values, values, color="g")
        if args.bar_values:
            add_bar_values(ax=ax, bars=bars)
        plt.xticks(x_axis_values, x_ticks)
        plt.tick_params(labelsize=20)
        plt.xlabel("SVs", fontsize=20)
        plt.ylabel("# of SVs", fontsize=20)
        plt.title(args.title, fontsize=20)
        plot_path = os.path.abspath(os.path.join(args.output_dir, args.cnt_output_subdir))
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        plot_file_name = os.path.join(plot_path, "cnt.png")
        plt.savefig(plot_file_name)
        if args.vis_interactive:
            plt.show()
        plt.gcf()
    elif args.command == "lr":
        nas = read_adjacencies_from_source(source=args.rck_nas)
        reads_to_nas = defaultdict(list)
        for na in nas:
            reads_str = na.extra.get(args.lr_field, "")
            reads = reads_str.split(",")
            for read in reads:
                if len(read) == 0:
                    continue
                reads_to_nas[read].append(na)
        nas_cnts = {len(reads_to_nas[read]) for read in reads_to_nas.keys()}
        x_axis_values = [i for i in range(2, max(nas_cnts))]
        values = defaultdict(int)
        for read in reads_to_nas.keys():
            values[len(reads_to_nas[read])] += 1
        values = [values[i] for i in x_axis_values]
        x_ticks = sorted(x_axis_values)
        fig, ax = plt.subplots()
        fig.set_figheight(20)
        fig.set_figwidth(20)
        bars = plt.bar(x_axis_values, values, color="g")
        if args.bar_values:
            add_bar_values(ax=ax, bars=bars)
        plt.xticks(x_axis_values, x_ticks)
        plt.xlabel("# of SVs")
        plt.ylabel("# of reads supporting x SVs")
        plt.title(args.title)
        plot_path = os.path.abspath(os.path.join(args.output_dir, args.lr_output_subdir))
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        plot_file_name = os.path.join(plot_path, "lr_cnt.png")
        plt.savefig(plot_file_name)
        if args.vis_interactive:
            plt.show()
        plt.gcf()
    elif args.command == "merged":
        nas = read_adjacencies_from_source(source=args.rck_nas)
        source_pattern = re.compile(args.origin_regex)
        source_groups_cnt = defaultdict(int)
        for na in nas:
            naid = na.extra.get(EXTERNAL_NA_ID, na.idx)
            if args.origin_field not in na.extra:
                if args.no_origin_field == "skip":
                    continue
                elif args.no_origin_field == "self":
                    origin = naid
                else:
                    raise Exception("Unknown strategy {no_origin_field} for Adjacency {naid} that misses the origin field"
                                    "".format(no_origin_field=args.no_origin_field, naid=naid))
            else:
                origin = na.extra[args.origin_field]
            origin_strings = origin.split(args.origin_sep)
            origins = []
            for origin_string in origin_strings:
                origin_match = source_pattern.match(origin_string)
                if origin_match is None:
                    continue
                origins.append(origin_match.group("source"))
            if len(origins) == 0:
                continue
            origins = tuple(sorted(origins))
            source_groups_cnt[origins] += 1
        source_group_set_cnt = defaultdict(int)
        groups = set()
        pairwise_cnts = defaultdict(int)
        print("\n----- Quantitative subgroups ----\n")
        for key, value in sorted(source_groups_cnt.items(), key=lambda entry: entry[1]):
            for element in key:
                groups.add(element)
            counter = Counter(key)
            if len(set(key)) == 1:
                group = set(key).pop()
                pairwise_cnts[(group, group)] += value
            else:
                for g1, g2 in itertools.combinations(sorted(set(key)), r=2):
                    pairwise_cnts[tuple(sorted([g1, g2]))] += value
                    pairwise_cnts[tuple(reversed(sorted([g1, g2])))] += value
            source_group_set_cnt[tuple(sorted(counter.keys()))] += value
            counter_str = ", ".join(["{name} ({cnt})".format(name=name, cnt=cnt) for name, cnt in sorted(counter.items(), key=lambda entry: (entry[1], entry[0]))])
            print("{key} :: {value}".format(key=counter_str, value=value))

        print("\n----- Groups ----- \n")
        print(", ".join(sorted(groups)))

        print("\n----- By group size (then by size)------\n")
        for key, value in sorted(source_group_set_cnt.items(), key=lambda entry: (len(entry[0]), entry[1])):
            print("{key} :: {value}".format(key=",".join(key), value=value))

        print("\n----- By size (then by group size)------\n")
        for key, value in sorted(source_group_set_cnt.items(), key=lambda entry: (entry[1], len(entry[0]))):
            print("{key} :: {value}".format(key=",".join(key), value=value))

        print("\n----- Pairsize ------\n")
        for key, value in sorted(pairwise_cnts.items()):
            if key[0] >= key[1]:
                print("{key} :: {value}".format(key=",".join(key), value=value))

        if args.vis:
            ser = pd.Series(list(pairwise_cnts.values()),
                            index=pd.MultiIndex.from_tuples(pairwise_cnts.keys()))
            df = ser.unstack().fillna(0)
            df = df.astype('int64')
            sns.set(font_scale=1.5)
            sns.heatmap(df, annot=True, fmt="d")
            plt.title(args.title)
            plot_path = os.path.abspath(os.path.join(args.output_dir, args.merged_output_subdir))
            if not os.path.exists(plot_path):
                os.makedirs(plot_path)
            plot_file_name = os.path.join(plot_path, "pairwise_merged_cnts.png")

            plt.savefig(plot_file_name, bbox_inches='tight')
        if args.vis_interactive:
            plt.show()
    elif args.command == "support":
        max_support = 30
        import statistics
        nas = read_adjacencies_from_source(source=args.rck_nas)
        source_nas = []
        if args.sources is not None:
            for source in args.sources:
                source_nas.extend(read_adjacencies_from_source(source=source))
        source_nas_by_ids = {na.extra.get(EXTERNAL_NA_ID, na.idx): na for na in source_nas}
        counts = defaultdict(int)
        nas_read_cnt = defaultdict(int)
        for na in nas:
            if "support_read_names" in na.extra:
                support_cnt = len(na.extra["support_read_names"].split(","))
                if support_cnt > max_support:
                    support_cnt = max_support
                counts[support_cnt] += 1
                nas_read_cnt[na] += support_cnt
            elif "origin_ids" in na.extra:
                origin_nas = [source_nas_by_ids[origin_id] for origin_id in na.extra["origin_ids"].split(",") if origin_id in source_nas_by_ids]
                if len(origin_nas) == 0:
                    continue
                supports = []
                for origin_na in origin_nas:
                    if "support_read_names" in origin_na.extra:
                        supports.append(len(origin_na.extra["support_read_names"].split(",")))
                support_cnt = int(statistics.mean(supports))
                if support_cnt > max_support:
                    support_cnt = max_support
                counts[support_cnt] += 1
                nas_read_cnt[na] += support_cnt
            else:
                counts[-1] += 1
        x_axis_values = [i for i in range(1, max(counts) + 1)]
        values = [counts[i] for i in x_axis_values]
        x_ticks = sorted(x_axis_values)
        fig, ax = plt.subplots()
        fig.set_figheight(20)
        fig.set_figwidth(20)
        bars = plt.bar(x_axis_values, values, color="g")
        if args.bar_values:
            add_bar_values(ax=ax, bars=bars)
        plt.xticks(x_axis_values, x_ticks)
        plt.xlabel("# of supporting reads")
        plt.ylabel("# of SVs")
        plt.title(args.title)
        plot_path = os.path.abspath(os.path.join(args.output_dir, args.support_output_subdir))
        if not os.path.exists(plot_path):
            os.makedirs(plot_path)
        plot_file_name = os.path.join(plot_path, "support_cnt.png")
        plt.savefig(plot_file_name)
        if args.vis_interactive:
            plt.show()
        plt.gcf()
        top_svs = {na: cnt for na, cnt in nas_read_cnt.items() if cnt >= 28 and na.distance_non_hap > 300}
        result = random.sample(top_svs.keys(), 10)
        for entry in result:
            print(entry, entry.distance_non_hap)


if __name__ == "__main__":
    main()
