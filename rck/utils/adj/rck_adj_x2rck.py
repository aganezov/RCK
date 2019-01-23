import argparse

import sys
import os

current_file_level = 3
current_dir = os.path.dirname(os.path.realpath(__file__))
for _ in range(current_file_level):
    current_dir = os.path.dirname(current_dir)
sys.path.append(current_dir)

import rck
from rck.core.io import get_logging_cli_parser, get_standard_logger_from_args
from rck.utils.adj.convert import *
from rck.utils.adj.process import filter_adjacencies_by_chromosomal_regions, get_shared_nas_parser, processed_gundem2015_adjacencies, get_chromosome_strip_parser


def main():
    parser = argparse.ArgumentParser(prog="RCK-UTILS-ADJ-x2rck")
    parser.add_argument('--version', action='version', version=rck.version)
    ####
    shared_parser = get_shared_nas_parser()
    shared_parser.add_argument("--output", "-o", dest="rck_adj_file", type=argparse.FileType("wt"), default=sys.stdout)
    cli_logging_parser = get_logging_cli_parser()
    chr_strip_parser = get_chromosome_strip_parser()
    ####
    subparsers = parser.add_subparsers(title="commands", dest="command")
    subparsers.required = True
    ####
    lumpy_parser = subparsers.add_parser("lumpy", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert Lumpy VCF SV calls into RCK NAS format")
    lumpy_parser.add_argument("--id-suffix", dest="id_suffix", default="lumpy")
    lumpy_parser.add_argument("lumpy_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    longranger_parser = subparsers.add_parser("longranger", parents=[shared_parser, cli_logging_parser, chr_strip_parser],
                                              help="Convert LongRanger VCF SV calls into RCK NAS format")
    longranger_parser.add_argument("--id-suffix", dest="id_suffix", default="longranger")
    longranger_parser.add_argument("longranger_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    naibr_parser = subparsers.add_parser("naibr", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert NAIBR NAS calls into RCK NAS format")
    naibr_parser.add_argument("--id-suffix", dest="id_suffix", default="naibr")
    naibr_parser.add_argument("naibr_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    manta_parser = subparsers.add_parser("manta", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert Manta VCF SV calls into RCK NAS format")
    manta_parser.add_argument("--id-suffix", dest="id_suffix", default="manta")
    manta_parser.add_argument("manta_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    sniffles_parser = subparsers.add_parser("sniffles", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert Sniffles VCF SV calls into RCK NAS format")
    sniffles_parser.add_argument("--id-suffix", dest="id_suffix", default="sniffles")
    sniffles_parser.add_argument("sniffles_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    grocsv = subparsers.add_parser("grocsv", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert GROCSVS VCF SV calls into RCK NAS format")
    grocsv.add_argument("--id-suffix", dest="id_suffix", default="grocsv")
    grocsv.add_argument("grocsv_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    grocsv.add_argument("--samples")
    grocsv.add_argument("--samples-all-any", choices=["all", "any"], default="any")
    grocsv.add_argument("--samples-only", action="store_true", dest="samples_only")
    ####
    delly = subparsers.add_parser("delly", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert Delly VCF SV calls into RCK NAS format")
    delly.add_argument("--id-suffix", dest="id_suffix", default="delly")
    delly.add_argument("delly_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    delly.add_argument("--stream", action="store_true", dest="delly_force_stream")
    ####
    pbsv = subparsers.add_parser("pbsv", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert PBSV VCF SV calls into RCK NAS format")
    pbsv.add_argument("--id-suffix", dest="id_suffix", default="pbsv")
    pbsv.add_argument("pbsv_vcf_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    remixt = subparsers.add_parser("remixt", parents=[shared_parser, cli_logging_parser, chr_strip_parser], help="Convert ReMixT Novel adjacencies calls into RCK format")
    remixt.add_argument("--i-separator", default="\t")
    remixt.add_argument("--id-suffix", dest="id_suffix", default="remixt")
    remixt.add_argument("--clone-ids", choices=["1", "2", "1,2"], default="1,2")
    remixt.add_argument("--skip-absent", action="store_true", dest="skip_absent")
    remixt.add_argument("--no-remixt-na-correction", action="store_false", dest="remixt_correction")
    remixt.add_argument("remixt_file", type=argparse.FileType("rt"), default=sys.stdin)
    ####
    gundem2015_parser = subparsers.add_parser("gundem2015", parents=[shared_parser, cli_logging_parser, chr_strip_parser],
                                              help="Convert SV calls from Gundem et al (2015) (BRASS2???) "
                                                   "into RCK NAS format")
    gundem2015_parser.add_argument("--id-suffix", dest="id_suffix", default="gundem2015")
    gundem2015_parser.add_argument("gundem2015_file", type=argparse.FileType("rt"), default=sys.stdin)
    gundem2015_parser.add_argument("--i-separator", default="\t")
    gundem2015_parser.add_argument("--samples", nargs="+", required=True)
    gundem2015_parser.add_argument("--min-sample-cnt", type=int, default=1)
    gundem2015_parser.add_argument("--no-flip-second-strand", action="store_false", dest="flip_second_strand")
    ####
    args = parser.parse_args()
    setup = build_setup(args=args)
    logger = get_standard_logger_from_args(args=args, program_name="RCK-UTILS-ADJ-x2rck")
    nas = []
    if args.o_extra_fields is None or len(args.o_extra_fields) == 0 or args.o_extra_fields == ",":
        extra = None
    elif args.o_extra_fields != "all":
        extra = args.o_extra_fields.split(",")
    else:
        extra = args.o_extra_fields
    if args.command == "lumpy":
        logger.info("Starting converting adjacencies from the Lumpy VCF format to that of RCK")
        logger.info("Reading Lumpy VCF records from {file}".format(file=args.lumpy_vcf_file))
        lumpy_vcf_records = get_vcf_records_from_source(source=args.lumpy_vcf_file)
        logger.info("Converting Lumpy VCF records to RCK adjacencies")
        nas = get_nas_from_lumpy_vcf_records(lumpy_vcf_records=lumpy_vcf_records, setup=setup)
    elif args.command == "longranger":
        logger.info("Starting converting adjacencies from the LongRanger VCf format to that of RCK")
        logger.info("Reading LongRanger VCF records from {file}".format(file=args.longranger_vcf_file))
        longranger_vcf_records = get_vcf_records_from_source(source=args.longranger_vcf_file)
        logger.info('Converting LongRanger VCF records to RCK adjacencies')
        nas = get_nas_from_longranger_vcf_records(longranger_vcf_records=longranger_vcf_records, setup=setup)
    elif args.command == "naibr":
        logger.info("Starting converting adjacencies from NAIBR records to that of RCK")
        logger.info("Reading and converting NAIBR records from {file} to RCK adjacencies".format(file=args.naibr_file))
        nas = get_nas_from_naibr_source(source=args.naibr_file, setup=setup)
    elif args.command == "manta":
        logger.info("Starting converting adjacencies from Manta records to that of RCK")
        logger.info("Reading Manta VCF records from {file}".format(file=args.manta_vcf_file))
        manta_vcf_records = get_vcf_records_from_source(source=args.manta_vcf_file)
        logger.info("Converting Manta VCF records to RCK adjacencies")
        nas = get_nas_from_manta_vcf_records(manta_vcf_records=manta_vcf_records, setup=setup)
    elif args.command == "sniffles":
        logger.info("Starting converting adjacencies from Sinffles records to that of RCK")
        logger.info("Reading Sniffles VCF records from {file}".format(file=args.sniffles_vcf_file))
        sniffles_vcf_records = get_vcf_records_from_source(source=args.sniffles_vcf_file)
        logger.info("Converting Sniffles VCF records to RCK adjacencies")
        nas = get_nas_from_sniffles_vcf_records(sniffles_vcf_records=sniffles_vcf_records, setup=setup)
    elif args.command == "grocsv":
        logger.info("Starting converting adjacencies from GROCSVS records to that of RCK")
        logger.info("Reading GROCSVS VCF records from {file}".format(file=args.grocsv_vcf_file))
        samples = args.samples.split(",") if args.samples is not None else args.samples
        grocsv_vcf_records = get_vcf_records_from_source(source=args.grocsv_vcf_file)
        logger.info("Converting GROCSVS VCF records to RCK adjacencies")
        nas = get_nas_from_grocsv_vcf_records(grocsv_vcf_records=grocsv_vcf_records, setup=setup, samples=samples, sample_all_any=args.samples_all_any,
                                              samples_only=args.samples_only)
    elif args.command == "delly":
        logger.info("Starting converting adjacencies from Delly records to that of RCK")
        if args.delly_force_stream:
            logger.info("Forced stream is enabled")
            logger.info("Streamlining reading and converting Delly VCF from {i_file} to RCK adjacencies into {o_file}"
                        "".format(i_file=args.delly_vcf_file, o_file=args.rck_nas_file))
            delly_vcf_to_nas_stream(source=args.delly_vcf_file, dest=args.rck_nas_file, setup=setup, extra=extra)
            sys.exit(0)
        else:
            logger.info("Reading Delly VCF records from {file}".format(file=args.delly_vcf_file))
            delly_vcf_records = get_vcf_records_from_source(source=args.delly_vcf_file)
            logger.info("Converting Delly VCF records to rCK adjacencies")
            nas = get_nas_from_delly_vcf_records(delly_vcf_records=delly_vcf_records, setup=setup)
    elif args.command == "pbsv":
        logger.info("Starting converting adjacencies from PBSV records to that of RCK")
        logger.info("Reading PBSV VCF records from {file}".format(file=args.pbsv_vcf_file))
        pbsv_vcf_records = get_vcf_records_from_source(source=args.pbsv_vcf_file)
        logger.info("Converting PBSV VCF records to RCK adjacencies")
        nas = get_nas_from_pbsv_vcf_records(pbsv_vcf_records=pbsv_vcf_records, setup=setup)
    elif args.command == "gundem2015":
        logger.info("Starting converting adjacencies from Gundem et al 2015 (BRASS2???) to that of RCK")
        logger.info("Reading Gundem 2015 et al (BRASS???) records from {file}".format(file=args.gundem2015_file))
        nas = get_nas_from_gundem2015_source(source=args.gundem2015_file, setup=setup, separator=args.i_separator, flip_second_strand=args.flip_second_strand)
        logger.info("Extracting adjacencies for sample {samples} with a minimum cnt of {min_cnt}".format(samples=",".join(args.samples), min_cnt=args.min_sample_cnt))
        nas = processed_gundem2015_adjacencies(adjacencies=nas, sample_names=args.samples, min_per_sample_cnt=args.min_sample_cnt)
    elif args.command == "remixt":
        logger.info("Starting converting adjacencies and their (haploid) copy numbers from ReMixT to that of RCK")
        logger.info("Reading and converting ReMixT resotds from {file} to RCK adjacencies".format(file=args.remixt_file))
        clone_ids = args.clone_ids.split(",")
        nas = get_nas_from_remixt_source(source=args.remixt_file, setup=setup, separator=args.i_separator, clone_ids=clone_ids, skip_absent=args.skip_absent,
                                         remixt_na_correction=args.remixt_correction)
    logger.info("A total of {cnt} adjacencies were obtained.".format(cnt=len(nas)))
    logger.debug("Output extra fields were identified as {o_extra}".format(o_extra=",".join(extra)))
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
    logger.debug("Include chromosomes : {include_chromosomes}".format(include_chromosomes=",".join(map(str, include_regions))))
    logger.debug("Exclude chromosomes : {exclude_chromosomes}".format(exclude_chromosomes=",".join(map(str, exclude_regions))))
    logger.info("Filtering adjacencies based on input/exclude chromosomes")
    nas = filter_adjacencies_by_chromosomal_regions(adjacencies=nas, include=include_regions, exclude=exclude_regions, include_both=args.include_both, exclude_both=args.exclude_both)
    nas = list(nas)
    logger.info("A total of {cnt} adjacencies were retained after filtering".format(cnt=len(nas)))
    logger.info("Writing RCK adjacencies to {file}".format(file=args.rck_adj_file))
    write_adjacencies_to_destination(destination=args.rck_adj_file, adjacencies=nas, extra=extra)
    logger.info("Success")


if __name__ == "__main__":
    main()
