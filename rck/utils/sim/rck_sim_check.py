import argparse
import os
import sys
from typing import Set, List, Dict

from pysam.libcvcf import defaultdict

from rck.core.structures import Genome, AdjacencyType, Adjacency


def genome_name_from_file_path(path: str) -> str:
    return str(os.path.basename(path).split(".")[0])


def sim_check_as_strs(ref_genome: Genome, rearr_genomes: Dict[str, Genome], ref_genome_name: str = "R") -> List[str]:
    result: List[str] = []
    ref_genome.propagate_hap_labeles_to_positions()
    ref_segments = {s for s in ref_genome.iter_segments(allow_ins=False, ins_prefix="INS")}
    ref_adjacencies = {adj for adj in ref_genome.iter_adjacencies()}
    ref_telomeres = set(ref_genome.iter_telomeres())
    diploid_novel_adjacencies_per_position = defaultdict(set)
    haploid_novel_adjacencies_per_position = defaultdict(set)
    for rearr_genome_name, rearr_genome in rearr_genomes.items():
        rearr_genome.propagate_hap_labeles_to_positions()
        rearr_segments = {s for s in rearr_genome.iter_segments(allow_ins=False, ins_prefix="INS")}
        if not rearr_segments.issubset(ref_segments):
            result.append(f"BAD: some non-ins segments in rearranged genome {rearr_genome_name} are not present in the refined reference: {rearr_segments - ref_segments}")
        else:
            result.append(f"GOOD: all segment in rearranged genome {rearr_genome_name} are present in the reference")
        rearr_adjacencies: List[Adjacency] = [adj for adj in rearr_genome.iter_adjacencies()]
        rearr_adjacencies_ref: Set[Adjacency] = {adj for adj in rearr_adjacencies if adj.adjacency_type == AdjacencyType.REFERENCE}
        if not rearr_adjacencies_ref.issubset(ref_adjacencies):
            result.append(f"BAD: some adjacencies annotated as 'reference' in the rearranged genome {rearr_genome_name} are not"
                          f" present in a refined reference: {rearr_adjacencies_ref - ref_adjacencies}")
        else:
            result.append(f"GOOD: all adjacencies annotated as 'reference' in the rearranged genomes {rearr_genome_name} are present in the refined reference")
        rearr_adjacencies_novel: Set[Adjacency] = {adj for adj in rearr_adjacencies if adj.adjacency_type == AdjacencyType.NOVEL}
        for adj in rearr_adjacencies_novel:
            diploid_novel_adjacencies_per_position[adj.position1].add(adj.position2)
            diploid_novel_adjacencies_per_position[adj.position2].add(adj.position1)
            p1_hap = adj.position1.get_non_hap_copy()
            p2_hap = adj.position2.get_non_hap_copy()
            haploid_novel_adjacencies_per_position[p1_hap].add(adj.position2)
            haploid_novel_adjacencies_per_position[p2_hap].add(adj.position1)
        rearr_genome_telomeres = set(rearr_genome.iter_telomeres())
        if not rearr_genome_telomeres.issubset(ref_telomeres):
            result.append(f"BAD: some telomeres in rearranged genome {rearr_genome_name} are not telomeres in the refined reference: {rearr_genome_telomeres - ref_telomeres}")
        else:
            result.append(f"GOOD: all telomeres in the rearranged genome {rearr_genome_name} are also telomeres in the refined reference")
    diploid_is_violations = set()
    haploid_is_violations = set()
    for key, values in diploid_novel_adjacencies_per_position.items():
        if len(values) > 1:
            diploid_is_violations.add(key)
    for key, values in haploid_novel_adjacencies_per_position.items():
        if len(values) > 1:
            haploid_is_violations.add(key)
    if len(diploid_is_violations) > 1:
        result.append(f"BADish: some haplotype-specific positions are involved in more than 1 novel adjacencies: {','.join(map(str, diploid_is_violations))}")
    else:
        result.append(f"GOOD: no diploid IS violations")
    if len(haploid_is_violations) > 1:
        result.append(f"BADish: some haploid positions are involved in more than 1 novel adjacencies: {','.join(map(str, haploid_is_violations))}")
    else:
        result.append("GOOD: no haploid IS violations")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ref", type=argparse.FileType("rt"))
    parser.add_argument("r_genomes", type=argparse.FileType("rt"), nargs="+")
    parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()
    ref_genome = Genome.from_strs(args.ref)
    ref_genome_name = genome_name_from_file_path(args.ref.name)
    rearr_genomes = {genome_name_from_file_path(source.name): Genome.from_strs(source) for source in args.r_genomes}
    result = sim_check_as_strs(ref_genome, rearr_genomes, ref_genome_name)
    print("\n".join(result))


if __name__ == "__main__":
    main()
