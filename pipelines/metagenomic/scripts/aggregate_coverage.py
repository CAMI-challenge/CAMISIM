#!/usr/bin/env python3
import sys
import os
from typing import Dict, Set, Tuple, TextIO


def parse_fasta_sequence_ids(fasta_path: str) -> Set[str]:
    """Extracts sequence identifier tokens from FASTA headers of a genome file."""
    seq_ids: Set[str] = set()
    if not os.path.exists(fasta_path):
        # return seq_ids
        raise FileNotFoundError(f"Genome FASTA listed in genome_locations.tsv does not exist: {fasta_path}")

    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                # Extract first word after '>' as the sequence ID
                seq_id = line.split()[0][1:]
                seq_ids.add(seq_id)
    return seq_ids


def load_genome_mappings(locations_path: str) -> Tuple[Dict[str, str], Set[str]]:
    """
    Parses genome locations mapping and maps every sequence identifier
    to its respective FASTA file basename.
    """
    seq_to_fasta: Dict[str, str] = {}
    fasta_basenames: Set[str] = set()

    with open(locations_path, 'r', encoding='utf-8') as f:
        for line in f:
            line_stripped = line.strip()
            if not line_stripped or line_stripped.startswith('#'):
                continue

            parts = line_stripped.split('\t')
            if len(parts) < 2 or parts[0] == 'genome':
                continue

            gpath = parts[1]
            fbname = os.path.basename(gpath)
            fasta_basenames.add(fbname)

            for seq_id in parse_fasta_sequence_ids(gpath):
                seq_to_fasta[seq_id] = fbname

    return seq_to_fasta, fasta_basenames


def aggregate_samtools_coverage(
    stdin_stream: TextIO,
    seq_to_fasta: Dict[str, str],
    fasta_basenames: Set[str]) -> Tuple[Dict[str, float], Dict[str, int]]:
    """Aggregates length-weighted mean depth per genome from samtools coverage records."""
    genome_sum: Dict[str, float] = {fbname: 0.0 for fbname in fasta_basenames}
    genome_len: Dict[str, int] = {fbname: 0 for fbname in fasta_basenames}

    for line in stdin_stream:
        line_stripped = line.strip()
        if line_stripped.startswith('#') or not line_stripped:
            continue

        parts = line_stripped.split('\t')
        if len(parts) < 7:
            continue

        rname = parts[0]
        startpos = int(parts[1])
        endpos = int(parts[2])
        region_len = endpos - startpos + 1
        meandepth = float(parts[6])

        fbname = seq_to_fasta.get(rname)
        if fbname:
            genome_sum[fbname] += meandepth * region_len
            genome_len[fbname] += region_len
        else:
            print(f"WARNING: reference sequence '{rname}' not found in genome_locations-derived FASTA mapping", file=sys.stderr)

    return genome_sum, genome_len


def write_aggregated_coverage(
    output_path: str,
    fasta_basenames: Set[str],
    genome_sum: Dict[str, float],
    genome_len: Dict[str, int]) -> None:
    """Writes the length-weighted coverage per genome to a TSV file."""
    with open(output_path, 'w', encoding='utf-8') as outf:
        for fbname in sorted(fasta_basenames):
            glen = genome_len[fbname]
            if glen <= 0 or genome_sum[fbname] <= 0.0:
                continue
            cov = genome_sum[fbname] / glen
            outf.write(f"{fbname}\t{cov}\n")


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: aggregate_coverage.py <genome_locations.tsv> <output_coverage.tsv>", file=sys.stderr)
        sys.exit(1)

    genome_locations_path = sys.argv[1]
    output_cov_path = sys.argv[2]

    seq_to_fasta, fasta_basenames = load_genome_mappings(genome_locations_path)
    genome_sum, genome_len = aggregate_samtools_coverage(sys.stdin, seq_to_fasta, fasta_basenames)
    write_aggregated_coverage(output_cov_path, fasta_basenames, genome_sum, genome_len)


if __name__ == '__main__':
    main()
