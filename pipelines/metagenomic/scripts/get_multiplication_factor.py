#!/usr/bin/env python
from Bio import SeqIO
import os
import sys


def main(total_size, genome_locations, file_path_distribution, project_dir):
    abundances = {}
    # get the abundances from the distribution files
    with open(file_path_distribution, 'r') as ab:
        for line in ab:
            genome_id, abundance = line.strip().split('\t')
            abundances[genome_id] = float(abundance)
    total = sum(abundances.values())
    # normalise to 1
    abundances = { x : abundances[x]/total for x in abundances }
    # match abundances with genomes and normalise by genome size
    total_relative_size = 0
    with open(genome_locations, 'r') as loc:
        for line in loc:
            genome_id, location = line.strip().split('\t')
            relative_size = 0
            if os.path.isabs(location):
                for record in SeqIO.parse(location,"fasta"):
                    relative_size += abundances[genome_id] * len(record.seq)
            else:
                path = os.path.join(project_dir,location)
                for record in SeqIO.parse(path,"fasta"):
                    relative_size += abundances[genome_id] * len(record.seq)
            total_relative_size += relative_size
    print(total_size / float(total_relative_size))


if __name__ == "__main__":
    main (total_size = float(sys.argv[1]),
          genome_locations = sys.argv[2],
          file_path_distribution = sys.argv[3],
          project_dir = sys.argv[4])
