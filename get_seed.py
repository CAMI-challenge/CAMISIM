#!/usr/bin/env python

import random
import sys

    



if __name__ == "__main__":

    seed = int(sys.argv[1])
    count_samples = int(sys.argv[2])
    file_genome_locations = sys.argv[3]

    genome_id_list = []

    location_file = open(file_genome_locations, "r")

    for line in location_file:
        genome_id = line.split('\t')[0]
        genome_id_list.append(genome_id)

    location_file.close()

    text = "used_initial_seed" + '\t' + str(seed) + '\n'
    text = text + "genome_id" + '\t' + "sample_id" + '\t' + "seed" + '\n'

    random.seed(seed)
    
    f = open("seed.txt", "w")
    

    for i in range(count_samples):
        for genome in genome_id_list:
            sample_seed = random.randint(0, sys.maxsize)
            text = text + genome + '\t' + str(i) + '\t' + str(sample_seed) + '\n'
        
    f.write(text)
    f.close()

