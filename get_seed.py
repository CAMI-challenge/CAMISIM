#!/usr/bin/env python

import random
import sys
import argparse
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
		"-seed",
		help="the initial seed",
		action='store',
		default="")
    parser.add_argument(
		"-count_samples",
		help="the sample count",
		action='store',
		default="")  
    parser.add_argument(
		"-file_genome_locations",
		help="the file containing the genome locations",
		action='store',
		default="")      
    parser.add_argument(
		"-anonym_seed",
		help="whether seeds for anonymization should be generated",
		action="store_true",
		default=False)
    options = parser.parse_args()    

    seed = int(options.seed)
    count_samples = int(options.count_samples)
    file_genome_locations = options.file_genome_locations
    anonym_seed = options.anonym_seed

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

    if(anonym_seed):

        # file with seeds for read anonymization
        text = "used_initial_seed" + '\t' + str(seed) + '\n'
        text = text + "sample_id" + '\t' + "seed" + '\n'

        f = open("seed_read_anonymisation.txt", "w")

        for i in range(count_samples):
            sample_seed = random.randint(0, sys.maxsize)
            text = text + str(i) + '\t' + str(sample_seed) + '\n'
        
        f.write(text)
        f.close()

        # file with seeds for gsa anonymization
        text = "used_initial_seed" + '\t' + str(seed) + '\n'
        text = text + "sample_id" + '\t' + "seed" + '\n'

        f = open("seed_gsa_anonymisation.txt", "w")

        for i in range(count_samples):
            sample_seed = random.randint(0, sys.maxsize)
            text = text + str(i) + '\t' + str(sample_seed) + '\n'
        
        f.write(text)
        f.close()

