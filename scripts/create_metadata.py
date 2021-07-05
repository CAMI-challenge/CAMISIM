#!/usr/bin/env python

"""
This script creates a metadata file in json format
"""

import json
import argparse
import os
import sys
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0

def parse_options():
    """
    parse command line options
    """
    parser = argparse.ArgumentParser()

    helptext="Root path of input run for which metadata should be created, should contain metadata.tsv and genome_to_id.tsv"
    parser.add_argument("-i", "--input-run", type=str, help=helptext)

    helptext="output file to write metadata to"
    parser.add_argument("-o", "--output", type=str, help=helptext)
    
    helptext="Name of the data set"
    parser.add_argument("-n", "--name", type=str, help=helptext)

    if not len(sys.argv) > 1:
        parser.print_help()
        return None
    args = parser.parse_args()

    return args

def read_metadata(path):
    """
    Reads the metadata files for creating some of the mappings
    """
    metadata = {}
    if not os.path.exists(path):
        raise IOError("No such file or directory: %s" % path)
    metadata_path = os.path.join(path, "metadata.tsv")
    genome_path = os.path.join(path, "genome_to_id.tsv")
    with open(metadata_path,'r') as md: 
        for line in md: 
            if line.startswith("genome"):
                continue
            genome, otu, ncbi, novelty = line.strip().split('\t')
            metadata[genome] = [otu, ncbi, novelty]
    with open(genome_path, 'r') as gid:
        for line in gid:
            genome, path = line.strip().split('\t')
            metadata[genome].append(path)
    return metadata
   
def read_config(path):
    """
    Reads the config file used for the creation of the data set
    """
    config_path = os.path.join(path, "config.ini")
    config = ConfigParser()
    config.read(config_path)

    return config

def get_sequencing_technology(read_simulator):
    """
    Given the read simulator type from the config file, returns the underlying sequencing technology (this has to be edited if new sequencing technologies/reads simulator are added)
    """
    dict_of_simulators = {
        "art" : "Illumina HiSeq2500",
        "wgsim" : "Illumina",
        "pbsim" : "Pacific Biosciences",
        "nanosim" : "Oxford Nanopore",
        "nanosim" : "Oxford Nanopore, R9 chemistry"
    }
    return dict_of_simulators[read_simulator]

def get_bam_list(path,sample):
    bam_path = os.path.join(path,sample,"bam")
    bam_files = os.listdir(bam_path)
    bam_list = []
    for f in bam_files:
        if f.endswith(".bam"):
            bam_list.append(os.path.join(bam_path,f))
    return bam_list

def get_sample_dicts(path):
    files = os.listdir(path)
    sample_dicts = []
    for f in files:
        try:
            date, time, sample, nr = f.split("_")
        except ValueError:
            continue 
        sample_dict = {}
        sample_dict["Number"] = nr
        sample_dict["Name"] = f
        sample_dict["Reads"] = os.path.join(f,"reads","anonymous_reads.fq.gz")
        sample_dict["Gold_Standard_Taxonomic_Profile"] = "taxonomic_profile_%s.txt" % nr
        sample_dict["Gold_Standard_Assembly"] = os.path.join(f,"contigs","anonymous_gsa.fasta.gz")
        sample_dict["Gold_Standard_Read_Binning"] = os.path.join(f,"reads","reads_mapping.tsv.gz")
        sample_dict["Gold_Standard_Contig_Binning"] = os.path.join(f,"contigs", "gsa_mapping.tsv.gz")
        bam_list = get_bam_list(path,f)
        sample_dict["Read_to_genome_mappings"] = bam_list
        sample_dicts.append(sample_dict)
    return sample_dicts

def list_genomes(metadata):
    genome_dicts = []
    for genome in metadata: #otu, ncbi, novelty, path
        genome_dict = {}
        otu, ncbi, novelty, path = metadata[genome]
        genome_dict["OTU_Name"] = genome
        genome_dict["OTU_NCBI_ID"] = otu
        genome_dict["Genome_NCBI_ID"] = ncbi
        genome_dict["Novelty_category"] = novelty
        genome_dict["Path"] = path
        genome_dicts.append(genome_dict)
    return genome_dicts

def create_json(path, metadata, config, name):
    json = {}
    json["Name"] = name
    samples = config.get('CommunityDesign', 'number_of_samples')
    json["Nr_Samples"] = samples
    size = config.get('ReadSimulator', 'size')
    json["Total_size_bp"] = float(size) * int(samples)
    simulator = config.get('ReadSimulator', 'type')
    json["Sequencing_technology"] = get_sequencing_technology(simulator)
    if simulator != "art" and simulator != "wgsim": # these use insert sizes
        json["Average_read_length_in_bp"] = 150
        json["Read_length_standard_deviation_in_bp"] = 0
        json["Average_insert_size_in_bp"] = config.get('ReadSimulator', 'fragments_size_mean')
        json["Insert_size_standard_deviation_in_bp"] = config.get('ReadSimulator', 'fragment_size_standard_deviation')
    else: # these use read sizes
        json["Average_read_length_in_bp"] = config.get('ReadSimulator', 'fragments_size_mean')
        json["Read_length_standard_deviation_in_bp"] = config.get('ReadSimulator', 'fragment_size_standard_deviation')
        json["Average_insert_size_in_bp"] = config.get('ReadSimulator', 'fragments_size_mean')
        json["Insert_size_standard_deviation_in_bp"] = config.get('ReadSimulator', 'fragment_size_standard_deviation')
    if simulator == "art":
       json["Paired-end"] = True
    else:
       json["Paired-end"] = False
    sample_dicts = get_sample_dicts(path)
    json["Samples"] = sample_dicts
    file_dict = {}
    file_dict["Pooled_Gold_Standard_Assembly"] = os.path.join(path, "anonymous_gsa_pooled.fasta.gz")
    file_dict["Pooled_Gold_Standard_Binning"] = os.path.join(path, "gsa_pooled_mapping.tsv.gz")
    file_dict["Genomes"] = list_genomes(metadata)
    file_dict["CAMISIM_config"] = os.path.join(path, "config.ini")
    file_dict["CAMISIM_metadata"] = os.path.join(path, "metadata.tsv")
    file_dict["CAMISIM_genome_mapping"] = os.path.join(path, "genome_to_id.tsv")
    json["Other_files"] = file_dict
    return json

def write_json(in_path, out_path, metadata, config, name):
    json_path = os.path.join(out_path, "metadata.json")
    json_dict = create_json(in_path, metadata, config, name)
    json_string = json.dumps(json_dict)
    with open(json_path, 'w') as json_file:
        json_file.write(json_string)

if __name__ == "__main__":
    args = parse_options()
    if not args is None:
        inpath = args.input_run
        outpath = args.output
        name = args.name
        metadata = read_metadata(inpath)
        config = read_config(inpath)
        write_json(inpath, outpath, metadata, config, name)

