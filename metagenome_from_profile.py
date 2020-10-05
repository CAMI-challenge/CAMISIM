#!/usr/bin/env python

import sys
from scripts.Validator.validator import Validator
from scripts.configfilehandler import ConfigFileHandler
from scripts.loggingwrapper import LoggingWrapper as logger
from configparser import ConfigParser
import scripts.get_genomes as GG
import shutil
import os
import argparse

def parse_options():
    parser = argparse.ArgumentParser()
    
    helptext="16S profile to create metagenome from. Can either be CAMI-format or biom-format."
    parser.add_argument("-p","--profile", default=None, type=str, help=helptext)

    helptext="Number of samples to be generated. If nothing is given, this defaults to 1 (CAMI format) or the number of samples present in the biom file. If a specific number is given, the samples are simulated using the first sample of the biom file"
    parser.add_argument("-s","--samples", default=None, type=int, help=helptext)

    #helptext="Whether the related genomes are supposed to be downloaded."
    # download the mapped full genomes?
    #parser.add_argument("-no-dl", "--dont-download-genomes", action='store_false', default=True, help=helptext)

    helptext="Output directory, default: out/"
    # out path
    parser.add_argument("-o", default="out/", type=str, help=helptext,metavar="OUT PATH")
    
    helptext="Path where temporary files are stored (gets deleted after pipeline is finished)"
    # temporary path
    parser.add_argument("-tmp", default=None, type=str, help=helptext)

    default="tools/assembly_summary_complete_genomes.txt"
    # file pointing to reference genomes
    parser.add_argument("-ref","--reference-genomes", default=default, help="File pointing to reference genomes of the format: NCBI id\tScientific name\tNCBI ftp address of full genome. Default: %s" % default)

    default=None
    # additional files containing genomes which should be considered for mapping
    parser.add_argument("-ar", "--additional-references", default=default, help="File containing additional reference genomes, mapped to OTUs from the input profile")

    default = "defaults/default_config.ini"
    # optional config file (out_path will get overwritten if it is set in config file)
    parser.add_argument("-c","--config",default=default,help="Path to config file. Careful when setting \"metadata\", \"id_to_genome_file\", \"distribution_file_paths\"(they will be set by the pipeline) and the out path differently from the command line out path, default: %s" % default,metavar="CONFIG FILE")

    default = "tools/ncbi-taxonomy_20170222.tar.gz"
    parser.add_argument("--ncbi",default=default,help="Path to the NCBI taxdump for finding corresponding reference genomes, default = %s" % default)

    parser.add_argument("-nr", "--no-replace", action='store_false',default=True, help="Use sampling without replacing, so genomes are used for exactly one OTU only (decreases accuracy)")
    
    helptext = "If no genomes are found for certain OTUs, fill up with previously unused genomes"
    parser.add_argument("-f", "--fill-up", action='store_true',default=False,help=helptext)

    helptext = "Only perform community design, do not simulate"
    parser.add_argument("-d", "--community-only", action='store_true', default=False, help=helptext)
    
    helptext="Seed for the random generator"
    parser.add_argument("--seed",type=int,default=None,help=helptext)

    parser.add_argument("--debug",action='store_true',default=False,help="get more debug information")
    if not len(sys.argv) > 1:
        parser.print_help()
        return None
    args = parser.parse_args()
    
    return args

def create_config(args,cfg):
    config = ConfigParser()
    config.read(cfg)

    config.set('Main', 'output_directory', os.path.join(args.o,''))
    if args.tmp is not None:
        config.set('Main', 'temp_directory', os.path.join(args.tmp,''))
    else:
        config.set('Main', 'temp_directory', "/tmp")

    if args.seed is not None:
        config.set('Main', "seed", args.seed)
    name = os.path.join(args.o,"config.ini")
    with open(name,'w+') as cfg_path:
        config.write(cfg_path)
    return name

if __name__ == "__main__":
    args = parse_options()
    if not args is None:
        log = logger(verbose = args.debug)
        if args.debug:
            log.info("Using commands:")
            for arg in vars(args):
                log.info("-%s: %s" % (arg, getattr(args,arg)))
        if not os.path.exists(args.o):
            os.mkdir(args.o)
        config = GG.generate_input(args) # total number of genomes and path to updated config
        c = create_config(args,config)
        if (not args.community_only):
            if args.debug:
                os.system("./metagenomesimulation.py %s --debug" % c)
            else:
                os.system("./metagenomesimulation.py %s" % c)
