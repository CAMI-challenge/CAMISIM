#!/usr/bin/env python

"""
This script can create (pooled) gold standards for assembly and binning of
a) a subset of samples of a single CAMISIM run
b) a subset of samples from two CAMISIM runs with different sequencing technologies (requires same seed)
"""

import sys
import os
import subprocess
import shutil
import argparse

def parse_options():
    """
    parse the command line options
    """
    parser = argparse.ArgumentParser()

    helptext = "Root path of input runs to be considered, can be one or multiple CAMISIM runs (if more than one, they are required to have the same random seed/genome mapping)\nSample folder names are expected to follow this schema: yyyy.mm.dd_hh.mm.ss_sample_"
    parser.add_argument("-i", "--input-runs", type=str, help=helptext, nargs='+') 

    helptext = "Samples to be considered for pooled gold standards. If none are provided, pooled gold standard is created over all samples"
    parser.add_argument("-s", "--samples", type=int, help=helptext, nargs='*')

    helptext = "Output directory for all gold standards and files"
    parser.add_argument("-o", "--output-directory", type=str, help=helptext)

    helptext = "Number of threads to be used, default 1"
    parser.add_argument("-t", "--threads", type=int, default=1,help=helptext)

    helptext = "Path to the bamToGold perl script"
    parser.add_argument("-b", "--bamToGold", type=str, help=helptext)

    if not len(sys.argv) > 1:
        parser.print_help()
        return None
    args = parser.parse_args()

    return args


def get_samples(root_paths, samples):
    """
    Given the root paths  of the CAMISIM runs and the subset of samples, returns a dict from sample number to folders
    Assumes the sample folders to be in the format YYYY.MM.DD_HH.MM.SS_sample_#
    """
    used_samples = {}
    for path in root_paths:
        if not os.path.exists(path):
            raise IOError("No such file or directory: %s" % path)
        files = os.listdir(path)
        for f in files:
            try:
                date, time, sample, nr = f.split("_")
            except ValueError:
                continue
            if samples is None or nr in samples:
                if nr in used_samples:
                    used_samples[nr].append(os.path.join(path,f))
                else:
                    used_samples[nr] = [os.path.join(path,f)]
    return used_samples

def read_metadata(root_paths):
    """
    Reads the metadata files of the runs to create binning gold standards later on
    """
    metadata = {}
    for path in root_paths:
        if not os.path.exists(path):
            raise IOError("No such file or directory: %s" % path)
        metadata_path = os.path.join(path, "metadata.tsv")
        if not os.path.exists(metadata_path):
            raise IOError("Metadata file not found in %s" % path)
        with open(metadata_path,'r') as md:
            for line in md:
                if line.startswith("genome"):
                    continue
                genome, otu, ncbi, novelty = line.strip().split('\t')
                if genome in metadata:
                    set_otu, set_ncbi, set_novelty = metadata[genome][:3]
                    if otu != set_otu or ncbi != set_ncbi or novelty != set_novelty:
                        raise IOError("Metadata between runs differs, different environments and/or seeds have been used")
                else:
                    metadata[genome] = [otu, ncbi, novelty]
        genome_to_id_path = os.path.join(path, "genome_to_id.tsv")
        with open(genome_to_id_path, 'r') as gid:
            for line in gid:
                genome, path = line.strip().split('\t')
                if genome in metadata:
                    set_path = metadata[genome][-1]
                    if len(metadata[genome]) > 3 and set_path.rsplit("/",1)[-1] != path.rsplit("/",1)[-1]: 
                        # genome path has been set and differs (only fasta file)
                        raise IOError("genome_to_id between runs differs, different environments and/or seeds have been used: %s | %s" % (path, set_path))
                    if len(metadata[genome]) == 3:
                        metadata[genome].append(path)
                else: # this should not happen
                    raise IOError("Genome found in genome_to_id without metadata, check your CAMISIM run")
    return metadata

def bamToGold(bamtogold, merged, out, metadata, threads):
    """
    Calls the bamToGold script for all of the merged bam files, creating the gold standard
    """
    out_name = os.path.join(out, "anonymous_gsa.fasta")
    all_files = os.listdir(merged)
    bams = []
    for f in all_files:
        if f.endswith(".bam"):
            bams.append(f)
    for bam in bams:
        genome = bam.rstrip(".bam")
        otu, ncbi, novelty, path = metadata[genome]
        cmd = "{bamToGold} -r {path} -b {bam} -l 1 -c 1 >> {gsa}".format(
            bamToGold = bamtogold,
            path = path,
            bam = os.path.join(out,"bam",bam),
            gsa = out_name
        )
        subprocess.call([cmd],shell=True)

def fix_headers(genome, bams, out):
    """
    Sometimes short sequences are not present in one or the other bam file since no reads were created for them, the header of the merged file needs to have the union of all sequences in all bam files
    """
    sequences = set()
    header = ""
    header_name = os.path.join(out, genome + "_header.sam")
    for bam in bams:
        cmd = "samtools view -H {bam} > {hname}".format(
            bam = bam,
            hname = header_name
        )
        subprocess.call([cmd],shell=True)
        with open(header_name, 'r') as header_file:
            for line in header_file:
                if line.startswith("@SQ"): # sequences 
                    sq, sn, ln = line.strip().split('\t')
                    sequence_name = sn.split(":",1)[1]
                    if sequence_name not in sequences:
                        sequences.add(sequence_name)
                        header += line
                elif line.startswith("@HD") and header == "": #primary header, use arbitrary one
                    header += line 
    with open(header_name,'w+') as header_file:
        header_file.write(header)
    return header_name

def merge_bam_files(bams_per_genome, out, threads):
    """
    Merges (+sort +index)  all given bam files per genome (exact paths, single sample/multiple runs or multiple samples)
    """
    out_path = os.path.join(out,"bam")
    os.mkdir(out_path)
    for genome in bams_per_genome:
        list_of_bam = " ".join(bams_per_genome[genome]) # can be used as input to samtools immediately
        header = fix_headers(genome, bams_per_genome[genome], out_path)
        if header is not None:
            for bam in bams_per_genome[genome]: # add new header to all bam files
                cmd = "samtools reheader {header} {bam} >> {out}/out.bam; mv {out}/out.bam {bam}".format(
                    header = header,
                    out = out_path,
                    bam = bam
                )
                subprocess.call([cmd],shell=True)
        cmd = "samtools merge -@ {threads} - {bam_files} | samtools sort -@ {threads} - {path}/{genome}; samtools index {path}/{genome}.bam".format(
            threads = threads,
            bam_files = list_of_bam,
            path = out_path,
            genome = genome
        )
        subprocess.call([cmd],shell=True) # this runs a single command at a time (but that one multi threaded)
    return out_path

def create_gold_standards(bamtogold, used_samples, metadata, out, threads):
    """
    Creation of the gold standards per sample. Uses the helper script bamToGold and merges all bam files of the same genome per sample across runs
    """
    for sample in used_samples:
        runs = used_samples[sample]
        bam_per_genome = {}
        for run in runs:
            bam_dir = os.path.join(run,"bam")
            all_files = os.listdir(bam_dir)
            bam_files = []
            for f in all_files:
                if f.endswith(".bam"):
                    bam_files.append(f)
            for bam_file in bam_files:
                genome = bam_file.rstrip(".bam")
                if genome in bam_per_genome:
                    bam_per_genome[genome].append(os.path.join(run,"bam",bam_file))
                else:
                    bam_per_genome[genome] = [os.path.join(run,"bam",bam_file)]
        sample_path = os.path.join(out,"sample_%s" % sample) # creating a folder for every sample
        os.mkdir(sample_path)
        merged = merge_bam_files(bam_per_genome, sample_path, threads)
        bamToGold(bamtogold, merged, sample_path, metadata, threads)

def create_pooled_gold_standard(bamtogold, used_samples, metadata, out, threads):
    bam_per_genome = {}
    for sample in used_samples:
        runs = used_samples[sample]
        for run in runs:
            bam_dir = os.path.join(run,"bam")
            all_files = os.listdir(bam_dir)
            bam_files = []
            for f in all_files:
                if f.endswith(".bam"):
                    bam_files.append(f)
            for bam_file in bam_files:
                genome = bam_file.rstrip(".bam")
                if genome in bam_per_genome:
                    bam_per_genome[genome].append(os.path.join(run,"bam",bam_file))
                else:
                    bam_per_genome[genome] = [os.path.join(run,"bam",bam_file)]
    bam_pooled = os.path.join(out, "pooled")
    os.mkdir(bam_pooled)
    merged = merge_bam_files(bam_per_genome, bam_pooled, threads)
    bamToGold(bamtogold, merged, bam_pooled, metadata, threads)

if __name__ == "__main__":
    args = parse_options()
    if not args is None:
        root_paths = args.input_runs # list of input paths
        samples = args.samples
        out = args.output_directory
        if not os.path.exists(out):
            os.mkdir(out)
        threads = args.threads
        bamtogold = args.bamToGold
        used_samples = get_samples(root_paths, samples)
        metadata = read_metadata(root_paths)
        if len(root_paths) > 1: # do create individual gold standards per sample
            create_gold_standards(bamtogold, used_samples, metadata, out, threads)
        create_pooled_gold_standard(bamtogold, used_samples, metadata, out, threads) # in any case, create pooled gold standard
