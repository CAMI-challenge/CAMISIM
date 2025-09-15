#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

def get_new_name(name, set_of_sequence_names):
    index = 0
    new_name = name
    while new_name in set_of_sequence_names:
        new_name = name + "_{}".format(index)
        index += 1
    return new_name

def add_sequences_to_map(stream_input, stream_map, genome_id, sequence_min_length, set_of_sequence_names, file_format="fasta"):
    for seq_record in SeqIO.parse(stream_input, file_format):
        # remove description, else art illumina messes up sam format
        seq_record.description = ''
        if len(seq_record.seq) < sequence_min_length:
            #self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
            #	os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
            continue
        if seq_record.id in set_of_sequence_names:
            new_id = get_new_name(seq_record.id, set_of_sequence_names)
            stream_map.write("{}\t{}\t{}\n".format(genome_id, seq_record.id, new_id))
            seq_record.id = new_id
        set_of_sequence_names.add(seq_record.id)

def cleanup_and_filter_sequences(
    stream_input, stream_output, stream_map,
    genome_id, sequence_min_length, set_of_sequence_names, file_format="fasta"):

    total_base_pairs = 0
    for seq_record in SeqIO.parse(stream_input, file_format):
        # remove description, else art illumina messes up sam format
        seq_record.description = ''
        if len(seq_record.seq) < sequence_min_length:
            #self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
            #	os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
            continue
        if seq_record.id in set_of_sequence_names:
            new_id = get_new_name(seq_record.id, set_of_sequence_names)
            stream_map.write("{}\t{}\t{}\n".format(genome_id, seq_record.id, new_id))
            seq_record.id = new_id
        set_of_sequence_names.add(seq_record.id)
        # file_handler.write(">{}\n".format(sequence_id))
        # file_handler.writelines("{}\n".format(seq_record.seq))
        stream_output.write(seq_record.format(file_format))
        total_base_pairs += len(seq_record.seq)
    return total_base_pairs

def move_genome_file(
    file_path_input, file_path_output,
    stream_map, genome_id, sequence_min_length=1, set_of_sequence_names=None, file_format="fasta"):

    #assert self.validate_file(file_path_input)
    assert file_format == "fasta", "'{}' is not supported, yet.".format(file_format)
    #assert not self.validate_file(file_path_output, silent=True), "Overwriting files prohibited: '{}'".format(
    #	file_path_output)
    if set_of_sequence_names is None:
        set_of_sequence_names = []
        
    #if (self.validate_file(file_path_output, silent=True)):
        #self._logger.warning("File %s existing, skipping" % file_path_output)
    #    with open(file_path_input, 'r') as stream_input:
    #        add_sequences_to_map(stream_input, stream_map, genome_id, sequence_min_length, set_of_sequence_names)
    #    return
    with open(file_path_input, 'r') as stream_input, open(file_path_output, 'w') as stream_output:
        total_base_pairs = cleanup_and_filter_sequences(
            stream_input, stream_output, stream_map, genome_id, sequence_min_length, set_of_sequence_names, file_format)
    if total_base_pairs == 0:
        msg = "No valid sequences in '{}'".format(stream_input.name)
        #self._logger.error(msg)
        raise Exception(msg)
    
def parse_tsv_to_dict(file):
    result_dict = {}
    with open(file, 'r') as f:
        for line in f:
            key, value = line.strip().split('\t')
            result_dict[key] = value
    return result_dict

def write_genome_id_to_path_map(genome_id_to_path_map, file_path_output, outdir):
    """
    Write mapping of genome id to genome file path to a file.

    @param file_path_output: File path
    @type file_path_output: file | FileIO | StringIO
    @param genome_id_to_path_map:
    @type genome_id_to_path_map: dict[str|unicode, str|unicode]
    """
    with open(file_path_output, 'w') as stream_out:
        stream_genome_id_to_path_map(stream_out, genome_id_to_path_map, outdir)

def stream_genome_id_to_path_map(stream_out, genome_id_to_path_map, outdir):
    """
    Write mapping of genome id to genome file path to a stream.

    @param stream_out: Stream like object
    @type stream_out: file | FileIO | StringIO
    @param genome_id_to_path_map:
    @type genome_id_to_path_map: dict[str|unicode, str|unicode]
    """
    #assert self.is_stream(stream_out)
    for genome_id, file_path in genome_id_to_path_map.items():
        out_filepath = outdir + file_path.rsplit("/")[-1]
        stream_out.write("{}\t{}\n".format(genome_id, out_filepath))

if __name__ == "__main__":
    
    genome_id_to_file_path_file = sys.argv[1]
    outdir = sys.argv[2]
    internal_genome_id_to_file_path_file = sys.argv[3]

    filename_seq_map = "sequence_id_map.txt"
    sequence_min_length = 1
    set_of_sequence_names = set()
    directory_output = "./out_genomes/"
    genome_id_to_path_map = dict()

    # Create the output directory if it does not exist
    os.makedirs(directory_output, exist_ok=True)
    
    genome_id_to_file_path = parse_tsv_to_dict(genome_id_to_file_path_file)

    file_path_sequence_map = filename_seq_map
    with open(file_path_sequence_map, 'w') as stream_map:

        for genome_id, genome_file_path in genome_id_to_file_path.items():
            file_name = os.path.basename(genome_file_path)

            new_genome_file_path = os.path.join(directory_output, file_name)

            move_genome_file(
                file_name, new_genome_file_path, stream_map, genome_id, sequence_min_length, set_of_sequence_names)
            genome_id_to_path_map[genome_id] = new_genome_file_path
        
            
    write_genome_id_to_path_map(genome_id_to_path_map, internal_genome_id_to_file_path_file, outdir)
