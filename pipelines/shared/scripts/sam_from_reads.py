#!/usr/bin/env python

import os
import re
from Bio import SeqIO
import argparse


def edit_fastq_transcriptome(read_id_to_seq_id_dict, aligned_reads_path, unaligned_reads_path):
    file_path = os.path.basename(aligned_reads_path).split("_",1)[1]
    with open(file_path, 'a+') as write_handler:

        for record in SeqIO.parse(aligned_reads_path, "fastq-sanger"):

            name, start, align_status, index, strand, soffset, align_length, eoffset = record.id.strip().replace(';','_').rsplit('_', 7)

            record.id = read_id_to_seq_id_dict[name][0] + "-" + index
            record.description = ""
            write_handler.write(record.format("fastq-sanger"))

    file_path = os.path.basename(unaligned_reads_path).split("_",1)[1]
    with open(file_path, 'a+') as write_handler:

        for record in SeqIO.parse(unaligned_reads_path, "fastq-sanger"):

            name, start, align_status, index, strand, soffset, align_length, eoffset = record.id.strip().replace(';','_').rsplit('_', 7)

            record.id = read_id_to_seq_id_dict[name][0] + "-" + index
            record.description = ""
            write_handler.write(record.format("fastq-sanger"))


def sam_from_reads():
    parser = argparse.ArgumentParser(description='Generate SAM from reads.')
    parser.add_argument('error_profile_path', help='Path to error profile file')
    parser.add_argument('aligned_reads_path', help='Path to aligned reads file')
    parser.add_argument('unaligned_reads_path', help='Path to unaligned reads file')
    parser.add_argument('reference_path', help='Path to reference file')
    parser.add_argument('--fastq', action='store_true', help='Read files are in fastq format')
    parser.add_argument('--stdout', action='store_true', help='Write to stdout instead of a file')
    parser.add_argument('--transcriptome', action='store_true', help='Read files are in transcriptomic reads')
    parser.add_argument('--transcript_seq_id_map', help='Path to file holding transcript id to seqid')
    args = parser.parse_args()

    if args.transcriptome and not args.transcript_seq_id_map:
        parser.error("--transcript_seq_id_map is required when --transcriptome is set")

    error_profile_path = args.error_profile_path
    aligned_reads_path = args.aligned_reads_path
    unaligned_reads_path = args.unaligned_reads_path
    reference_path = args.reference_path

    cigars = get_cigars_nanosim(error_profile_path)

    prefix_aligned = os.path.splitext(aligned_reads_path)[0]
    prefix_unaligned = os.path.splitext(unaligned_reads_path)[0]

    if args.transcriptome:
        transcript_seq_id_map_file = args.transcript_seq_id_map

        write_sam_from_transcriptome(aligned_reads_path, cigars, reference_path, prefix_aligned, transcript_seq_id_map_file, True, args.stdout)
        write_sam_from_transcriptome(unaligned_reads_path, cigars, reference_path, prefix_unaligned, transcript_seq_id_map_file, not args.stdout, args.stdout)

        # Read the TSV file into a dictionary
        read_id_to_seq_id_dict = read_tsv_to_dict(transcript_seq_id_map_file)
        edit_fastq_transcriptome(read_id_to_seq_id_dict, aligned_reads_path, unaligned_reads_path)

    else:
        if args.fastq:
            write_sam(aligned_reads_path, cigars, reference_path, prefix_aligned, True, args.stdout)
            write_sam(unaligned_reads_path, cigars, reference_path, prefix_unaligned, not args.stdout, args.stdout)
        else:
            write_sam_from_fasta(aligned_reads_path, cigars, reference_path, prefix_aligned, True, args.stdout)
            convert_fasta(aligned_reads_path, reference_path)

            write_sam_from_fasta(unaligned_reads_path, cigars, reference_path, prefix_unaligned, not args.stdout, args.stdout)
            convert_fasta(unaligned_reads_path, reference_path)


def get_cigars_nanosim(error_profile):

    errors = {}
    with open(error_profile, 'r') as ep:
        for line in ep:
            if line.startswith("Seq"):
                continue # header
            name, pos, error_type, length, refseq, qseq = re.split(r'\t|\s{2,}',line) # split at tab and multiple whitespace
            if error_type == "mis":
                continue # this version ignores mismatches

            # Parse NanoSim identifier safely from the right:
            # <ref>_<start>_<aligned/unaligned>_<index>_<strand>_<soffset>_<align_length>_<eoffset>
            parts = name.rsplit("_", 7)
            ref_name = parts[0]
            start = parts[1]
            index = parts[3]

            # Use a collision-free key: ref-start-index
            seqname = f"{ref_name}-{start}-{index}"

            # It is not sufficient to use "<seqname>-<start_position>" as key (like used above).
            # This is because there are multiple entries in the generated read file with the same combination of sequence name and
            # start position. For more information see https://github.com/bcgsc/NanoSim/issues/151 .
            # This results in a problem with the CIGAR creation and a truncated sam file, because the CIGAR length does not match
            # the sequence length. With using the whole sequence identifier as key and the nanosim version 3.1.0 this can be fixed.
            # seqname = name.replace("_", "-")
            # This fix still triggered the error message:
            # [E::sam_parse1] CIGAR and query sequence are of different length

            if seqname in errors:
                errors[seqname].append((int(pos),error_type,int(length)))
            else:
                errors[seqname] = [(int(pos),error_type,int(length))]
    cigars = {}
    # if ins at pos x with length y, then that ins started at x-y! (see test_error_profile)
    for sequence in errors:
        sorted_errors = sorted(errors[sequence],key=lambda x: x[0])
        ref_len = 0
        CIGAR = ""
        for pos, etype, length in sorted_errors:
            if etype == 'ins':
                if (int(pos) - ref_len > 0):
                    CIGAR += str(int(pos) - ref_len) + "M"
                CIGAR += str(length) + "I"
                ref_len = int(pos)
                length = 0 # only relevant if insertion at the end
            elif etype == 'del':
                if (int(pos) - ref_len > 0):
                    CIGAR += str(int(pos) - ref_len) + "M"
                CIGAR += str(length) + "D"
                ref_len = int(pos) + int(length)
        # if deletion at the end, the number of matches has to be reduces
        cigars[sequence] = (CIGAR, int(pos) + int(length))
    return cigars

# this function writes the sam file from fasta files simulated with nanosim
def write_sam_from_fasta(read_file, id_to_cigar_map, reference_path, orig_prefix, write_header_flag, stdout=False):

    references = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"))
    fixed_names = {x.split('.',1)[0].replace("_","-"):x for x in references}

    sam_file = os.path.join(orig_prefix) + ".sam"

    if write_header_flag:
        write_header(sam_file, references, stdout)
    with open(read_file, 'r') as reads:
        for line in reads:
            if line.startswith('>'):
                name, start, align_status, index, strand, soffset, align_length, eoffset = line.strip().replace(';','_').rsplit('_', 7)
                ref_name = name[1:] # first sign of name is ">"
                ref_name_fixed = fixed_names[ref_name]

                # Must match get_cigars_nanosim() key
                query = f"{ref_name}-{start}-{index}"

                # It is not sufficient to use "<seqname>-<start_position>" as query (like used above).
                # This is because there are multiple entries in the generated read file with the same combination of sequence name and
                # start position. For more information see https://github.com/bcgsc/NanoSim/issues/151 .
                # This results in a problem with the CIGAR creation and a truncated sam file, because the CIGAR length does not matches
                # the sequence length. With using the whole sequence identifier as query and the nanosim version 3.1.0 this can be fixed.
                # query = ref_name + "-" + start  + "-" + align_status  + "-" + index  + "-" + strand  + "-" + soffset  + "-" + align_length  + "-" + eoffset
                # This fix still triggered the error message:
                # [E::sam_parse1] CIGAR and query sequence are of different length

                QNAME = ref_name_fixed + "-" + index
                if strand == 'R':
                    FLAG = str(16)
                else:
                    FLAG = str(0)
                if align_status == "unaligned": #special cigar/no pos for non-mapping reads
                    POS = str(0)
                    CIGAR = "*"
                    RNAME = "*" # treated as unmapped
                    FLAG = str(4)
                else:
                    # Convert 0-based NanoSim start to 1-based SAM POS
                    POS = str(int(start) + 1)
                    try:
                        CIGAR, pos = id_to_cigar_map[query]
                    except KeyError: #sequence did not have any errors
                        CIGAR, pos = "%sM" % align_length, align_length
                    RNAME = ref_name_fixed
                MAPQ = str(255)
                RNEXT = '*'
                PNEXT = '0'
                QUAL = '*'  # no quality for nanosim
            else:
                SEQ = line.strip()
                TLEN = str(len(SEQ))
                if CIGAR != '*': # unmapped bases counted as insertions in read
                    # use real CIGAR for sam file
                    pre_cigar = soffset + "I" if int(soffset) > 0 else ""
                    post_cigar = eoffset + "I" if int(eoffset) > 0 else ""

                    cigar_len = get_cigar_length(pre_cigar + CIGAR + post_cigar)
                    remaining = len(SEQ) - cigar_len
                    if remaining > 0:
                        mid_cigar = str(remaining) + "M"
                        CIGAR = pre_cigar + CIGAR + mid_cigar + post_cigar
                    elif remaining == 0:
                        CIGAR = pre_cigar + CIGAR + post_cigar
                    else:
                        CIGAR = str(len(SEQ)) + "M"
                    ### temporarily disabled ###

                sam_line = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
                if stdout:
                    print("\t".join(sam_line))
                else:
                    with open(sam_file, 'a+') as samfile:
                        samfile.write("\t".join(sam_line) + "\n")
    return references

# this function writes the sam file from fastq files simulated with nanosim
def write_sam(read_file, id_to_cigar_map, reference_path, orig_prefix, write_header_flag, stdout=False):

    references = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"))
    fixed_names = {x.split('.',1)[0].replace("_","-"):x for x in references}
    sam_file = os.path.join(orig_prefix) + ".sam"

    if write_header_flag:
        write_header(sam_file, references, stdout)
    for record in SeqIO.parse(read_file, 'fastq'):
        name, start, align_status, index, strand, soffset, align_length, eoffset = record.id.strip().replace(';','_').rsplit('_', 7)

        ref_name = name
        ref_name_fixed = fixed_names[ref_name]

        # Must match get_cigars_nanosim() key
        query = f"{ref_name}-{start}-{index}"

        # It is not sufficient to use "<seqname>-<start_position>" as query (like used above).
        # This is because there are multiple entries in the generated read file with the same combination of sequence name and
        # start position. For more information see https://github.com/bcgsc/NanoSim/issues/151 .
        # This results in a problem with the CIGAR creation and a truncated sam file, because the CIGAR length does not matches
        # the sequence length. With using the whole sequence identifier as query and the nanosim version 3.1.0 this can be fixed.
        # query = ref_name + "-" + start  + "-" + align_status  + "-" + index  + "-" + strand  + "-" + soffset  + "-" + align_length  + "-" + eoffset
        # This fix still triggered the error message:
        # [E::sam_parse1] CIGAR and query sequence are of different length

        QNAME = ref_name_fixed + "-" + index
        if strand == 'R':
            FLAG = str(16)
        else:
            FLAG = str(0)
        if align_status == "unaligned":
            POS = str(0)
            CIGAR = "*"
            RNAME = "*"
            FLAG = str(4)
        else:
            # Convert 0-based NanoSim start to 1-based SAM POS
            POS = str(int(start) + 1)
            try:
                CIGAR, pos = id_to_cigar_map[query]
            except KeyError:
                CIGAR, pos = "%sM" % align_length, align_length
            RNAME = ref_name_fixed
        MAPQ = str(255)
        RNEXT = '*'
        PNEXT = '0'
        SEQ = str(record.seq)
        QUAL = "".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])
        TLEN = str(len(SEQ))

        if CIGAR != '*': # unmapped bases counted as insertions in read

            # CIGAR = str(len(SEQ)) + "M"

            # use real CIGAR for sam file
            pre_cigar = soffset + "I" if int(soffset) > 0 else ""
            post_cigar = eoffset + "I" if int(eoffset) > 0 else ""

            cigar_len = get_cigar_length(pre_cigar + CIGAR + post_cigar)
            remaining = len(SEQ) - cigar_len
            if remaining > 0:
                mid_cigar = str(remaining) + "M"
                CIGAR = pre_cigar + CIGAR + mid_cigar + post_cigar
            elif remaining == 0:
                CIGAR = pre_cigar + CIGAR + post_cigar
            else:
                CIGAR = str(len(SEQ)) + "M"
            ### temporarily disabled ###

        sam_line = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
        if stdout:
            print("\t".join(sam_line))
        else:
            with open(sam_file, 'a+') as samfile:
                samfile.write("\t".join(sam_line) + "\n")

    return references

def write_header(sam_file, sequence_ids, stdout=False):
    header_lines = ["@HD\tVN:1.4\tSQ:unsorted"]
    header_lines.extend(f"@SQ\tSN:{prefix}\tLN:{len(sequence_ids[prefix])}" for prefix in sequence_ids)

    if stdout:
        print("\n".join(header_lines))
    else:
        with open(sam_file, "w") as samfile:
            samfile.write("\n".join(header_lines) + "\n")


def get_cigar_length(cigar):
    length = 0
    decimals = []
    for char in cigar:
        try:
            x = int(char)
            decimals.append(x)
        except ValueError:
            if char != "D":
                for i in range(len(decimals)):
                    length += 10**(len(decimals) - i - 1) * decimals[i]
            decimals = []
    return length

# This function converts the fasta files simulated with nanosim to fastq
def convert_fasta(read_file, reference_path):
    references = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"))
    fixed_names = {x.split('.',1)[0].replace("_","-"):x for x in references}
    out_name = read_file.rsplit("_",2)[0] + ".fq" # /path/to/genomeid_reads.fasta
    records = SeqIO.parse(read_file, "fasta")
    with open(out_name, 'a+') as fastq:
        for record in records:
            record.letter_annotations["phred_quality"] = [40] * len(record)
            record.id = record.id.replace(";","_")
            record.id = fixed_names[record.id.split("_")[0]] + "-" + record.id.split("_")[3] # this is the index of the read

            # NEW in CAMISIM 2:
            # Skip the description at this point.
            # In CAMISIM 1 the desctiption would be written into the fastq file and removed later in the anonymization (fastastramer.py).

            #record.description = fixed_names[record.description.split("_")[0]] + "-" + record.description.split("_")[3]
            record.description = ''

            SeqIO.write(record, fastq, "fastq")

 # this function writes the sam file from fastq files simulated with nanosim
def write_sam_from_transcriptome(read_file, id_to_cigar_map, reference_path, orig_prefix, transcript_seq_id_map_file, write_header_flag, stdout=False):
    # Read the TSV file into a dictionary
    read_id_to_seq_id_dict = read_tsv_to_dict(transcript_seq_id_map_file)

    references = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"))
    # fixed_names = {x.split('.',1)[0].replace("_","-"):x for x in references}
    sam_file = os.path.join(orig_prefix) + ".sam"

    if write_header_flag:
        write_header(sam_file, references, stdout)
    for record in SeqIO.parse(read_file, 'fastq'):
        name, start, align_status, index, strand, soffset, align_length, eoffset = record.id.strip().replace(';','_').rsplit('_', 7)

        ref_name = name

        ref_name_fixed, start_gene = read_id_to_seq_id_dict[ref_name]
        # ref_name_fixed = read_id_to_seq_id_dict[ref_name]

        # Must match get_cigars_nanosim() key
        query = f"{ref_name}-{start}-{index}"

        # It is not sufficient to use "<seqname>-<start_position>" as query (like used above).
        # This is because there are multiple entries in the generated read file with the same combination of sequence name and
        # start position. For more information see https://github.com/bcgsc/NanoSim/issues/151 .
        # This results in a problem with the CIGAR creation and a truncated sam file, because the CIGAR length does not matches
        # the sequence length. With using the whole sequence identifier as query and the nanosim version 3.1.0 this can be fixed.
        # query = ref_name + "-" + start  + "-" + align_status  + "-" + index  + "-" + strand  + "-" + soffset  + "-" + align_length  + "-" + eoffset
        # This fix still triggered the error message:
        # [E::sam_parse1] CIGAR and query sequence are of different length

        QNAME = ref_name_fixed + "-" + index

        if strand == 'R':
            FLAG = str(16)
        else:
            FLAG = str(0)
        if align_status == "unaligned":
            POS = str(0)
            CIGAR = "*"
            RNAME = "*"
            FLAG = str(4)
        else:
            POS = str(int(start_gene) + int(start))
            try:
                CIGAR, pos = id_to_cigar_map[query]
            except KeyError:
                CIGAR, pos = "%sM" % align_length, align_length
            RNAME = ref_name_fixed
        MAPQ = str(255)
        RNEXT = '*'
        PNEXT = '0'
        SEQ = str(record.seq)
        QUAL = "".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])
        TLEN = str(len(SEQ))

        if CIGAR != '*': # unmapped bases counted as insertions in read
            # use real CIGAR for sam file
            pre_cigar = soffset + "I" if int(soffset) > 0 else ""
            post_cigar = eoffset + "I" if int(eoffset) > 0 else ""

            cigar_len = get_cigar_length(pre_cigar + CIGAR + post_cigar)
            remaining = len(SEQ) - cigar_len
            if remaining > 0:
                mid_cigar = str(remaining) + "M"
                CIGAR = pre_cigar + CIGAR + mid_cigar + post_cigar
            elif remaining == 0:
                CIGAR = pre_cigar + CIGAR + post_cigar
            else:
                CIGAR = str(len(SEQ)) + "M"

        sam_line = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
        if stdout:
            print("\t".join(sam_line))
        else:
            with open(sam_file, 'a+') as samfile:
                samfile.write("\t".join(sam_line) + "\n")

    return references

def read_tsv_to_dict(tsv_file_path):
    """
    Reads a TSV file and stores its contents in a dictionary.

    Parameters:
    tsv_file_path (str): The file path of the TSV file to read.

    Returns:
    dict: A dictionary containing the key-value pairs from the TSV file.
    """

    # Initialize an empty dictionary to store the key-value pairs
    key_value_map = {}

    # Open the TSV file for reading
    with open(tsv_file_path, 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Strip any leading/trailing whitespace from the line and split it by the tab character
            key, value1, value2 = line.strip().split('\t')

            # Add the key-value pair to the dictionary
            key_value_map[key.replace("_","-")] = [value1, value2]

    return key_value_map


if __name__ == "__main__":
    sam_from_reads()
