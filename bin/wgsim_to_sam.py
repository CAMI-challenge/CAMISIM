#!/usr/bin/env python3

from Bio import SeqIO
import os
import sys

def create_sam_file(fastq_forward, fastq_reverse, out_sam_file, reference_genome, create_cigar):

    sam_file = open(out_sam_file, 'w')

    counter = 0

    # Paired end reads
    fastq_records1 = SeqIO.parse(fastq_forward, "fastq")
    fastq_records2 = SeqIO.parse(fastq_reverse, "fastq")

    for record1, record2 in zip(fastq_records1, fastq_records2):

        # Process first read of the pair
        read_id1, read_seq1, read_qual1 = record1.id, str(record1.seq), record1.letter_annotations["phred_quality"]
        seq_id, start1, end1, _, _, _ = read_id1.split("_")
        start1_int = int(start1)
        start1 = str(start1_int)
        end1_int = int(end1)
        end1 = str(end1_int)

        # Process second read of the pair
        read_id2, read_seq2, read_qual2 = record2.id, str(record2.seq), record2.letter_annotations["phred_quality"]
        _, end2, start2, _, _, _ = read_id2.split("_")
        start2 = str(int(start2))
        end2 = str(int(end2))

        #Write SAM header in first run
        if (counter == 0):
            # Load the reference genome
            genome = SeqIO.to_dict(SeqIO.parse(reference_genome, "fasta"))
            # Calculate total length
            total_length = sum(len(seq) for seq in genome.values())
            sam_file.write("@SQ\tSN:" + seq_id + "\tLN:" + str(total_length) + "\n")


        if (create_cigar):

            forward_sequence = extract_sequence(reference_genome, seq_id, start1_int, end1_int)
            reverse_sequence = extract_sequence(reference_genome, seq_id, start1_int, end1_int, reverse=True)

            cigar_string1 = get_cigar(forward_sequence, record1)
            cigar_string2 = get_cigar(reverse_sequence, record2)
        else:
            cigar_string1 = str(len(read_seq1)) + "M"
            cigar_string2 = str(len(read_seq2)) + "M"
        

        # Write SAM line for each read.
        sam_file.write(read_id1 + "\t99\t" + seq_id + "\t" + start1 + "\t60\t" + cigar_string1 + "\t=\t" + end1 + "\t" + str((int(end1)+1)-int(start1)) + "\t" + read_seq1 + "\t" + "".join([chr(q + 33) for q in read_qual1]) + "\n")
        sam_file.write(read_id2 + "\t147\t" + seq_id + "\t" + start2 + "\t60\t" + cigar_string2 + "\t=\t" + end2 + "\t" + str(int(end2)-(int(start2)+1)) + "\t" + read_seq2 + "\t" + "".join([chr(q + 33) for q in read_qual2]) + "\n")

        counter = counter + 1

    sam_file.close()


def extract_sequence(genome_file, seq_id, start, end, reverse=False):
    # Load the reference genome
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    # Retrieve the sequence from the reference genome
    sequence = genome[seq_id].seq[start:end]
    
    # If this is a reverse read, return the reverse complement of the sequence
    if reverse:
        sequence = sequence.reverse_complement()
    
    return sequence  

def get_cigar(ref, read):
    cigar = ""
    count = 0
    prev = ""
    for r, s in zip(ref, read):
        if r == s: # Match
            cur = "M"
        elif s == "-": # Deletion
            cur = "D"
        elif r == "-": #Insertion
            cur = "I"
        else:
            cur = "M"  # mismatches are also represented by M

        # Same operation in current and previous base
        if cur == prev or prev == "":
            count += 1
        # Different operation in current and previous base
        else:
            # Add previous operation and the count to CIGAR
            cigar += str(count) + prev
            # Reset count
            count = 1
        prev = cur
    cigar += str(count) + prev
    return cigar

# Main method and entry point of this script
if __name__ == "__main__":

    fastq_forward = sys.argv[1]
    fastq_reverse = sys.argv[2]
    out_sam_file = sys.argv[3]
    reference_genome = sys.argv[4]
    create_cigar_string = str(sys.argv[5])

    if(create_cigar_string == "True"):
        create_cigar = True
    else:
        create_cigar = False    

    create_sam_file(fastq_forward, fastq_reverse, out_sam_file, reference_genome, create_cigar)
