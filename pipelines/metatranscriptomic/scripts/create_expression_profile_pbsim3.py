import argparse
import csv
import random
import gffutils
import pyfaidx
from sequence_extractor import extract_sequence

# Set up argument parsing
parser = argparse.ArgumentParser(description="Generate reads and modify SAM files.")
parser.add_argument("--seed", type=int, required=True, help="Random seed.")
parser.add_argument("--db", type=str, required=True, help="Path to the db file.")
parser.add_argument("--fasta_file", type=str, required=True, help="Path to the FASTA file.")
parser.add_argument("--sample_id", type=str, required=True, help="Sample ID.")
parser.add_argument("--genome_id", type=str, required=True, help="Genome ID.")
parser.add_argument("--size", type=float, required=True, help="Size parameter for read count calculation.")
parser.add_argument("--read_length", type=float, required=True, help="Read length.")
parser.add_argument("--child_feature_type", type=str, required=True, help="The child feature type.")
parser.add_argument("--fasta_distribution_file", type=str, required=True, help="Path to the FASTA distribution file.")

args = parser.parse_args()

random.seed(args.seed)

# db_file = 'db_file.db'
# db = gffutils.create_db('${gff_file}', dbfn=db_file)
db_file = gffutils.FeatureDB(args.db)

fasta = pyfaidx.Fasta(args.fasta_file)

total_read_count = 0

exp_f_name = f"sample_{args.sample_id}_{args.genome_id}_expression_profile.tsv"

with open(exp_f_name, 'w') as exp_f:

    with open("transcript_id_to_seq_id.tsv", 'w') as transcript_id_to_seq_id_file:
        # Iterate through each file
        with open(args.fasta_distribution_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            # Iterate through each row in the TSV file
            for row in reader:
                gene_identifier, abundance = row

                read_count = (args.size*(10**9)) * float(abundance) / args.read_length

                if read_count < 1:
                    read_count = 1 if random.random() < read_count else 0

                # gene_identifier_modified = gene_identifier.split("transcript:")[1]

                read_count = round(read_count)

                total_read_count = total_read_count + read_count

                if read_count > 0:

                    seqid, start, combined_sequence = extract_sequence(db_file, fasta, gene_identifier, args.child_feature_type)

                    # exp_f.write(gene_identifier_modified+"\\t0\\t"+str(read_count)+"\\n")
                    exp_f.write(gene_identifier+"\t"+str(read_count)+"\t0\t"+str(combined_sequence)+"\n")

                    # transcript_id_to_seq_id_file.write(gene_identifier_modified+"\\t"+str(seqid)+"\\n")
                    # transcript_id_to_seq_id_file.write(gene_identifier+"\\t"+str(seqid)+"\\n")
                    transcript_id_to_seq_id_file.write(gene_identifier+"\t"+str(seqid)+"\t"+str(start)+"\n")

with open('total_read_count.txt', 'w') as count_file:
    count_file.write(str(total_read_count))