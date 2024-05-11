import argparse
import csv
import random
import gffutils
import pyfaidx
import sys

# Set up argument parsing
parser = argparse.ArgumentParser(description="Generate reads and modify SAM files.")
parser.add_argument("--seed", type=int, required=True, help="Random seed.")
parser.add_argument("--fasta_distribution_file", type=str, required=True, help="Path to the FASTA distribution file.")
parser.add_argument("--gff_file", type=str, required=True, help="Path to the GFF file.")
parser.add_argument("--fasta_file", type=str, required=True, help="Path to the FASTA file.")
parser.add_argument("--sample_id", type=str, required=True, help="Sample ID.")
parser.add_argument("--genome_id", type=str, required=True, help="Genome ID.")
parser.add_argument("--size", type=float, required=True, help="Size parameter for read count calculation.")
parser.add_argument("--read_length", type=float, required=True, help="Read length.")
parser.add_argument("--fragment_size_mean", type=float, required=True, help="Mean fragment size.")
parser.add_argument("--fragment_size_sd", type=float, required=True, help="Fragment size standard deviation.")
parser.add_argument("--profile", type=str, required=True, help="Base profile name.")

args = parser.parse_args()

random.seed(args.seed)

db_file = 'db_file.db'
db = gffutils.create_db(args.gff_file, dbfn=db_file)

fasta = pyfaidx.Fasta(args.fasta_file)

# with open(f"{args.sample_id}_{args.genome_id}_read_counts.tsv", 'w') as expression_file, \
with open(f"{args.sample_id}_{args.genome_id}_commands.sh", 'w') as bash_script:

    with open(args.fasta_distribution_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            gene_identifier, abundance = row
            read_count = ((args.size * (10**9)) * float(abundance) / args.read_length) / 2
            # read_count = max(1, round(read_count * random.random())) ?

            if read_count < 1:
                read_count = 1 if random.random() < read_count else 0

            read_count = round(read_count)

            # expression_file.write(f"{gene_identifier}\t{read_count}\n")

            if read_count >= 1:
                gene = db[gene_identifier]
                output_file = f"sample{args.sample_id}_{args.genome_id}_{gene_identifier}.fa"

                with open(output_file, 'w') as out_f:
                    # out_f.write(f">{gene.seqid}:{gene.start}-{gene.end}\n{gene.sequence(fasta)}\n")
                    out_f.write(">"+str(gene.seqid)+"\n")
                    out_f.write(str(gene.sequence(fasta)))

                seed_for_art = random.randint(0, sys.maxsize)

                bash_script.write(f"art_illumina --paired -sam -na -i {output_file} -l {args.read_length} -c {read_count} -m {args.fragment_size_mean} -s {args.fragment_size_sd} -o sample{args.sample_id}_{args.genome_id}_{gene_identifier} -1 {args.profile}1.txt -2 {args.profile}2.txt -rs {seed_for_art}\n")
                bash_script.write(f"awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}} !/^@/{{ $4=$4+{gene.start}; $8=$8+{gene.start} }}1' sample{args.sample_id}_{args.genome_id}_{gene_identifier}.sam > tmp.sam && mv tmp.sam sample{args.sample_id}_{args.genome_id}_{gene_identifier}.sam\n")
