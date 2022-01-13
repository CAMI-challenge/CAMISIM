__author__ = 'majda'

import glob
import os
from Bio import SeqIO


def main(directory="./"):
	dict_sequence_mapping = write_sam(directory)
	merge_fastq(directory, dict_sequence_mapping)
	return dict_sequence_mapping


def merge_fastq(directory, dict_sequence_mapping):
	fastq = os.path.join(directory, "*.fq")
	file_list = sorted(glob.glob(fastq))
	for file_path in file_list:
		orig_file_prefix = os.path.basename(file_path).rsplit("_",1)[0]
		file_path = os.path.join(directory, orig_file_prefix + ".fq")
		with open(file_path, 'a+') as write_handler:
			for record in SeqIO.parse(file_path, "fastq-sanger"):
				record.id = dict_sequence_mapping[orig_file_prefix][record.id]
				record.description = ""
				write_handler.write(record.format("fastq-sanger"))


def write_sam(directory):
	maf = os.path.join(directory, "*.maf")
	list_of_maf_file_path = glob.glob(maf)
	dict_sequence_mapping = {}
	for file_path in list_of_maf_file_path:
		# write sam header
		orig_file_prefix = os.path.basename(file_path).rsplit("_",1)[0]
		sam_file = os.path.join(directory, orig_file_prefix + ".sam")
		with open(sam_file, "w") as samfile:
			samfile.write("@HD\tVN:1.4\tSQ:unsorted\n")
	prefix_to_true_sid = {}
	for file_path in list_of_maf_file_path:
		# get seq_ID
		orig_file_prefix = os.path.basename(file_path).rsplit("_",1)[0]
		prefix = file_path.rsplit(".", 1)[0]
		record = SeqIO.read(prefix + ".ref", "fasta")
		sequence_id = record.id
		prefix_to_true_sid[prefix] = record.id
		sam_file = os.path.join(directory, orig_file_prefix + ".sam")
		# write sam sequence header
		with open(sam_file, "a") as samfile:
			samfile.write("@SQ\tSN:{name}\tLN:{len}\n".format(name=sequence_id, len=len(record.seq)))
	for file_path in list_of_maf_file_path:
		# get seq_ID
		prefix = file_path.rsplit(".", 1)[0]
		# read fastq
		dict_seq = {}
		dict_seq_quality = {}
		fname = prefix + ".fastq"
		if not os.path.isfile(fname):
			fname = prefix + ".fq" # only allow fastq and fq ending
		for record in SeqIO.parse(fname, "fastq-sanger"):
			dict_seq[record.id] = str(record.seq)
			dict_seq_quality[record.id] = SeqIO.QualityIO._get_sanger_quality_str(record)
		sequence_id = prefix_to_true_sid[prefix]
		# write sam sequences
		orig_file_prefix = os.path.basename(file_path).rsplit("_",1)[0]
		sam_file = os.path.join(directory, orig_file_prefix + ".sam")
		with open(sam_file, "a") as samfile:
			dict_of_maf_sid = read_maf(samfile, sequence_id, dict_seq, dict_seq_quality, file_path)
			if orig_file_prefix not in dict_sequence_mapping:
				dict_sequence_mapping[orig_file_prefix] = {}
			for maf_sid in dict_of_maf_sid:
				dict_sequence_mapping[orig_file_prefix][maf_sid] = dict_of_maf_sid[maf_sid]
	return dict_sequence_mapping


def read_maf(samfile, sequence_id, dict_seq, dict_seq_quality, file_path):
	dict_of_maf_sid = {}
	with open(file_path, "r") as maffile:
		n = 0
		index = 0
		for line in maffile:
			if line.startswith("#"):
				continue
			if line.startswith("a"):
				n = 0
				# new multiple_alignment begins
			elif line.startswith("s"):
				# sequence lines begin
				n += 1
				maf_s = line.strip().split()
				# maf_s = line.strip().split(" ")
				if n % 2 == 1:
					POS = maf_s[2]
					SEQ_ref = maf_s[6]
				else:
					SEQ_read = maf_s[6]
					index += 1
					maf_sid = maf_s[1]
					if maf_s[4] == "+":
						FLAG = str(0)
					else:
						FLAG = str(16)
					# if n == 2:
					MAPQ = str(255)
					CIGAR = cigar_code_creation(SEQ_ref, SEQ_read)
					RNEXT = "*"
					PNEXT = "0"
					SEQ = dict_seq[maf_sid]
					QUAL = dict_seq_quality[maf_sid]
					TLEN = str(len(SEQ))
					QNAME = "{sid}-{index}".format(sid=sequence_id, index=index)
					RNAME = sequence_id
					dict_of_maf_sid[maf_sid] = QNAME

					sam_parameter = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
					samfile.write("\t".join(sam_parameter)+"\n")
	return dict_of_maf_sid

# does not count mismatches (X)
def cigar_code_creation(char_ref, char_read):
	cigar = ""
	debug = False
	i = 0
	while (i < len(char_ref) and i < len(char_read)):
		length = 0
		while (i < len(char_ref) and char_ref[i] == '-'):
			i += 1
			length += 1
		if length > 0:
			cigar += "%sI" % length
			continue
		length = 0
		while (i < len(char_read) and char_read[i] == '-'):
			i += 1
			length += 1
		if length > 0:
			cigar += "%sD" % length
			continue
		length = 0
		while (i < len(char_ref) and i < len(char_read) and char_read[i] != '-' and char_ref[i] != '-'):
			i += 1
			length += 1
		if length > 0:
			cigar += "%sM" % length
	return cigar

if __name__ == "__main__":
	main()
