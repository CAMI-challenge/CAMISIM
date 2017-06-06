import os

def read_reference(reference_path):
	refseq = ""
	with open(reference_path, 'r') as ref:
		for line in ref:
			if not line.startswith('>'): # seq name
				refseq += line.strip() # drop newlines
			else:
				prefix = line[1:].strip()
	return refseq, prefix

def write_header(sam_file, length, sequence_id):
		with open(sam_file, "w") as samfile:
			samfile.write("@HD\tVN:1.4\tSQ:unsorted\n")
			samfile.write("@SQ\tSN:{name}\tLN:{len}\n".format(name=sequence_id, len=length))

# writes all sam in directory and converts fasta reads to fq
def write_all(read_path, dict_id_filepath):
	files = os.listdir(read_path)
	id_to_cigar_map = {}
	for f in files:
		if f.endswith("_error_profile"):	
			prefix = f.rsplit("_",2)[0] # path/to/genomeid_error_profile
			id_to_cigar_map[prefix] = get_cigars(os.path.join(read_path,f))
	for f in files:
		if f.endswith("_reads.fasta"):
			prefix = f.rsplit(".",1)[0].rsplit("_",1)[0] # path/to/genomeid_reads.fasta
			prefix = write_sam(os.path.join(read_path,f),id_to_cigar_map[prefix], dict_id_filepath[prefix], prefix)
			convert_fasta(os.path.join(read_path,f),prefix)

def write_sam(read_file, id_to_cigar_map, reference_path, orig_prefix):
	reference, prefix = read_reference(reference_path) # orig_prefix is prefix without _ in name
	write_sam = os.path.join(read_file.rsplit("/",1)[0], orig_prefix) + ".sam"
	write_header(write_sam, len(reference), prefix)
	with open(read_file, 'r') as reads:
		for line in reads:
			if line.startswith('>'):
				name, start, align_status, index, strand, soffset, align_length, eoffset = line.strip().split('_')
				QNAME = prefix + "-" + index # first sign of name is ">"
				query = name[1:] + "-" + index # nanosim replaces _ by -
				if strand == 'R':
					FLAG = str(16)
				else:
					FLAG = str(0)
				RNAME = prefix
				if align_status == "unaligned": #special cigar/no pos for non-mapping reads
					POS = str(0)
					CIGAR = "*"
				else:
					POS = start
					CIGAR, pos = id_to_cigar_map[query]
				MAPQ = str(255)
				RNEXT = '*'
				PNEXT = '0'
				QUAL = '*' # no quality is given for nanosim
			else:
				SEQ = line.strip()
				TLEN = str(len(SEQ)) 
				if CIGAR != '*': # unmapped bases counted as insertions in read
					CIGAR = soffset + "I" + CIGAR + str(int(align_length) - int(pos)) + "M" + eoffset + "I"
				sam_line = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
				clen = get_cigar_length(CIGAR)
				with open(write_sam, 'a+') as samfile:
					samfile.write("\t".join(sam_line) + "\n")
	return prefix

def get_cigars(error_profile):
	errors = {}
	slen = {}
	with open(error_profile, 'r') as ep:
		for line in ep:
			if line.startswith("Seq"):
				continue # header
			name, pos, error_type, length, refseq, qseq = line.split('\t')
			if error_type == "mis":
				continue # this version ignores mismatches
			seqname = name.split("_")
			seqname = seqname[0] + "-" + seqname[-1] # later on used as sequence name
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

def get_cigar_length(cigar):
	length = 0
	decimals = []
	for char in cigar:
		try:
			x = int(char)
			decimals.append(x)
		except:
			if char != "D":
				for i in xrange(len(decimals)):
					length += 10**(len(decimals) - i - 1) * decimals[i]
			decimals = []
	return length

def convert_fasta(fasta_reads, name):
	out_name = fasta_reads.rsplit("_",1)[0] + ".fq" # /path/to/genomeid_reads.fasta
	with open(fasta_reads,'r') as reads:
		with open(out_name, 'w') as out:
			for line in reads:
				if line.startswith(">"):
					spl = line.strip().split("_")
					index = spl[3]
					out.write("@" + name + "-" + index + '\n')
				else:
					out.write(line)
					out.write("+" + name + "-" + index + '\n')
					out.write("I" * (len(line) - 1) + '\n') # no quality information is available
