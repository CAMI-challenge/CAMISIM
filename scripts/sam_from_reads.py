import os
import re
from Bio import SeqIO

def write_header(sam_file, sequence_ids):
        with open(sam_file, "w") as samfile:
            samfile.write("@HD\tVN:1.4\tSQ:unsorted\n")
            for prefix in sequence_ids:
                samfile.write("@SQ\tSN:{name}\tLN:{len}\n".format(name=prefix, len=len(sequence_ids[prefix])))

def write_sam(read_file, id_to_cigar_map, reference_path, orig_prefix):
    references = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"))
    fixed_names = {x.split('.',1)[0].replace("_","-"):x for x in references}
    write_sam = os.path.join(read_file.rsplit("/",1)[0], orig_prefix) + ".sam"
    if (not os.path.exists(write_sam)):
        write_header(write_sam, references)
    with open(read_file, 'r') as reads:
        for line in reads:
            if line.startswith('>'):
                name, start, align_status, index, strand, soffset, align_length, eoffset = line.strip().replace(';','_').split('_')
                ref_name = name[1:] # first sign of name is ">"
                ref_name_fixed = fixed_names[ref_name]
                query = ref_name + "-" + start 
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
                    POS = start
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
                    CIGAR = str(len(SEQ)) + "M"
                #    CIGAR = soffset + "I" + CIGAR + str(int(align_length) - int(pos)) + "M" + eoffset + "I"
                #    ###temporarily disabled###
                sam_line = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]
                clen = get_cigar_length(CIGAR)
                with open(write_sam, 'a+') as samfile:
                    samfile.write("\t".join(sam_line) + "\n")
    return references

def get_cigars_nanosim(error_profile):
    errors = {}
    slen = {}
    with open(error_profile, 'r') as ep:
        for line in ep:
            if line.startswith("Seq"):
                continue # header
            name, pos, error_type, length, refseq, qseq = re.split(r'\t|\s{2,}',line) # split at tab and multiple whitespace
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
                for i in range(len(decimals)):
                    length += 10**(len(decimals) - i - 1) * decimals[i]
            decimals = []
    return length

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
            record.description = fixed_names[record.description.split("_")[0]] + "-" + record.description.split("_")[3]
            SeqIO.write(record, fastq, "fastq")
