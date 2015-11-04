#! /usr/bin/env python
import os
import sys
import string
import optparse
import fasta
import math


def format(seq, N=60):
	nseg = int(math.ceil(len(seq)/(N+0.0)))
	return '\n'.join([seq[i * N:(i + 1) * N] for i in range(nseg)])
# write into fasta format file 

parser = optparse.OptionParser(version="%prog ")

parser.add_option(
	"-i", "--input", dest="input_fasta", action="store",
	help="name of input file in fasta format")
parser.add_option(
	"-L", "--LibHmm", dest="hmm_path", action="store",
	default="HMM3", help="path of hmm database")
parser.add_option(
	"-o", "--output", dest="out_fname", action="store",
	help="name of output file")
parser.add_option(
	"-k", "--kingdoms", dest="kingdoms", action="store",
	default="arc,bac", help="kingdom used")
#        default="arc,bac,euk", help="kingdom used")

parser.add_option(
	"-r", "--hmmsearch", dest="hmmsearch", action="store", default="hmmsearch",
	help="path to rnammer executable")

parser.add_option(
	"-m", "--moltypes", dest="moltypes", action="store",
	default="lsu,ssu,tsu", help="molecule type detected")
#        default="lsu,ssu,tsu,lsurnammer,ssurnammer,tsurnammer", help="molecule type detected")

parser.add_option(
	"-e", "--Evalue", dest="evalue", action="store", type="float",
	default=0.01, help="evalue cut-off for hmmsearch")
parser.add_option(
	"-p", "--pThreads", dest="p", action="store", type="int",
	default=1, help="number of threads for hmmsearch")


try:
	(options, args) = parser.parse_args()
except:
	parser.print_help()
	sys.exit(1)

if options.input_fasta is None or options.hmm_path is None:
	parser.print_help()
	sys.exit(1)

# os.environ["HMMERDB"] += ":"+os.path.abspath(options.hmm_path)
# print os.environ["HMMERDB"]
out_fname = os.path.abspath(options.out_fname)
out_dir = os.path.dirname(out_fname)
fname = os.path.abspath(options.input_fasta)

tr = string.maketrans("gatcryswkmbdhvnGATCRYSWKMBDHVN", "ctagyrswmkvhdbnCTAGYRSWMKVHDBN")


def rev_record(record):
	return ">" + record.header + "|rev\n" + format(record.sequence[::-1].translate(tr))


records = [rec for rec in fasta.fasta_itr(fname)]
headers = [[rec.header, len(rec.sequence)] for rec in records]

ff = open(out_fname + '.fa', 'w')
for (i, rec) in enumerate(records):
	ff.write('>s' + str(i) + '\n' + format(rec.sequence) + '\n')
	ff.write('>s' + str(i) + '|rev\n' + format(rec.sequence[::-1].translate(tr)) + '\n')
ff.close()
# sys.exit(1)
# a temporary fasta file, use s(int) to easy the parsing


def parse_hmmsearch(kingdom, moltype, src):
	# function to parse hmmsearch output
	resu = []
	data = open(src).readlines()
	# inds = [-1] + [i for (i, x) in enumerate(data[2]) if x == " "]
	# inds = [(inds[j] + 1, inds[j + 1]) for j in range(len(inds) - 1)]
	data = [line for line in data if line[0] != "#"]
	for line in data:
		if not len(line.strip()):
			continue
		[
			read, acc, tlen, qname, qaccr, qlen, seq_evalue, seq_score, seq_bias,
			seq_num, seq_of, dom_cEvalue, dom_iEvalue, dom_score, dom_bias,
			hmm_start, hmm_end, dom_start, dom_end, env_start, env_end] = line.split()[:21]
		# [line[x[0]:x[1]].strip() for x in inds[:21]]
		if string.atof(dom_iEvalue) < options.evalue:
			#            resu.append("\t".join([read, acc, tlen, qname, qaccr, \
			#                    qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of, \
			#                    dom_cEvalue, dom_iEvalue, dom_score, dom_bias, hmm_start, \
			#                    hmm_end, dom_start, dom_end, env_start, env_end]))
			resu.append("\t".join([qname, dom_start, dom_end, read, dom_iEvalue]))

#    print resu[0]
#    print resu[-1]
	return resu


hmm_resu = []
for kingdom in options.kingdoms.split(','):
	for moltype in options.moltypes.split(','):
		# print kingdom, moltype
		hmm_out_fname = "%s.%s_%s.out" % (out_fname, kingdom, moltype)
		dom_out_fname = "%s.%s_%s.dom" % (out_fname, kingdom, moltype)
		cmd = '%s --cpu %d -o %s --domtblout %s -E %g %s/%s_%s.hmm %s' % \
			(options.hmmsearch, options.p, hmm_out_fname, dom_out_fname,
				options.evalue, os.path.abspath(options.hmm_path), kingdom, moltype, out_fname + '.fa')
		#        print cmd
		#        hmm_resu += parse_hmmsearch(os.popen(cmd))
		#        print cmd
		os.system(cmd)
		hmm_resu += parse_hmmsearch(kingdom, moltype, dom_out_fname)
		os.remove(hmm_out_fname)
		os.remove(dom_out_fname)

dict_read2kingdom = {}
for line in hmm_resu:
	[feature_type, r_start, r_end, read, evalue] = line.strip().split('\t')
	read = read.split('|')[0]
	evalue = string.atof(evalue)
	kingdom = feature_type.split('_')[0]
	if read in dict_read2kingdom:
		if evalue < dict_read2kingdom[read][1]:
			dict_read2kingdom[read] = [kingdom, evalue]
	else:
		dict_read2kingdom[read] = [kingdom, evalue]

header = ['##seq_name', 'method', 'feature', 'start', 'end', 'evalue', 'strand', 'frame', 'attribute']
ff = open(out_fname, "w")
dict_rRNA = {
	'arc_lsu': '23S_rRNA', 'arc_ssu': '16S_rRNA', 'arc_tsu': '5S_rRNA',
	'bac_lsu': '23S_rRNA', 'bac_ssu': '16S_rRNA', 'bac_tsu': '5S_rRNA',
	'euk_lsu': '28S_rRNA', 'euk_ssu': '18S_rRNA', 'euk_tsu': '8S_rRNA'}
#             'arc_lsurnammer': '23S_rRNA', 'arc_ssu.rnammer': '16S_rRNA', 'arc_tsurnammer': '5S_rRNA',
#             'bac_lsurnammer': '23S_rRNA', 'bac_ssu.rnammer': '16S_rRNA', 'bac_tsurnammer': '5S_rRNA',
#             'euk_lsurnammer': '28S_rRNA', 'euk_ssu.rnammer': '18S_rRNA', 'euk_tsurnammer': '8S_rRNA'

ff.write('\t'.join(header)+'\n')
for line in hmm_resu:
	#    [kingdom, moltype, read, acc, tlen, qname, qaccr, \
	#         qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of, \
	#         dom_cEvalue, dom_iEvalue, dom_score, dom_bias, hmm_start, \
	#         hmm_end, dom_start, dom_end, env_start, env_end] = line.strip().split('\t')
	[feature_type, r_start, r_end, read, evalue] = line.strip().split('\t')
	if dict_read2kingdom[read.split('|')[0]][0] != feature_type.split('_')[0]:
		continue
	feature_type = dict_rRNA[feature_type]
	if read.endswith('|rev'):
		strand = '-'
		tmp = map(string.atoi, [r_start, r_end])
		pos = string.atoi(read[1:-4])
		header = headers[pos][0]
		L = headers[pos][1]
		[r_end, r_start] = [str(L + 1 - x) for x in tmp]
	else:
		strand = '+'
		pos = string.atoi(read[1:])
		header = headers[pos][0]
	ff.write('\t'.join([header, 'rna_hmm3', 'rRNA', r_start, r_end, evalue, strand, 'NA', feature_type]) + '\n')
ff.close()

os.remove(out_fname + '.fa')
