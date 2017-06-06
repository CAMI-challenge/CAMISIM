#! /usr/bin/env python
import os
import subprocess
import sys
import string
import optparse
import fasta
import math
import tempfile
import shutil


def format(seq, N=60):
	nseg = int(math.ceil(len(seq)/(N+0.0)))
	return '\n'.join([seq[i * N:(i + 1) * N] for i in range(nseg)])
# write into fasta format file

parser = optparse.OptionParser(version="%prog ")

parser.add_option(
	"-i", "--input", dest="input_fasta", action="store",
	help="name of input file in fasta format")

parser.add_option(
	"-r", "--rnammer", dest="rnammer", action="store",
	help="path to rnammer executable")

parser.add_option(
	"-o", "--output", dest="out_prefix_fname", action="store",
	help="name of output file")

parser.add_option(
	"-k", "--kingdoms", dest="kingdoms", action="store",
	default="arc,bac,euk", help="kingdom used")
#        default="arc,bac", help="kingdom used")

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

parser.add_option(
	"-T", "--tmp", dest="temp_folder", action="store", default=None, type=str,
	help="temporary folder")


try:
	(options, args) = parser.parse_args()
except:
	parser.print_help()
	sys.exit(1)

# or options.hmm_path is None
if options.input_fasta is None:
	parser.print_help()
	sys.exit(1)

rnammer = options.rnammer
if rnammer is None:
	print "no executable for rnammer given"
	parser.print_help()
	sys.exit(1)

out_prefix_fname = options.out_prefix_fname
if out_prefix_fname is None:
	parser.print_help()
	sys.exit(1)

# os.environ["HMMERDB"] += ":"+os.path.abspath(options.hmm_path)
# print os.environ["HMMERDB"]
out_prefix_fname = os.path.abspath(out_prefix_fname)
out_dir = os.path.dirname(out_prefix_fname)
fname = os.path.abspath(options.input_fasta)

tr = string.maketrans("gatcryswkmbdhvnGATCRYSWKMBDHVN", "ctagyrswmkvhdbnCTAGYRSWMKVHDBN")


def rev_record(record):
	return ">" + record.header + "|rev\n" + format(record.sequence[::-1].translate(tr))

records = [rec for rec in fasta.fasta_itr(fname)]
headers = [[rec.header, len(rec.sequence)] for rec in records]

in_fasta = out_prefix_fname + '.fa'
ff = open(in_fasta, 'w')
for (i, rec) in enumerate(records):
	ff.write('>s' + str(i) + '\n' + format(rec.sequence) + '\n')
	# ff.write('>s' + str(i) + '|rev\n' + format(rec.sequence[::-1].translate(tr)) + '\n')
ff.close()
# sys.exit(1)
# a temporary fasta file, use s(int) to easy the parsing


def parse_hmmsearch(kingdom, moltype, src):
	# function to parse hmmsearch output
	resu = []
	data = open(src).readlines()
	inds = [-1] + [i for (i, x) in enumerate(data[2]) if x == " "]
	inds = [(inds[j] + 1, inds[j + 1]) for j in range(len(inds) - 1)]
	data = [line for line in data if line[0] != "#"]
	for line in data:
		if not len(line.strip()):
			continue
		[read, acc, tlen, qname, qaccr, qlen, seq_evalue, seq_score, seq_bias,
			seq_num, seq_of, dom_cEvalue, dom_iEvalue, dom_score, dom_bias,
			hmm_start, hmm_end, dom_start, dom_end, env_start, env_end] = line.split()[:21]
		#            [line[x[0]:x[1]].strip() for x in inds[:21]]
		if string.atof(dom_iEvalue) < options.evalue:
			#            resu.append("\t".join([read, acc, tlen, qname, qaccr, \
			#                    qlen, seq_evalue, seq_score, seq_bias, seq_num, seq_of, \
			#                    dom_cEvalue, dom_iEvalue, dom_score, dom_bias, hmm_start, \
			#                    hmm_end, dom_start, dom_end, env_start, env_end]))
			resu.append("\t".join([qname, dom_start, dom_end, read, dom_iEvalue]))

	return resu

dict_rRNA = {
	'arc_lsu': '23S_rRNA', 'arc_ssu': '16S_rRNA', 'arc_tsu': '5S_rRNA',
	'bac_lsu': '23S_rRNA', 'bac_ssu': '16S_rRNA', 'bac_tsu': '5S_rRNA',
	'euk_lsu': '28S_rRNA', 'euk_ssu': '18S_rRNA', 'euk_tsu': '8S_rRNA'}

temp_dir_path = options.temp_folder
if temp_dir_path is None:
	temp_dir_path = tempfile.mkdtemp(suffix="rnammer", dir=None)

hmm_resu = []
for kingdom in options.kingdoms.split(','):
	for moltype in options.moltypes.split(','):
		# print kingdom, moltype
		# hmm_out_fname = "%s.%s_%s.out" % (out_prefix_fname, kingdom, moltype)
		# dom_out_fname = "%s.%s_%s.dom" % (out_prefix_fname, kingdom, moltype)
		out_fasta = out_prefix_fname + "." + dict_rRNA['{0}_{1}'.format(kingdom, moltype)] + ".fna"
		out_gff = out_prefix_fname + "." + dict_rRNA['{0}_{1}'.format(kingdom, moltype)] + ".gff"
		# cmd = "{0} -S {1} -m {2} -gff {3} -f {4} - < {5}".format(rnammer, kingdom, moltype, out_gff, out_fasta, in_fasta)
		cmd = "{0} -S {1} -T {6} -m {2} -gff {3} -f {4} - < {5}".format(rnammer, kingdom, moltype, out_gff, out_fasta, in_fasta, temp_dir_path)
		# subprocess.call(cmd)
		# cmd = "{0} {1} {2} {3} {4} {5} {6} {7}".format(
		# str(bash_rnammer), str(rnammer), str(kingdom), str(moltype), str(out_gff), str(out_fasta), str(in_fasta), str(temp_dir_path))
		# print cmd
		hmm_proc = subprocess.call(cmd, shell=True, bufsize=-1, cwd=temp_dir_path)
		# cmd = 'hmmsearch --cpu %d -o %s --domtblout %s -E %g %s/%s_%s.hmm %s' % \
		#    (options.p, hmm_out_fname, dom_out_fname,
		#        options.evalue, os.path.abspath(options.hmm_path), kingdom, moltype, out_fname + '.fa')
		# hmm_resu += parse_hmmsearch(os.popen(cmd))
		# print cmd
		# os.system(cmd)
		# hmm_resu += parse_hmmsearch(kingdom, moltype, dom_out_fname)
		# os.remove(hmm_out_fname)
		# os.remove(dom_out_fname)

if options.temp_folder is None:
	shutil.rmtree(temp_dir_path)
	os.remove(in_fasta)
