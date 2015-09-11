__author__ = 'hofmann'

import sys
import os


def read_distr(file_path_distribution, directory_bam=None):
	gid_to_sid_to_length = {}
	gid_to_file_path_bam = {}
	with open(file_path_distribution) as read_handler:
		for line in read_handler:
			explode = line.strip().split('\t')
			gid = explode[0]
			sid = explode[1]
			length = explode[4]
			if gid not in gid_to_sid_to_length:
				gid_to_sid_to_length[gid] = {}
			gid_to_sid_to_length[gid][sid] = long(length) / 6
			file_name = "{}.bam".format(gid)
			if directory_bam:
				gid_to_file_path_bam[gid] = os.path.join(directory_bam, file_name)
	return gid_to_sid_to_length, gid_to_file_path_bam


def main(file_path_distribution, directory_bam, dir_out):
	gid_to_sid_to_length, gid_to_file_path_bam = read_distr(file_path_distribution, directory_bam)
	for gid, file_path in gid_to_file_path_bam.iteritems():
		output = os.path.join(dir_out, "{}".format(gid))
		cmd = "samtools view -h '{bam}' | python '{parser}' '{distr}' '{gid}'| samtools sort - '{output}'; samtools index '{output}.bam'".format(
			bam=file_path,
			parser=os.path.realpath(sys.argv[0]),
			distr=file_path_distribution,
			gid=gid,
			output=output)
		os.system(cmd)


def parse_sam(file_path_distribution, gid):
	gid_to_sid_to_length, gid_to_file_path_bam = read_distr(file_path_distribution)
	for line in sys.stdin:
		if line.startswith('@'):
			sys.stdout.write(line)
			continue
		explode = line.strip().split('\t')
		sid = explode[0].rsplit('-', 1)[0]
		pos = long(explode[3])
		explode[3] = str(((pos-1) % gid_to_sid_to_length[gid][sid]) + 1)
		sys.stdout.write('\t'.join(explode) + '\n')
		sys.stdout.flush()


if __name__ == "__main__":
	if len(sys.argv) == 4:
		main(sys.argv[1], sys.argv[2], sys.argv[3])
	elif len(sys.argv) == 3:
		parse_sam(sys.argv[1], sys.argv[2])
