__author__ = 'hofmann'

import os
import sys
import shutil
import tempfile

import source.concat_fasta_on_fasta as concat_fasta_on_fasta
import source.parallel as parallel
from source.argumenthandler import ArgumentHandler


def main(options):
	assert isinstance(options, ArgumentHandler)
	gather_markergenes(detector_exe=os.path.join(options.pipeline_directory, "tools", "16SDetector", "run.py"),
					   mg_type="16S",
					   filename_query_genome_file_paths=options.input_genomes_file,
					   filename_reference_genome_file_paths=options.input_reference_file,
					   filename_reference_marker_genes=options.input_reference_fna_file,
					   out_file_path=os.path.join(options.output_directory, options.file_mg_16s),
					   hmmer=options.hmmer,
					   config_path=options.config_file_path,
					   max_processors=options.processors,
					   debug=options._debug_mode)


def gather_markergenes(detector_exe, mg_type, filename_query_genome_file_paths, filename_reference_genome_file_paths, filename_reference_marker_genes, out_file_path, hmmer, config_path, max_processors=1, debug=False):
	workingdir = tempfile.mkdtemp()
	query_genome_file_paths = parse_genome_file_path_file(filename_query_genome_file_paths)
	if filename_reference_marker_genes is not None:
		reference_genome_file_paths = parse_genome_file_path_file(filename_reference_genome_file_paths)
		query_genome_file_paths.update(reference_genome_file_paths)
	system_link_directory = os.path.join(workingdir, "sym_links")
	local_genome_file_paths = get_local_genome_paths(query_genome_file_paths, system_link_directory)

	# establish symbolic link to fasta files
	for genome_id, src in query_genome_file_paths.iteritems():
		dst = local_genome_file_paths[genome_id]
		if not os.path.exists(src):
			print "File does not exist: '{}".format(src)
			sys.exit(1)
		os.symlink(src, dst)

	#cmd_file_path = os.path.join(workingdir, "cmd_list.txt")
	#create_cmd_file(detector_exe, config_path=config_path, hmmer=hmmer, list_of_fasta=local_genome_file_paths.values(), out_dir=workingdir, cmd_file_path=cmd_file_path)
	cmd_task_list = create_cmd_task_list(detector_exe, config_path=config_path, hmmer=hmmer, list_of_fasta=local_genome_file_paths.values(), out_dir=workingdir)
	parallel.runCmdParallel(cmd_task_list, max_processors)

	tmp_out_file_path = tempfile.mktemp()
	tmp_out_file_bin_path = tmp_out_file_path + "_below_min_length.fna"
	merge_marker_genes_files(local_genome_file_paths, out_file_path, out_bin_file=tmp_out_file_bin_path, mg_type=mg_type, working_dir=workingdir)
	shutil.move(tmp_out_file_path, out_file_path)
	shutil.move(tmp_out_file_bin_path, out_file_path+".bin.fna")

	if filename_reference_marker_genes is not None:
		# append reference genome marker genes
		shutil.copy(out_file_path, out_file_path+".no_ref")
		with open(out_file_path, 'a') as write_handler, open(filename_reference_marker_genes) as read_handler:
			write_handler.writelines(read_handler)

	if not debug:
		shutil.rmtree(workingdir)


def create_cmd_task_list(detector_exe, config_path, hmmer, list_of_fasta, out_dir):
	cmd = "{exe} -c '{config}' -nn -hmmer {hmmer} -i '{input_file}' -out '{out_dir}'"
	cmd_list = [cmd.format(exe=detector_exe, config=config_path, hmmer=hmmer, input_file=file_path, out_dir=out_dir) for file_path in list_of_fasta]
	return [parallel.TaskCmd(cmd, out_dir) for cmd in cmd_list]


def get_local_genome_paths(dict_genome_id_to_path, workingdir="./"):
	if not os.path.exists(workingdir):
		os.makedirs(workingdir)
	dict_genome_id_to_local_path = {}
	for genome_id in dict_genome_id_to_path:
		basename = os.path.basename(dict_genome_id_to_path[genome_id])
		dict_genome_id_to_local_path[genome_id] = os.path.join(workingdir, basename)
	return dict_genome_id_to_local_path


def parse_genome_file_path_file(file_path):
	dict_genome_id_to_path = {}
	with open(file_path) as read_handler:
		for line in read_handler:
			if line.startswith('#'):
				continue
			line = line.strip()
			if len(line) == 0:
				continue
			genome_id, genome_path = line.split('\t')
			dict_genome_id_to_path[genome_id] = genome_path
	return dict_genome_id_to_path


def create_cmd_file(detector_exe, config_path, hmmer, list_of_fasta, out_dir, cmd_file_path):
	with open(cmd_file_path, 'w') as write_handler:
		cmd = "{exe} -c '{config}' -nn -hmmer {hmmer} -i '{input_file}' -out '{out_dir}'"
		for file_path in list_of_fasta:
			write_handler.write(cmd.format(exe=detector_exe,
										   config=config_path,
										   hmmer=hmmer,
										   input_file=file_path,
										   out_dir=out_dir))


# gather all marker genes into single files
def merge_marker_genes_files(dict_genome_id_to_path, out_file_path, out_bin_file, mg_type, working_dir='./'):
	suffixes = {"5S": "5S_rRNA",
				"16S": "16S_rRNA",
				"23S": "23S_rRNA",
				"8S": "8S_rRNA",
				"18S": "18S_rRNA",
				"28S": "28S_rRNA"}

	min_lengths = {"5S": "100",
				"16S": "900",
				"23S": "1000",
				"8S": "100",
				"18S": "900",
				"28S": "1000"}
	suffix = suffixes[mg_type]
	min_length = min_lengths[mg_type]
	assert isinstance(dict_genome_id_to_path, dict)
	for genome_id, genome_path in dict_genome_id_to_path.itervalues():
		input_filename = os.path.basename(genome_path)
		input_filepath = "{prefix}.ids.{suffix}.fna".format(prefix=input_filename,
															suffix=suffix)
		input_file = os.path.join(working_dir, "working", input_filepath)
		concat_fasta_on_fasta.merge(input_file, out_file_path, min_length, unique_id=genome_id, out_bin_file=out_bin_file)
