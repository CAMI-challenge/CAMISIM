__author__ = 'hofmann'


#!/bin/bash
#HOME_DIR="/home/user03/pyprojects/taxonomic_classification_plus_otu/"
#t1="/home/user03/output/nobackup/2014_09_24_random1000_ref_silva_hmmer3/16S_rRNA.fna"
#t2="/home/user03/output/nobackup/2014_09_24_random1000_ref_silva_hmmer3/mothur_cluster_list.txt"
#t3="/home/user03/pyprojects/taxonomic_classification_plus_otu/references/"
#t4="0.05"
#t5="no longer used"
#t6="45"


import os
import glob
import shutil
import tempfile

from source.argumenthandler import ArgumentHandler
from source.logger import Logger


def main(options):
	assert isinstance(options, ArgumentHandler)
	mg_cluster = MGCluster(mothur_executable=options.binary_mothur,
						   directory_silva_reference=options.silva_reference_directory,
						   max_processors=options.processors,
						   debug=options._debug_mode)
	return mg_cluster.cluster(marker_gene_fasta=os.path.join(options.project_directory, options.file_mg_16s),
						output_cluster_file=options.processors,
						distance_cutoff=options.distance_cutoff)


class MGCluster(object):
	def __init__(self, mothur_executable, directory_silva_reference, max_processors=1, debug=False, logger=None):
		self._logger = logger
		if not self._logger:
			self._logger = Logger("MGCluster")
		self._mothur_executable = mothur_executable
		self._working_dir = tempfile.mkdtemp()
		self._max_processors = max_processors
		self._debug = debug
		self._ref_silva_distances = os.path.join(directory_silva_reference, "mothur_ref_distances")
		self._ref_silva_names = os.path.join(directory_silva_reference, "mothur_ref_names")  # unique
		self._ref_silva_alignment = os.path.join(directory_silva_reference, "mothur_alignment_ref.fasta")
		self._local_distance = os.path.join(self._working_dir, "ref.align.dist")

	def cluster(self, marker_gene_fasta, output_cluster_file, distance_cutoff):
		local_marker_gene_fasta = self._get_symbolic_link_path(marker_gene_fasta)
		shutil.copy2(self._ref_silva_distances, self._local_distance)

		mothur_cmd = self._get_mothur_cmd(local_marker_gene_fasta, distance_cutoff)
		cmd = "echo -e \"{mothur_cmd}\" | {mothur_executable}".format(mothur_cmd=mothur_cmd,
																	  mothur_executable=self._mothur_executable)
		os.system(cmd)

		find_mask_list = os.path.join(self._working_dir, "*.list")
		list_of_files = glob.glob(find_mask_list)
		if len(list_of_files) == 0:
			self._logger.error("[MGCluster] Clustering failed")
			self._logger.warning("[MGCluster] Remove manually: '{}".format(self._working_dir))
			return False
		elif len(list_of_files) == 1:
			shutil.copy2(list_of_files[0], output_cluster_file)
			if not self._debug:
				shutil.rmtree(self._working_dir)
				self._logger.warning("[MGCluster] Remove manually: '{}".format(self._working_dir))
			self._logger.error("[MGCluster] Clustering success")
			return True
		else:
			self._logger.error("[MGCluster] Clustering with odd result, several files found!")
			self._logger.warning("[MGCluster] Remove manually: '{}".format(self._working_dir))
			return False

	def _get_symbolic_link_path(self, marker_gene_fasta):
		basename = os.path.basename(marker_gene_fasta)
		new_path = os.path.join(self._working_dir, basename)
		os.symlink(marker_gene_fasta, new_path)
		return new_path

	def _get_mothur_cmd(self, marker_gene_fasta, cutoff):
		basename = os.path.basename(marker_gene_fasta)
		filename, extension = os.path.splitext(basename)
		mothur_cmd = """unique.seqs(fasta={mg_fasta})
align.seqs(candidate=current, template={ref_align}, align=gotoh, flip=t, processors={processors})
remove.seqs(accnos={filename}.unique.flip.accnos, fasta=current, name=current)
dist.seqs(oldfasta={ref_align}, column={local_dist}, cutoff={cutoff}, processors={processors}, calc=onegap, countends=F)
merge.files(input={filename}.names-{ref_names}, output={filename}.merged.names)
merge.files(input={filename}.pick.names-{ref_names}, output={filename}.merged.names)
set.current(name={filename}.merged.names, column={local_dist})
cluster(cutoff={cutoff}, method=furthest, precision=100)"""
		return mothur_cmd.format(filename=filename,
								 mg_fasta=marker_gene_fasta,
								 ref_align=self._ref_silva_alignment,
								 ref_names=self._ref_silva_names,
								 local_dist=self._local_distance,
								 processors=self._max_processors,
								 cutoff=cutoff)
