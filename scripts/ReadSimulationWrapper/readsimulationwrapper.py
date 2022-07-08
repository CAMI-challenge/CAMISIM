#!/usr/bin/python

__original_author__ = 'majda'
__author__ = 'peter hofmann'
__version__ = '0.0.5'


import sys
import os
import random
import argparse
import tempfile
from scripts.parallel import TaskCmd, runCmdParallel, reportFailedCmd
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.GenomePreparation.genomepreparation import GenomePreparation
from scripts.ReadSimulationWrapper import sam_from_reads
from scripts.ReadSimulationWrapper import maf_converter
from tools.nanosim_profile import get_mean_fromkde


class ReadSimulationWrapper(GenomePreparation):
    """
    Default Class for all read simulation wrappers

    # TODO: validate genome: description still a problem for art illumina?
    """
    _label = "ReadSimulationWrapper"

    def __init__(
        self, file_path_executable,
        separator='\t', max_processes=1, logfile=None, verbose=True, debug=False, seed=None, tmp_dir=None):
        """
        Constructor

        @param file_path_executable:
        @type file_path_executable: str | unicode
        @param separator: separator to be expected in metadata files
        @type separator: str | unicode
        @param max_processes: Maximum number of processors simulating reads at the same time
        @type max_processes: int 
        @param logfile: file handler or file path to a log file
        @type logfile: file | FileIO | StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param debug: If true temporary files will be kept
        @type debug: bool
        @param seed: Seed used for read simulator, if option available
        @type seed: object
        @param tmp_dir: Directory for storage of temporary files
        @type tmp_dir: int 
        """
        assert self.validate_file(file_path_executable, executable=True)
        assert isinstance(separator, str)
        assert isinstance(max_processes, int)
        assert isinstance(verbose, bool)
        assert isinstance(debug, bool)
        assert seed is None or isinstance(seed, (int, float, str))
        assert tmp_dir is None or isinstance(tmp_dir, str)
        if tmp_dir is not None:
            assert self.validate_dir(tmp_dir)
        else:
            tmp_dir = tempfile.gettempdir()
        self._tmp_dir = self.get_full_path(tmp_dir)
        if seed is not None:
            random.seed(seed)
        super(ReadSimulationWrapper, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
        self._max_processes = max_processes
        self._separator = separator
        self._file_path_executable = file_path_executable
        self._read_length = 150
        self._temporary_files = set()

    def _close(self):
        """
        Remove temporary files
        """
        self._logger = None
        # delete temporary files
        self._remove_temporary_files()

    @staticmethod
    def _get_seed():
        return random.randint(0, sys.maxsize)

    def _remove_temporary_files(self):
        if self._debug:
            return
        while len(self._temporary_files) > 0:
            file_path = self._temporary_files.pop()
            if os.path.isfile(file_path):
                os.remove(file_path)

    # read genome location file
    def _read_genome_location_file(self, file_path):
        """
        Read file with the file paths of gnomes

        @param file_path: File genome id associated with the file path of a genome
        @type file_path: str | unicode

        @return: Dictionary of genome id to file path
        @rtype: dict[str|unicode, str|unicode]
        """
        self._logger.info('Reading genome location file')
        assert self.validate_file(file_path)
        dict_id_file_path = {}
        metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose, separator=self._separator)
        iterator_distributions = metadata_table.parse_file(file_path, as_list=True)
        for genome_id, file_path_genome in iterator_distributions:
            assert genome_id != '', "Invalid genomid: '{}'".format(genome_id)
            assert file_path_genome != '', "Invalid file path: '{}'".format(genome_id)
            assert self.validate_file(file_path_genome), "Invalid file path: '{}'".format(genome_id)

            # check uniqueness
            assert genome_id not in dict_id_file_path, "Genome '{}' not unique in the distribution file!".format(genome_id)
            dict_id_file_path[genome_id] = file_path_genome
        return dict_id_file_path

    def _read_distribution_file(self, file_path):
        """
        Read file with the distribution of a sample

        @param file_path: File genome id associated with the abundance of a genome
        @type file_path: str | unicode

        @return: Dictionary of genome id to file path
        @rtype: dict[str|unicode, float]
        """
        self._logger.info('Reading distribution file')
        assert self.validate_file(file_path)
        dict_id_abundance = {}
        # dict_id_file_path = {}
        metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose, separator=self._separator)
        iterator_distributions = metadata_table.parse_file(file_path, as_list=True)
        # for genome_id, abundance, genome_length, file_path_genome in iterator_distributions:
        abundance_sum = 0.
        for genome_id, abundance in iterator_distributions:
            assert genome_id != '', "Invalid genom id: '{}'".format(genome_id)
            assert abundance != '', "Invalid abundance: '{}'".format(genome_id)
            abundance = float(abundance)
            assert self.validate_number(abundance, zero=True), "Invalid abundance: '{}'".format(genome_id)

            assert genome_id not in dict_id_abundance, "Genome '{}' not unique in the distribution file!".format(genome_id)
            dict_id_abundance[genome_id] = abundance
            abundance_sum += abundance
        dict_id_abundance = {x : dict_id_abundance[x]/abundance_sum for x in dict_id_abundance} # normalise to 1
        return dict_id_abundance

    def get_multiplication_factor(
        self, dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
        file_format="fasta", sequence_type="dna", ambiguous=True):
        """
        A factor is calculated based on total size of a sample to calculate the required covered
        # coverage = abundance * factor
        Files will be validated while the length of sequences are determined.

        @attention min_sequence_length: Sequences that are shorter than the expected fragment size are removed

        @param dict_id_file_path: Dictionary of genome id to file path
        @type dict_id_file_path: dict[str|unicode, str|unicode]
        @param dict_id_abundance: Dictionary of genome id to abundance
        @type dict_id_abundance: dict[str|unicode, float]
        @param total_size: Size of sample in base pairs
        @type total_size: int 
        @param min_sequence_length: Minimum length of a sequence in base pairs
        @type min_sequence_length: int 
        @param file_format: fasta or fastq
        @type file_format: str | unicode
        @param sequence_type: dna or rna or protein
        @type sequence_type: str | unicode
        @param ambiguous: DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
        @type ambiguous: bool

        @return: Factor abundances will be multiplied by
        @rtype: float
        """
        assert isinstance(dict_id_file_path, dict), "Expected dictionary, genome id as key, file path as value"
        assert isinstance(dict_id_abundance, dict), "Expected dictionary, genome id as key, abundance as value"
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(min_sequence_length, int), "Expected natural digit"
        assert isinstance(file_format, str), "Expected file format 'fasta'"
        assert isinstance(sequence_type, str), "Expected sequence type 'rna' or 'dna' or 'protein'"
        assert isinstance(ambiguous, bool)

        relative_size_total = 0
        for genome_id, abundance in dict_id_abundance.items():
            try:
                min_seq_length, genome_length = self.get_sequence_lengths(
                    file_path=dict_id_file_path[genome_id],
                    file_format=file_format,
                    sequence_type=sequence_type,
                    ambiguous=ambiguous,
                    key=None,
                    silent=False)

                if min_seq_length < min_sequence_length:
                    self._logger.info("Genome '{}' has sequences below minimum, creating filtered copy.".format(genome_id))
                    new_file_path = self._remove_short_sequences(
                        dict_id_file_path[genome_id], min_sequence_length, file_format="fasta")
                    dict_id_file_path[genome_id] = new_file_path
                    self._temporary_files.add(new_file_path)

            except IOError as e:
                self._remove_temporary_files()
                raise e

            relative_size = abundance * genome_length
            relative_size_total += relative_size
        return total_size / float(relative_size_total)

    def _remove_short_sequences(self, file_path, min_sequence_length, file_format="fasta"):
        """
        Copies a genome with sequences shorter than a minimum removed.

        @param file_path: File genome id associated with the abundance of a genome
        @type file_path: str | unicode
        @param min_sequence_length: Minimum length of a sequence in base pairs
        @type min_sequence_length: int 
        @param file_format:
        @type file_format: str | unicode

        @return: File path of the genome with removed short sequences
        @rtype: str | unicode
        """
        assert self.validate_file(file_path)
        assert isinstance(min_sequence_length, int), "Expected natural digit"
        assert isinstance(file_format, str), "Expected file format 'fasta'"

        file_path_output = tempfile.mktemp(dir=self._tmp_dir)
        with open(file_path) as stream_input, open(file_path_output, 'w') as stream_output:
            total_base_pairs = self._stream_sequences_of_min_length(
                stream_input, stream_output,
                sequence_min_length=min_sequence_length,
                file_format=file_format
                )
            if total_base_pairs == 0:
                msg = "No valid sequences > {} found!".format(min_sequence_length)
                self._logger.error(msg)
                raise IOError(msg)
        return file_path_output
    
    def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
        """
        Parallel simulation of reads

        @param dict_id_abundance: Dictionary of genome id to abundance
        @type dict_id_abundance: dict[str|unicode, float]
        @param dict_id_file_path: Dictionary of genome id to file path
        @type dict_id_file_path: dict[str|unicode, str|unicode]
        @param factor: Factor abundances will be multiplied by
        @type factor: float | int 
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        """
        self._logger.info("Simulating reads using %s readsimulator..." % self._label)
        assert isinstance(dict_id_file_path, dict), "Expected dictionary, genome id as key, file path as value"
        assert isinstance(dict_id_abundance, dict), "Expected dictionary, genome id as key, abundance as value"
        assert isinstance(factor, (int, float)), "Factor must be numerical"
        assert self.validate_dir(directory_output)

        # add commands to a list of tasks to run them in parallel instead of calling them sequentially
        tasks = []
        for genome_id in dict_id_abundance.keys():
            file_path_input = dict_id_file_path[genome_id]
            abundance = dict_id_abundance[genome_id]
            if abundance == 0:
                continue
            if self._label == "ReadSimulationWgsim" or "ReadSimulationNanosim" in self._label:
                # name "fold_coverage" is misleading for wgsim/nanosim, which use number of reads as input
                fold_coverage = int(round(abundance * factor / self._fragment_size_mean))
            else:
                fold_coverage = abundance * factor
            file_path_output_prefix = os.path.join(directory_output, str(genome_id))
            self._logger.debug("{id}\t{fold_coverage}".format(id=genome_id, fold_coverage=fold_coverage))
            system_command = self._get_sys_cmd(
                file_path_input=file_path_input,
                fold_coverage=fold_coverage,
                file_path_output_prefix=file_path_output_prefix)
            self._logger.debug("SysCmd: '{}'".format(system_command))
            self._logger.info("Simulating reads from {}: '{}'".format(genome_id, file_path_input))
            tasks.append(TaskCmd(system_command))
        list_of_fails = runCmdParallel(tasks, maxProc=self._max_processes)

        if list_of_fails is not None:
            self._logger.error("{} commands returned errors!".format(len(list_of_fails)))
            reportFailedCmd(list_of_fails)
        self._logger.info("Simulating reads finished")

    def _get_sys_cmd(file_path_input, fold_coverage, file_path_output_prefix):
        """
        Abstract class, implement this in your read simulators inheriting from the wrapper
        """
        return
                
class ReadSimulationPBsim(ReadSimulationWrapper):
    """
    Simulate long (PacBio) reads using pbsim
    """

    _label = "ReadSimulationPBsim"


    def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
        super(ReadSimulationPBsim, self).__init__(file_path_executable, **kwargs)
        self._directory_error_profiles = directory_error_profiles
        self._profile = 'standard'
    
    def simulate(
        self, file_path_distribution, file_path_genome_locations, directory_output,
        total_size, profile, fragment_size_mean, fragment_size_standard_deviation):
        """
        Simulate reads based on a given sample distribution

        @param file_path_distribution: File genome id associated with the abundance of a genome
        @type file_path_distribution: str | unicode
        @param file_path_genome_locations: File genome id associated with the file path of a genome
        @type file_path_genome_locations: str | unicode
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        @param total_size: Size of sample in base pairs
        @type total_size: int 
        @param profile: wgsim options: 'errorfree', 'standard'
        @type profile: str | unicode
        @param fragment_size_mean: Size of the fragment of which the ends are used as reads in base pairs
        @type fragment_size_mean: int 
        @param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
        @type fragment_size_standard_deviation: int 
        """
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(fragment_size_mean, int), "Expected natural digit"
        assert isinstance(fragment_size_standard_deviation, int), "Expected natural digit"
        assert total_size > 0, "Total size needs to be a positive number"
        assert fragment_size_mean > 0, "Mean fragments size needs to be a positive number"
        assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
        assert self.validate_dir(directory_output)

        if fragment_size_mean and fragment_size_standard_deviation:
            assert self.validate_number(fragment_size_mean, minimum=1)
            assert self.validate_number(fragment_size_standard_deviation, minimum=0)
            self._fragment_size_mean = fragment_size_mean
            self._fragment_size_standard_deviation = fragment_size_standard_deviation
        # else use pbsim automatic option

        dict_id_abundance = self._read_distribution_file(file_path_distribution)
        dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
        locs = set(dict_id_abundance.keys()) - set(dict_id_file_path.keys())
        assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location %s" % locs
    
        min_sequence_length = 100 # TODO ???
        
        factor = self.get_multiplication_factor(
            dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
            file_format="fasta", sequence_type="dna", ambiguous=True)

        self._logger.debug("Multiplication factor: {}".format(factor))
        self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
        sequence_map = maf_converter.main(directory_output)
        self._fix_extensions(directory_output, sequence_map)

    def _fix_extensions(self, directory_output, sequence_map): # rename fastq to fq
        files = os.listdir(directory_output)
        for f in files:
            if (f.endswith("fastq")):
                oldname = "%s/%s" % (directory_output,f)
                prefix = f.rsplit('_',1)[0] # original name
                with open(oldname,'r') as reads:
                    newname = "%s/%s.fq" % (directory_output,"".join(f.split(".")[:-1]))
                    with open(newname, 'w') as fq: # rename file to fq and rename sequence names
                        for line in reads:
                            if len(line) < 1:
                                continue
                            seq_name = line[1:].strip()
                            if seq_name in sequence_map[prefix]:
                                newline = line[0] + sequence_map[prefix][seq_name] + '\n'
                                fq.write(newline)
                            else:
                                fq.write(line)
    
    def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
        """
        Build system command to be run.

        @param file_path_input: Path to genome fasta file
        @type file_path_input: str | unicode
        @param fold_coverage: coverage of a genome
        @type fold_coverage: int  | float
        @param file_path_output_prefix: Output prefix used by art illumina
        @type file_path_output_prefix: str | unicode

        @return: System command to run art illumina
        @rtype: str | unicode
        """
        assert self.validate_file(file_path_input)
        assert isinstance(fold_coverage, (int, float))
        assert self.validate_dir(file_path_output_prefix, only_parent=True)

        error_profile = os.path.join(self._directory_error_profiles)

        arguments = [
            '--data-type', "CLR",
            '--model_qc', os.path.join(error_profile + "/model_qc_clr"),
            '--depth', str(fold_coverage),
            '--seed', str(self._get_seed()),
            '--prefix', file_path_output_prefix
            ]
        if self._fragment_size_mean is not None:
            arguments.extend([
                '--length-mean', str(self._fragment_size_mean),
                ])
        if self._fragment_size_standard_deviation is not None:
            arguments.extend([
                '--length-sd', str(self._fragment_size_standard_deviation),
            ])

        arguments.extend([
            file_path_input,
            ])
            
        if self._logfile:
            arguments.append(">> '{}'".format(self._logfile))

        cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
        return cmd

class ReadSimulationNanosim3(ReadSimulationWrapper):
    """
    Simulate long reads(Oxford Nanopore) using nanosim
    """

    _label = "ReadSimulationNanosim"

    def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
        super(ReadSimulationNanosim3, self).__init__(file_path_executable, **kwargs)
        self._directory_error_profiles = directory_error_profiles
        # set to None because it has to be calculated from the training data

    def simulate(
        self, file_path_distribution, file_path_genome_locations, directory_output,
        total_size, profile, fragment_size_mean, fragment_size_standard_deviation):
        """
        Simulate reads based on a given sample distribution

        @param file_path_distribution: File genome id associated with the abundance of a genome
        @type file_path_distribution: str | unicode
        @param file_path_genome_locations: File genome id associated with the file path of a genome
        @type file_path_genome_locations: str | unicode
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        @param total_size: Size of sample in base pairs
        @type total_size: int 
        @param profile: wgsim options: 'errorfree', 'standard'
        @type profile: str | unicode
        @param fragment_size_mean: Size of the fragment of which the ends are used as reads in base pairs
        @type fragment_size_mean: int 
        @param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
        @type fragment_size_standard_deviation: int 
        """
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(fragment_size_mean, int), "Expected natural digit"
        assert isinstance(fragment_size_standard_deviation, int), "Expected natural digit"
        assert total_size > 0, "Total size needs to be a positive number"
        assert fragment_size_mean > 0, "Mean fragments size needs to be a positive number"
        assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
        assert self.validate_dir(directory_output)

        dict_id_abundance = self._read_distribution_file(file_path_distribution)
        dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
        locs = set(dict_id_abundance.keys()) - set(dict_id_file_path.keys())
        assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location %s" % locs
        if profile is not None:
            self._profile = profile
        # TODO: this is calculated for every sample due to...reasons right now
        read_length_file = os.path.join(self._directory_error_profiles,profile) + "_aligned_reads.pkl"
        # this is the default filename for Nanosim3 training data
        self._fragment_size_mean = get_mean_fromkde.integrate_mean(read_length_file)
        # nanosim does not use fragment size, this is for getting the correct number of reads
        # this value has been calculated using the script tools/nanosim_profile/get_mean from the values in nanosim_profile

        factor = total_size  # nanosim needs number of reads as input not coverage

        self._logger.debug("Multiplication factor: {}".format(factor))
        self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
        self._sam_from_reads(directory_output, dict_id_file_path)

    def _sam_from_reads(self, directory_output, dict_id_file_path):
        files = os.listdir(directory_output)
        id_to_cigar_map = {}
        for f in files:
            if f.endswith("_error_profile"): # these are the introduced errors by Nanosim
                prefix = f.rsplit("_",3)[0] # get basename (changed in NanoSim3)
                id_to_cigar_map[prefix] = sam_from_reads.get_cigars_nanosim(os.path.join(directory_output,f))
                #os.remove(os.path.join(directory_output,f)) # error_profile files are huge (TODO temporary requirement is still high)
        for f in files:
            if f.endswith("_reads.fasta"):
                prefix = f.rsplit(".",1)[0].rsplit("_",2)[0] #_aligned
                read_file = os.path.join(directory_output,f)
                cigars = id_to_cigar_map[prefix]
                reference_path = dict_id_file_path[prefix]
                sam_from_reads.write_sam(read_file, cigars, reference_path, prefix)
                sam_from_reads.convert_fasta(read_file, reference_path)
                #os.remove(os.path.join(directory_output,f)) # do not store read file twice

    def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
        """
        Build system command to be run.

        @param file_path_input: Path to genome fasta file
        @type file_path_input: str | unicode
        @param fold_coverage: coverage of a genome
        @type fold_coverage: int  | float
        @param file_path_output_prefix: Output prefix used by art illumina
        @type file_path_output_prefix: str | unicode

        @return: System command to run art illumina
        @rtype: str | unicode
        """
        assert self.validate_file(file_path_input)
        assert isinstance(fold_coverage, (int, float))
        assert self.validate_dir(file_path_output_prefix, only_parent=True)

        arguments = [
            'genome',
            '-n', str(fold_coverage),  # rename this, because its not the fold_coverage for wgsim
            '-rg', file_path_input,
            '-o', file_path_output_prefix,
            '-c', os.path.join(self._directory_error_profiles,self._profile),
            '--seed', str(self._get_seed() % 2**32 - 1), # nanosim seed cannot be > 2**32 -1
            '-dna_type linear'
            ]
            
        if self._logfile:
            arguments.append(">> '{}'".format(self._logfile))

        cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
        return cmd

class ReadSimulationNanosim(ReadSimulationWrapper):
    """
    Simulate long reads(Oxford Nanopore) using nanosim
    """

    _label = "ReadSimulationNanosim"

    def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
        super(ReadSimulationNanosim, self).__init__(file_path_executable, **kwargs)
        self._profile = 'standard'

    def simulate(
        self, file_path_distribution, file_path_genome_locations, directory_output,
        total_size, profile, fragment_size_mean, fragment_size_standard_deviation):
        """
        Simulate reads based on a given sample distribution

        @param file_path_distribution: File genome id associated with the abundance of a genome
        @type file_path_distribution: str | unicode
        @param file_path_genome_locations: File genome id associated with the file path of a genome
        @type file_path_genome_locations: str | unicode
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        @param total_size: Size of sample in base pairs
        @type total_size: int 
        @param profile: wgsim options: 'errorfree', 'standard'
        @type profile: str | unicode
        @param fragment_size_mean: Size of the fragment of which the ends are used as reads in base pairs
        @type fragment_size_mean: int 
        @param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
        @type fragment_size_standard_deviation: int 
        """
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(fragment_size_mean, int), "Expected natural digit"
        assert isinstance(fragment_size_standard_deviation, int), "Expected natural digit"
        assert total_size > 0, "Total size needs to be a positive number"
        assert fragment_size_mean > 0, "Mean fragments size needs to be a positive number"
        assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
        assert self.validate_dir(directory_output)
        if profile is not None:
            self._profile = profile

        dict_id_abundance = self._read_distribution_file(file_path_distribution)
        dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
        locs = set(dict_id_abundance.keys()) - set(dict_id_file_path.keys())
        assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location %s" % locs

        self._fragment_size_mean = 7408 # nanosim does not use fragment size, this is for getting the correct number of reads
        # this value has been calculated using the script tools/nanosim_profile/get_mean from the values in nanosim_profile
        factor = total_size  # nanosim needs number of reads as input not coverage

        self._logger.debug("Multiplication factor: {}".format(factor))
        self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
        self._sam_from_reads(directory_output, dict_id_file_path)

    def _sam_from_reads(self, directory_output, dict_id_file_path):
        files = os.listdir(directory_output)
        id_to_cigar_map = {}
        for f in files:
            if f.endswith("_error_profile"): # these are the introduced errors by Nanosim
                prefix = f.rsplit("_",2)[0] # get basename
                id_to_cigar_map[prefix] = sam_from_reads.get_cigars_nanosim(os.path.join(directory_output,f))
                os.remove(os.path.join(directory_output,f)) # error_profile files are huge (TODO temporary requirement is still high)
        for f in files:
            if f.endswith("_reads.fasta"):
                prefix = f.rsplit(".",1)[0].rsplit("_",1)[0]
                read_file = os.path.join(directory_output,f)
                cigars = id_to_cigar_map[prefix]
                reference_path = dict_id_file_path[prefix]
                sam_from_reads.write_sam(read_file, cigars, reference_path, prefix)
                sam_from_reads.convert_fasta(read_file, reference_path)
                os.remove(os.path.join(directory_output,f)) # do not store read file twice

    def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
        """
        Build system command to be run.

        @param file_path_input: Path to genome fasta file
        @type file_path_input: str | unicode
        @param fold_coverage: coverage of a genome
        @type fold_coverage: int  | float
        @param file_path_output_prefix: Output prefix used by art illumina
        @type file_path_output_prefix: str | unicode

        @return: System command to run art illumina
        @rtype: str | unicode
        """
        assert self.validate_file(file_path_input)
        assert isinstance(fold_coverage, (int, float))
        assert self.validate_dir(file_path_output_prefix, only_parent=True)

        arguments = [
            'linear',
            '-n', str(fold_coverage),  # rename this, because its not the fold_coverage for wgsim
            '-r', file_path_input,
            '-o', file_path_output_prefix,
            '-c', "tools/nanosim_profile/ecoli",
            '--seed', str(self._get_seed() % 2**32 - 1) # nanosim seed cannot be > 2**32 -1
            ]
            
        if self._logfile:
            arguments.append(">> '{}'".format(self._logfile))

        cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
        return cmd

# #################
# ReadSimulationWgsim - wgsim Wrapper
# #################

class ReadSimulationWgsim(ReadSimulationWrapper):
    """
    Simulate reads using wgsim
    """
    _label = "ReadSimulationWgsim"
    
    def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
        super(ReadSimulationWgsim, self).__init__(file_path_executable, **kwargs)

    def simulate(
        self, file_path_distribution, file_path_genome_locations, directory_output,
        total_size, profile, fragment_size_mean, fragment_size_standard_deviation):
        """
        Simulate reads based on a given sample distribution

        @param file_path_distribution: File genome id associated with the abundance of a genome
        @type file_path_distribution: str | unicode
        @param file_path_genome_locations: File genome id associated with the file path of a genome
        @type file_path_genome_locations: str | unicode
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        @param total_size: Size of sample in base pairs
        @type total_size: int 
        @param profile: wgsim options: 'errorfree', 'standard'
        @type profile: str | unicode
        @param fragment_size_mean: Size of the fragment of which the ends are used as reads in base pairs
        @type fragment_size_mean: int 
        @param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
        @type fragment_size_standard_deviation: int 
        """
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(fragment_size_mean, int), "Expected natural digit"
        assert isinstance(fragment_size_standard_deviation, int), "Expected natural digit"
        assert total_size > 0, "Total size needs to be a positive number"
        assert fragment_size_mean > 0, "Mean fragments size needs to be a positive number"
        assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
        assert self.validate_dir(directory_output)
        if profile is not None:
            self._profile = profile
        else:
            self._profile = 0.0 # default

        if fragment_size_mean and fragment_size_standard_deviation:
            assert self.validate_number(fragment_size_mean, minimum=1)
            assert self.validate_number(fragment_size_standard_deviation, minimum=0)
            self._fragment_size_mean = fragment_size_mean
            self._fragment_size_standard_deviation = fragment_size_standard_deviation
        else:
            if fragment_size_standard_deviation:
                assert fragment_size_mean is not None, "Both, mean and sd are requires."
            if fragment_size_mean:
                assert fragment_size_standard_deviation is not None, "Both, mean and standard deviation, are required."
        self._logger.info("Simulating with '{}'% errors".format(profile))

        dict_id_abundance = self._read_distribution_file(file_path_distribution)
        dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
        locs = set(dict_id_abundance.keys()) - set(dict_id_file_path.keys())
        assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location %s" % locs

        # min_sequence_length = self._fragment_size_mean - self._fragment_size_standard_deviation
        factor = total_size  # wgsim needs number of reads as input not coverage

        self._logger.debug("Multiplication factor: {}".format(factor))
        self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)

    def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
        """
        Build system command to be run.

        @param file_path_input: Path to genome fasta file
        @type file_path_input: str | unicode
        @param fold_coverage: coverage of a genome
        @type fold_coverage: int  | float
        @param file_path_output_prefix: Output prefix used by art illumina
        @type file_path_output_prefix: str | unicode

        @return: System command to run art illumina
        @rtype: str | unicode
        """
        assert self.validate_file(file_path_input)
        assert isinstance(fold_coverage, (int, float))
        assert self.validate_dir(file_path_output_prefix, only_parent=True)

        arguments = [
            '-d', str(self._fragment_size_mean),
            '-s', str(self._fragment_size_standard_deviation),
            '-N', str(fold_coverage),  # rename this, because its not the fold_coverage for wgsim
            '-1', str(self._read_length),
            '-2', str(self._read_length),
            '-S', str(self._get_seed()),
            ]
        # errors
        arguments.extend([
            '-e', str(self._profile), # base error rate ("sequencing error")
            '-r', "0", # rate of mutations in the genome - unwanted
            '-R', "0", # no indels by default
            # 'X', "0", # this doesnt have to be set to 0 if R is 0 (p for extending indels is 0 if no indels are existent)
            # 'A', MAX_N_RATIO,
        ])
        arguments.extend([
            file_path_input,
            "{}".format(file_path_output_prefix + '1.fq'),
            "{}".format(file_path_output_prefix + '2.fq'),
            "{}".format(file_path_output_prefix + '.sam')
            ])
            
        if self._logfile:
            arguments.append(">> '{}'".format(self._logfile))

        cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
        return cmd
        

# #################
# ReadSimulationArt - Art-Illumina Wrapper
# #################


class ReadSimulationArt(ReadSimulationWrapper):
    """
    Simulate reads using art illumina
    Currently pair-end reads only!
    """
    _label = "ReadSimulationArtIllumina"

    _art_error_profiles = {
        "mi": "EmpMiSeq250R",
        "hi": "EmpHiSeq2kR",
        "hi150": "HiSeq2500L150R",
        "mbarc": "ART_MBARC-26_HiSeq_R"}

    _art_read_length = {
        "mi": 250,
        "hi": 100,
        "hi150": 150,
        "mbarc": 150}

    def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
        super(ReadSimulationArt, self).__init__(file_path_executable, **kwargs)
        # check availability of profiles
        file_names_of_error_profiles = [
            filename+file_end
            for ep, filename in self._art_error_profiles.items()
            for file_end in ['1.txt', '2.txt']
            ]
        assert self.validate_dir(directory_error_profiles, file_names=file_names_of_error_profiles)
        # set default profile
        self._profile = "mbarc"
        self._read_length = self._art_read_length["mbarc"]
        self._directory_error_profiles = directory_error_profiles

    def simulate(
        self, file_path_distribution, file_path_genome_locations, directory_output,
        total_size, profile, fragment_size_mean, fragment_size_standard_deviation,
        profile_filename=None, own_read_length=None):
        """
        Simulate reads based on a given sample distribution

        @param file_path_distribution: File genome id associated with the abundance of a genome
        @type file_path_distribution: str | unicode
        @param file_path_genome_locations: File genome id associated with the file path of a genome
        @type file_path_genome_locations: str | unicode
        @param directory_output: Directory for the sam and fastq files output
        @type directory_output: str | unicode
        @param total_size: Size of sample in base pairs
        @type total_size: int | long
        @param profile: Art illumina error profile: 'low', 'mi', 'hi', 'hi150', 'own'
        @type profile: str | unicode
        @param fragment_size_mean: Size of the fragment of which the ends are used as reads in base pairs
        @type fragment_size_mean: int 
        @param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
        @type fragment_size_standard_deviation: int | long
        @param profile_filename: Optional base name of user-supplied error profile files (without "[1/2].txt").
        @type profile_filename: str | unicode | None
        @param own_read_length: Optional read length for user-supplied error profile.
        @type own_read_length: int | long | None
        """
        assert isinstance(total_size, (float, int)), "Expected natural digit"
        assert isinstance(fragment_size_mean, int), "Expected natural digit"
        assert isinstance(fragment_size_standard_deviation, int), "Expected natural digit"
        assert total_size > 0, "Total size needs to be a positive number"
        assert fragment_size_mean > 0, "Mean fragments size needs to be a positive number"
        assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
        assert self.validate_dir(directory_output)
        # if user specifies own profile, add corresponding parameters
        if profile == "own":
            # sanity checks
            assert own_read_length, "Read length must be given when supplying own profile"
            assert isinstance(own_read_length, (int, long)), "Expected natural digit for read length"
            assert own_read_length > 0, "Read length must be a positive number"
            assert profile_filename, "Profile filename must be given when supplying own profile"
            # sanity check file name
            legal_for_filename = string.ascii_letters + string.digits + '_-./\\'
            assert self.validate_characters(profile_filename, legal_alphabet=legal_for_filename)
            # check if supplied files are present
            own_filenames = [
                profile_filename+file_end
                for file_end in ['1.txt', '2.txt']
            ]
            #assert self.validate_dir(self._directory_error_profiles, file_names=own_filenames)
            for own_file in own_filenames:
                assert self.validate_file(own_file)
            # add user-supplied profiles
            self._art_error_profiles["own"] = profile_filename
            self._art_read_length["own"] = own_read_length
        if profile is not None:
            assert profile in self._art_error_profiles, "Unknown art illumina profile: '{}'".format(profile)
            assert profile in self._art_read_length,  "Unknown art illumina profile: '{}'".format(profile)
            self._profile = profile
        if fragment_size_mean and fragment_size_standard_deviation:
            assert self.validate_number(fragment_size_mean, minimum=1)
            assert self.validate_number(fragment_size_standard_deviation, minimum=0)
            self._fragment_size_mean = fragment_size_mean
            self._fragment_size_standard_deviation = fragment_size_standard_deviation
        else:
            if fragment_size_standard_deviation:
                assert fragment_size_mean is not None, "Both, mean and sd are requires."
            if fragment_size_mean:
                assert fragment_size_standard_deviation is not None, "Both, mean and standard deviation, are required."
        self._logger.info("Using '{}' error profile.".format(profile))

        dict_id_abundance = self._read_distribution_file(file_path_distribution)
        dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
        locs = set(dict_id_abundance.keys()) - set(dict_id_file_path.keys())
        assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location %s" % locs

        min_sequence_length = self._fragment_size_mean - self._fragment_size_standard_deviation
        factor = self.get_multiplication_factor(
            dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
            file_format="fasta", sequence_type="dna", ambiguous=True)

        self._logger.debug("Multiplication factor: {}".format(factor))
        self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)

    def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
        """
        Build system command to be run.

        @param file_path_input: Path to genome fasta file
        @type file_path_input: str | unicode
        @param fold_coverage: coverage of a genome
        @type fold_coverage: int  | float
        @param file_path_output_prefix: Output prefix used by art illumina
        @type file_path_output_prefix: str | unicode

        @return: System command to run art illumina
        @rtype: str | unicode
        """
        assert self.validate_file(file_path_input)
        assert isinstance(fold_coverage, (int, float))
        assert self.validate_dir(file_path_output_prefix, only_parent=True)

        # TODO: mask 'N' default: '-nf 1'
        read_length = self._art_read_length[self._profile]
        error_profile = os.path.join(self._directory_error_profiles, self._art_error_profiles[self._profile])
        arguments = [
            "-sam", "-na",
            "-i '{}'".format(file_path_input),
            "-l", str(read_length),
            "-m", str(self._fragment_size_mean),
            "-s", str(self._fragment_size_standard_deviation),
            "-f", str(fold_coverage),
            "-o '{}'".format(file_path_output_prefix),
            "-1 '{}'".format(error_profile+'1.txt'),
            "-2 '{}'".format(error_profile+'2.txt'),
            ]

        if self._logfile:
            arguments.append(">> '{}'".format(self._logfile))

        # art illumina only accepts integer as seed!
        arguments.append("-rs '{}'".format(self._get_seed()))

        cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
        return cmd


# #################
# #################
#
# FUTURE WORK:
#   ReadSimulationPirs
#   ReadSimulationPBSIM
#
# #################
# #################

# #################
# ReadSimulationPirs
# #################


# class ReadSimulationPirs(ReadSimulationWrapper):
#   _label = "ReadSimulationPirs"
#
#   def simulate(
#       self, file_path_distributions, file_path_genome_locations, directory_output,
#       total_size, read_length, fragment_size_mean, fragment_size_standard_deviation):
#       raise Exception("Not fully implemented yet")
#       assert self.validate_number(read_length, minimum=1)
#       assert self.validate_number(fragment_size_mean, minimum=1)
#       assert self.validate_number(fragment_size_standard_deviation, minimum=0)
#       self._read_length = read_length
#       self._fragment_size_mean = fragment_size_mean
#       self._fragment_size_standard_deviation = fragment_size_standard_deviation
#
#       dict_id_abundance = self._read_distribution_file(file_path_distributions)
#       dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
#
#       # coverage = abundance * factor
#       # factor is calculated based on total size of sample
#       min_sequence_length = self._fragment_size_mean - self._fragment_size_standard_deviation
#       factor = self.get_multiplication_factor(
#           dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
#           file_format="fasta", sequence_type="dna", ambiguous=True)
#       self._logger.debug("Multiplication factor: {}".format(factor))
#       self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
#
#   # start pIRS readsimulator
#   def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
#       # tmp_directory = "./"  # just write files locally
#       # executable = os.path.join(self._directory_read_simulator, "pIRS/SimMetagenome.py")
#       executable = self._file_path_executable
#       # directory_temp = 'pIRS/nobackup/temp_abundance_file.csv'
#       temp_abundance_filename = "temp_abundance_file.csv"
#       for source_data_id in dict_id_abundance.keys():
#           file_path_input = dict_id_file_path[source_data_id]
#           abundance = dict_id_abundance[source_data_id]
#           new_abundance = float(abundance) * factor
#           with open(temp_abundance_filename, 'w') as temp_abundance_file:
#               temp_abundance_file.write(source_data_id+'\t'+str(new_abundance)+'\n')
#           with open(os.path.join(self._directory_read_simulator, 'pIRS/sampleConfig.cfg'), 'r') as config:
#               with open("new_sampleConfig.cfg", 'w') as new_config:
#                   for line in config:
#                       if line.startswith("referenceSeq="):
#                           line = "referenceSeq=" + file_path_input + '\n'
#                       elif line.startswith('frequenciesInfo='):
#                           line = "frequenciesInfo=" + temp_abundance_filename + '\n'
#                       new_config.write(line)
#               os.system("{} -c new_sampleConfig.cfg".format(executable))
#
#
# # #################
# # ReadSimulationPBSIM
# # #################
#
#
# class ReadSimulationPBSIM(ReadSimulationWrapper):
#   _label = "ReadSimulationPBSIM"
#
#   def simulate(
#       self, file_path_distributions, file_path_genome_locations, directory_output,
#       total_size, read_length, fragment_size_mean, fragment_size_standard_deviation):
#       raise Exception("Not fully implemented yet")
#       assert self.validate_number(read_length, minimum=1)
#       assert self.validate_number(fragment_size_mean, minimum=1)
#       assert self.validate_number(fragment_size_standard_deviation, minimum=0)
#       self._read_length = read_length
#       self._fragment_size_mean = fragment_size_mean
#       self._fragment_size_standard_deviation = fragment_size_standard_deviation
#
#       dict_id_abundance = self._read_distribution_file(file_path_distributions)
#       dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
#
#       # coverage = abundance * factor
#       # factor is calculated based on total size of sample
#       min_sequence_length = self._fragment_size_mean - self._fragment_size_standard_deviation
#       factor = self.get_multiplication_factor(
#           dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
#           file_format="fasta", sequence_type="dna", ambiguous=True)
#       self._logger.debug("Multiplication factor: {}".format(factor))
#       self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
#
#   # start PBSIM readsimulator
#   def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
#       # tmp_directory = "./"  # just write files locally
#       # os.chdir(tmp_directory)
#       # executable = os.path.join(self._directory_read_simulator, "pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim")
#       executable = self._file_path_executable
#       for source_data_id in dict_id_abundance.keys():
#           file_path_input = dict_id_file_path[source_data_id]
#           abundance = dict_id_abundance[source_data_id]
#           # file_path_input, abundance, tax_id_predict, genome_name = dict_id_abundance[source_data_id]
#           new_abundance = float(abundance) * factor
#           # if not os.path.exists(tmp_directory+seq_id):
#           #    os.mkdir(tmp_directory+seq_id)
#           prefix = os.path.join(directory_output, str(source_data_id))
#           arguments = [
#               "--data-type", "CLR",
#               "--depth", str(new_abundance),
#               "--prefix", prefix,
#               "--model_qc", os.path.join(self._directory_read_simulator, "pbsim-1.0.3-Linux-amd64/data/model_qc_clr"),
#               file_path_input]
#           # log_file.write("{} {}".format(executable, " ".join(arguments)))
#           os.system("{} {} 2> /dev/null".format(executable, " ".join(arguments)))
#
#       # if self._logger:
#       #    self._logger.error("[ReadSimulation] pbsim currently not active. Please ask developer for more information.")
#       # log_file.write(reader_folder+'pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim --data-type CLR
#       #  --depth '+str(new_abundance)+'
#       #  --model_qc '+reader_folder+'pbsim-1.0.3-Linux-amd64/data/model_qc_clr  '+address)
#       # os.system(reader_folder+'pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim --data-type CLR
#       #  --depth '+str(new_abundance)+'
#       #  --model_qc '+reader_folder+'pbsim-1.0.3-Linux-amd64/data/model_qc_clr  '+address)
#       # TODO: file_path_output = tempfile.mktemp(dir=self._tmp_dir)
#       maf_converter.main(directory_output)


dict_of_read_simulators = {
    "art": ReadSimulationArt,
    "wgsim": ReadSimulationWgsim,
    "nanosim": ReadSimulationNanosim,
    "nanosim3": ReadSimulationNanosim3,
    "pbsim": ReadSimulationPBsim
}
