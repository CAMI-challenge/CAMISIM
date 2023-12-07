import configparser
import sys

def convert_to_nextflow_config(input_file, output_file):
    config = configparser.ConfigParser()
    config.read(input_file)

    with open(output_file, 'w') as nf_config:
        for section in config.sections():
            for key, value in config.items(section):
                # skip commented lines
                if not key.startswith('#'):

                    if (key == 'phase'):
                        continue
                    elif (key == 'max_processors'):
                        continue
                    elif (key == 'dataset_id'):
                        continue
                    elif (key == 'output_directory'):
                        nf_config.write(f"  {'outdir'} = \"{value}\"\n")
                    elif (key == 'temp_directory'):
                        continue
                    elif (key == 'gsa'):
                        if(value == 'True'):
                            nf_config.write(f"  {key} = true\n")
                        else:
                            nf_config.write(f"  {key} = false\n")    
                    elif (key == 'pooled_gsa'):
                        if(not value == 'True'):
                            value = "[]"
                            nf_config.write(f"  {key} = \"{value}\"\n")
                        else:    
                            nf_config.write(f"  {key} = true\n")
                    elif (key == 'anonymous'):
                        if(value == 'True'):
                            nf_config.write(f"  {'anonymization'} = true\n")
                        else:
                            nf_config.write(f"  {'anonymization'} = false\n")
                    elif (key == 'compress'):
                        continue
                    elif (key == 'readsim'):
                        continue
                    elif (key == 'error_profiles'):
                        nf_config.write(f"  {'base_profile_name'} = \"{value}\"\n")
                    elif (key == 'samtools'):
                        continue
                    elif (key == 'profile'):
                        continue
                    elif (key == 'type'):
                        nf_config.write(f"  {key} = \"{value}\"\n")
                    elif (key == 'fragment_size_standard_deviation'):
                        nf_config.write(f"  {'fragment_size_sd'} = {value}\n")
                    elif (key == 'ncbi_taxdump'):
                        nf_config.write(f"  {'ncbi_taxdump_file'} = \"{value}\"\n")
                    elif (key == 'strain_simulation_template'):
                        continue        
                    elif (key == 'metadata'):
                        nf_config.write(f"  {'metadata_file'} = \"{value}\"\n")
                    elif (key == 'id_to_genome_file'):
                        nf_config.write(f"  {'genome_locations_file'} = \"{value}\"\n")
                    elif (key == 'id_to_gff_file'):
                        continue
                    elif (key == 'genomes_total'):
                        continue
                    elif (key == 'genomes_real'):
                        continue
                    elif (key == 'ratio'):
                        continue
                    elif (key == 'mode'):
                        nf_config.write(f"  {key} = \"{value}\"\n")
                    elif (key == 'view'):
                        if(value == 'True'):
                            nf_config.write(f"  {'verbose'} = true\n")
                        else:
                            nf_config.write(f"  {'verbose'} = false\n")
                    elif (key == 'distribution_file_paths'):
                        nf_config.write(f"  distribution_files = \"{value}\"\n")
                    else:
                        nf_config.write(f"  {key} = {value}\n")

        nf_config.write("  biom_profile=\"\"\n")
        nf_config.write("  reference_genomes=\"${projectDir}/tools/assembly_summary_complete_genomes.txt\"\n")
        nf_config.write("  no_replace = true\n")
        nf_config.write("  no_replace = false\n")
        nf_config.write("  additional_references=\"\"\n")
        nf_config.write("  conda.enabled = true\n")
        nf_config.write("  conda.useMamba = true\n")
        nf_config.write("  conda.cacheDir=\"/home/jfunk/conda_cache\"\n")
        nf_config.write("  read_length = 4508\n")
        nf_config.write("  simulate_fastq_directly = false\n")
        nf_config.write("  basecaller = \"guppy\"\n")
        nf_config.write("  profile_read_length=150\n")
        nf_config.write("  base_error_rate = 0\n")
        nf_config.write("  create_cigar = false\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]  # Get input file path from command line argument
    output_file = 'converted_nextflow.config'  # Output file name
    convert_to_nextflow_config(input_file, output_file)
    print(f"Nextflow config file '{output_file}' generated from '{input_file}'.")
