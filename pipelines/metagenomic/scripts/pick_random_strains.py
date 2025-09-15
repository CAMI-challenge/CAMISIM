import random
import os
from Bio import Phylo

def get_filenames_strains(file_path_template_newick_tree):
    list_of_filenames_strains = []
    tree = Phylo.read(file_path_template_newick_tree, 'newick')
    for leaf in tree.get_terminals():
        prefix = leaf.name
        if prefix.lower() == 'ancestor':
            continue
        list_of_filenames_strains.append(f'{prefix}.fasta')
    return list_of_filenames_strains

def get_empty_row(list_of_column_names, default_value='', as_list=False):
    assert isinstance(default_value, str)
    assert isinstance(as_list, bool)
    if as_list:
        return [default_value] * len(list_of_column_names)
    row = {column_name: default_value for column_name in list_of_column_names}
    return row

def pick_random_strains(amount, genome_id, filenames_strains, NCBI_ID, novelty_category, OTU, out_dir, keep_original=True, filename_prefix='simulated_'):
    column_name_source = 'source'
    column_name_gid = 'genome_ID'
    column_name_ncbi = 'NCBI_ID'
    column_name_novelty_category = 'novelty_category'
    column_name_otu = 'OTU'
    write_source = False

    genome_taxid = NCBI_ID
    
    added_genome_id_to_file_path_genome = {}
    added_meta_table = {}

    genome_id_to_file_path_path = f'genome_id_to_file_path_{genome_id}.tsv'
    meta_table_path = f'meta_table_{genome_id}.tsv'

    with open(genome_id_to_file_path_path, 'w') as genome_id_to_file_path, open(meta_table_path, 'w') as meta_data_file:
        if write_source:
            list_of_column_names = [column_name_gid, column_name_otu, column_name_ncbi, column_name_novelty_category,  column_name_source]
        else:
            list_of_column_names = [column_name_gid, column_name_otu, column_name_ncbi, column_name_novelty_category]

        #meta_data_file.write('\t'.join(list_of_column_names) + '\n')
        amount -= 1

        sample = random.sample(range(0, len(filenames_strains)), amount)
        for index in sample:
            filename = filenames_strains[index]
            name, ext = os.path.splitext(filename)
            new_id = f'{filename_prefix}{genome_id}.{name}'
            source = os.path.join('.', filename)  # Assume current directory
            destination = os.path.join('.', new_id + '.fna')
            os.rename(source, destination)
            added_genome_id_to_file_path_genome[new_id] = destination

            out_path = os.path.join(out_dir, "source_genomes")
            out_path = os.path.join(out_path, new_id + '.fna')
            
            genome_id_to_file_path.write(f'{new_id}\t{out_path}\n')

            row = get_empty_row(list_of_column_names)
            if write_source:
                row[column_name_source] = 'simulated'
            row[column_name_ncbi] = genome_taxid if genome_taxid is not None else 1
            row[column_name_novelty_category] = novelty_category
            row[column_name_otu] = OTU

            added_meta_table[new_id] = row

            line = f'{new_id}'
            for column_name in list_of_column_names:
                if not column_name == column_name_gid:
                    line += f'\t{added_meta_table[new_id][column_name]}'
            meta_data_file.write(line + '\n')

if __name__ == "__main__":
    import sys
    amount = int(sys.argv[1])
    genome_id = sys.argv[2]
    NCBI_ID = int(sys.argv[3])
    novelty_category = sys.argv[4]
    OTU = sys.argv[5]
    strain_simulation_template = sys.argv[6]
    out_dir = sys.argv[7]

    filenames_strains = get_filenames_strains(f"{strain_simulation_template}/template.tree")
    pick_random_strains(amount, genome_id, filenames_strains, NCBI_ID, novelty_category, OTU, out_dir)
