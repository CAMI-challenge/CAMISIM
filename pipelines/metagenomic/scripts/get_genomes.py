import sys
import biom
import os
import gzip
import urllib.request
import shutil
import re
from numpy import random as np_rand
import sqlite3

try:
    from ete4 import NCBITaxa
except ModuleNotFoundError:
    from ete3 import NCBITaxa

number_of_samples = 0
LEGACY_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
NEW_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'realm', 'domain', 'acellular root', 'cellular root']
RANKS = []
LINUX_FILENAME_MAX_BYTES = 255
MAX_GENOME_PATH_ERRORS_TO_SHOW = 10


"""
Reads a BIOM file and creates map of OTU: lineage, abundance
BIOM file format needs to have a taxonomy field in metadata which contains the taxonomy in the format:
RANK__SCINAME; LOWERRANK_LOWERSCINAME
"""
def read_taxonomic_profile(biom_profile, no_samples = None):
    table = biom.load_table(biom_profile)
    ids = table.ids(axis="observation")
    samples = table.ids()

    if no_samples is None:
        no_samples = len(samples)

    if no_samples is not None and no_samples != len(samples) and no_samples != 1:
        # log.warning("Number of samples (%s) does not match number of samples in biom file (%s)" % (no_samples, len(samples)))
        if no_samples > len(samples):
            no_samples = len(samples)
        #_log.warning("Using the first %s samples" % no_samples)

    profile = {}
    for otu in ids:
        lineage = table.metadata(otu,axis="observation")["taxonomy"]
        try:
            lineage = lineage.split(";") # if no spaces
        except AttributeError:
            pass
        abundances = []
        for sample in samples[:no_samples]:
            abundances.append(table.get_value_by_ids(otu,sample))
        profile[otu] = (lineage, abundances)
    
    return profile


def normalize_genome_key(genome_path):
    genome_path = genome_path.strip()
    if "://" in genome_path:
        return genome_path.rstrip("/")
    return os.path.normpath(genome_path)


def is_remote_path(path):
    return path.startswith("http://") or path.startswith("https://") or path.startswith("ftp://")


def normalize_mapping_path(path):
    if is_remote_path(path):
        return path
    return os.path.abspath(path)


def sanitize_filename_component(value):
    value = value.strip()
    # In taxonomy, square brackets often mark provisional or disputed placement;
    # drop the annotation brackets rather than turning them into separators.
    value = value.replace("[", "").replace("]", "")
    value = re.sub(r'[^A-Za-z0-9._-]+', '_', value)
    value = value.strip('._')
    return value or "genome"


def validate_linux_filename_component(file_name, genome_label):
    if file_name in ("", ".", ".."):
        raise ValueError(
            "Invalid genome output filename '{}' from label '{}'. Linux filenames cannot "
            "be empty, '.' or '..'.".format(file_name, genome_label)
        )
    if "/" in file_name or "\0" in file_name:
        raise ValueError(
            "Invalid genome output filename '{}' from label '{}'. Linux filenames cannot "
            "contain '/' or NUL bytes after sanitization.".format(file_name, genome_label)
        )

    file_name_length = len(file_name.encode("utf-8"))
    if file_name_length > LINUX_FILENAME_MAX_BYTES:
        raise ValueError(
            "Invalid genome output filename '{}' from label '{}'. The filename is {} bytes, "
            "but Linux filename components are limited to {} bytes; use a shorter unique "
            "genome label.".format(
                file_name, genome_label, file_name_length, LINUX_FILENAME_MAX_BYTES
            )
        )


def get_output_genome_file_name(genome_label):
    file_name = "{}.fa".format(sanitize_filename_component(genome_label))
    validate_linux_filename_component(file_name, genome_label)
    return file_name


def allocate_output_genome_file_name(used_output_names, output_entry):
    genome_label = output_entry["label"]
    output_name = get_output_genome_file_name(genome_label)
    stem, extension = os.path.splitext(output_name)
    suffix = 1

    while output_name in used_output_names:
        output_name = "{}_{}{}".format(stem, suffix, extension)
        validate_linux_filename_component(output_name, genome_label)
        suffix += 1

    used_output_names[output_name] = output_entry
    return output_name


def build_output_genome_path(output_dir, output_entry, used_output_names):
    file_name = allocate_output_genome_file_name(used_output_names, output_entry)
    return os.path.join(output_dir, file_name)


def format_genome_record_source(record):
    source_file = record.get("source_file")
    source_line = record.get("source_line")
    if source_file and source_line:
        return "{} line {}".format(source_file, source_line)
    if source_file:
        return source_file
    return "unknown input location"


def format_missing_local_genome_paths_error(missing_entries):
    lines = [
        "{} quality-passing local genome FASTA path(s) do not exist. Fix the "
        "input genome list paths or make these files visible to the worker before "
        "running CAMISIM.".format(len(missing_entries))
    ]
    for taxid, record in missing_entries[:MAX_GENOME_PATH_ERRORS_TO_SHOW]:
        lines.append(
            "- {}: taxid '{}', label '{}', path '{}'".format(
                format_genome_record_source(record),
                taxid,
                record.get("label", ""),
                record.get("path", "")
            )
        )
    if len(missing_entries) > MAX_GENOME_PATH_ERRORS_TO_SHOW:
        lines.append(
            "... {} more missing local genome path(s) not shown.".format(
                len(missing_entries) - MAX_GENOME_PATH_ERRORS_TO_SHOW
            )
        )
    return "\n".join(lines)


def validate_local_genome_paths(genomes_map):
    missing_entries = []
    for taxid in genomes_map:
        for record in genomes_map[taxid]["genomes"]:
            if not record.get("quality_pass", True):
                continue
            path = record["path"]
            if not is_remote_path(path) and not os.path.exists(path):
                missing_entries.append((taxid, record))

    if missing_entries:
        raise ValueError(format_missing_local_genome_paths_error(missing_entries))


def get_table_columns(header_parts, aliases):
    columns = {}
    lowered = [part.strip().lower() for part in header_parts]
    for column_name, options in aliases.items():
        for option in options:
            if option in lowered:
                columns[column_name] = lowered.index(option)
                break
    return columns


def iter_tabular_rows(file_path, aliases, default_indices, required_columns, optional_defaults = None):
    optional_defaults = optional_defaults or {}

    def build_row(parts, column_map, line, line_number):
        row = {}
        for column in required_columns:
            index = column_map[column]
            if len(parts) <= index:
                raise ValueError("Malformed row in '{}': {}".format(file_path, line.rstrip('\n')))
            row[column] = parts[index].strip()

        for column, default_value in optional_defaults.items():
            index = column_map.get(column)
            if index is not None and len(parts) > index:
                row[column] = parts[index].strip()
            else:
                row[column] = default_value
        row["__line_number"] = line_number
        return row

    with open(file_path, 'r') as stream:
        first_line = stream.readline()
        if not first_line:
            return

        first_parts = first_line.rstrip('\n').split('\t')
        header_columns = get_table_columns(first_parts, aliases)
        if all(column in header_columns for column in required_columns):
            column_map = {column: header_columns[column] for column in required_columns}
            for column in optional_defaults:
                if column in header_columns:
                    column_map[column] = header_columns[column]
        else:
            column_map = {column: default_indices[column] for column in required_columns}
            required_max_index = max(column_map.values())
            if len(first_parts) <= required_max_index:
                raise ValueError("Malformed row in '{}': {}".format(file_path, first_line.rstrip('\n')))
            for column in optional_defaults:
                index = default_indices.get(column)
                if index is not None and len(first_parts) > index:
                    column_map[column] = index
            yield build_row(first_parts, column_map, first_line, 1)

        for line_number, line in enumerate(stream, start=2):
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            yield build_row(parts, column_map, line, line_number)


def get_required_quality_column(columns, aliases):
    for alias in aliases:
        if alias in columns:
            return columns[alias]
    raise ValueError("Required column(s) missing from additional genome quality file: %s" % ", ".join(aliases))


def read_additional_genomes_quality(quality_file):
    quality_by_path = {}
    with open(quality_file, 'r') as quality_stream:
        header = quality_stream.readline().rstrip('\n').split('\t')
        if not header or header == ['']:
            raise ValueError("Additional genome quality file '%s' is empty" % quality_file)

        columns = {column.strip().lower(): idx for idx, column in enumerate(header)}
        idx_path = get_required_quality_column(columns, ['genome_path', 'path', 'genome'])
        idx_completeness = get_required_quality_column(columns, ['completeness'])
        idx_contamination = get_required_quality_column(columns, ['contamination'])
        idx_num_contigs = get_required_quality_column(columns, ['num_contigs', 'contigs'])
        max_idx = max(idx_path, idx_completeness, idx_contamination, idx_num_contigs)

        for line in quality_stream:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) <= max_idx:
                raise ValueError("Malformed row in additional genome quality file '%s': %s" % (quality_file, line))

            genome_key = normalize_genome_key(parts[idx_path])
            if genome_key in quality_by_path:
                raise ValueError("Duplicate quality entry for additional genome '%s'" % parts[idx_path])

            quality_by_path[genome_key] = {
                "completeness": float(parts[idx_completeness]),
                "contamination": float(parts[idx_contamination]),
                "num_contigs": int(float(parts[idx_num_contigs]))
            }
    return quality_by_path


def meets_additional_quality_requirements(quality_entry, max_contamination, min_completeness, max_num_contigs):
    return (
        quality_entry["contamination"] <= max_contamination and
        quality_entry["completeness"] >= min_completeness and
        quality_entry["num_contigs"] <= max_num_contigs
    )


"""
Reads list of available genomes in the (tsv) format:
NCBI_ID Genome_Label ftp_path
Additional files might be provided with either:
NCBI_ID Genome_Label genome_path
or:
NCBI_ID Genome_Label genome_path novelty_category
were path might either be online or offline/local
Copied genome files are named <sanitized_genome_label>.fa. If multiple selected
genomes produce the same output filename, later files get _1, _2, ...
suffixes before .fa.
"""
def read_genomes_list(genomes_path, additional_file = None, additional_quality_file = None,
                      max_contamination = 5, min_completeness = 95, max_num_contigs = 200):
    genomes_map = {}
    total_genomes = 0
    additional_quality = {}
    reference_aliases = {
        "taxid": ["ncbi_id", "ncbi_taxid", "taxid"],
        "name": ["scientific_name", "sci_name", "unique_taxon_name", "name"],
        "path": ["ftp_path", "genome_path", "path", "ftp"]
    }
    additional_aliases = {
        "taxid": ["ncbi_id", "ncbi_taxid", "taxid"],
        "name": ["scientific_name", "sci_name", "unique_taxon_name", "name"],
        "path": ["genome_path", "path", "ftp_path", "ftp"],
        "novelty": ["novelty_category", "novelty"]
    }

    if additional_file is not None:
        if additional_quality_file is None:
            raise ValueError("Parameter additional_genomes_quality_file is required when additional_references is set")
        additional_quality = read_additional_genomes_quality(additional_quality_file)
        for row in iter_tabular_rows(
            additional_file,
            additional_aliases,
            {"taxid": 0, "name": 1, "path": 2, "novelty": 3},
            ["taxid", "name", "path"],
            {"novelty": "known_strain"}
        ):
            ncbi_id = row["taxid"]
            sci_name = row["name"]
            path = row["path"]
            novelty = row["novelty"]
            path_key = normalize_genome_key(path)
            if path_key not in additional_quality:
                raise ValueError("Missing quality entry for additional genome '%s' in '%s'" % (path, additional_quality_file))
            quality_entry = additional_quality[path_key]
            genome_record = {
                "path": path,
                "original_path": path,
                "source_file": additional_file,
                "source_line": row["__line_number"],
                "novelty": novelty,
                "source": "additional",
                "label": sci_name,
                "quality_pass": meets_additional_quality_requirements(
                    quality_entry,
                    max_contamination,
                    min_completeness,
                    max_num_contigs
                )
            }
            if ncbi_id in genomes_map:
                genomes_map[ncbi_id]["genomes"].append(genome_record)
            else:
                genomes_map[ncbi_id] = {"scientific_name": sci_name, "genomes": [genome_record]}
            total_genomes += 1
    for row in iter_tabular_rows(
        genomes_path,
        reference_aliases,
        {"taxid": 0, "name": 1, "path": 2},
        ["taxid", "name", "path"]
    ):
        ncbi_id = row["taxid"]
        sci_name = row["name"]
        ftp = row["path"]
        http = ftp.replace("ftp://","http://") # not using ftp address but http (proxies)
        genome_record = {
            "path": http,
            "original_path": ftp,
            "source_file": genomes_path,
            "source_line": row["__line_number"],
            "novelty": 'known_strain',
            "source": "reference",
            "label": sci_name,
            "quality_pass": True
        }
        if ncbi_id in genomes_map:
            genomes_map[ncbi_id]["genomes"].append(genome_record)
        else:
            genomes_map[ncbi_id] = {"scientific_name": sci_name, "genomes": [genome_record]} # sci_name is always the same for same taxid (?)
        total_genomes += 1
    return genomes_map, total_genomes


def get_cached_rank(tax_id, ncbi, rank_cache):
    if tax_id not in rank_cache:
        rank_cache[tax_id] = ncbi.get_rank([tax_id]).get(tax_id)
    return rank_cache[tax_id]


def resolve_name_to_allowed_taxid(name, ncbi, name_taxid_cache, rank_cache):
    if name in name_taxid_cache:
        return name_taxid_cache[name]

    candidate_names = [name]
    fallback_name = name.split()[0]
    if fallback_name and fallback_name != name:
        candidate_names.append(fallback_name)

    resolved_taxid = None
    for candidate_name in candidate_names:
        if candidate_name in name_taxid_cache:
            resolved_taxid = name_taxid_cache[candidate_name]
        else:
            mapping = ncbi.get_name_translator([candidate_name])
            if candidate_name in mapping:
                taxid = mapping[candidate_name][0]
                rank = get_cached_rank(taxid, ncbi, rank_cache)
                resolved_taxid = taxid if rank in RANKS else None
            else:
                resolved_taxid = None
            name_taxid_cache[candidate_name] = resolved_taxid

        if resolved_taxid is not None:
            break

    name_taxid_cache[name] = resolved_taxid
    return resolved_taxid

"""
Given all available genomes, creates a map sorted by ranks of available genomes on that particular rank, ordered by their ncbi ids
"""
def get_genomes_per_rank(genomes_map, ncbi, rank_cache):
    per_rank_map = {}
    genome_bucket_map = {}
    for rank in RANKS:
        per_rank_map[rank] = {}
    for genome in genomes_map:
        try:
            lineage = ncbi.get_lineage(genome) # this might contain some others ranks than ranks
            ranks_lin = ncbi.get_rank(lineage)
            genome_records = []
            for genome_record in genomes_map[genome]["genomes"]:
                bucketed_genome_record = {
                    "path": genome_record["path"],
                    "original_path": genome_record["original_path"],
                    "source_file": genome_record["source_file"],
                    "source_line": genome_record["source_line"],
                    "genome_id": genome,
                    "novelty": genome_record["novelty"],
                    "source": genome_record["source"],
                    "label": genome_record["label"],
                    "quality_pass": genome_record["quality_pass"]
                }
                genome_records.append(bucketed_genome_record)
                genome_bucket_map[id(bucketed_genome_record)] = []
            for tax_id in lineage: # go over the lineage
                try:
                    check_rank = ranks_lin[tax_id]
                except KeyError:
                    continue
                rank_cache[tax_id] = check_rank
                if check_rank in per_rank_map: # if we are a legal rank
                    rank_map = per_rank_map[ranks_lin[tax_id]]
                    if tax_id not in rank_map: # tax id doesn't have a genome yet
                        rank_map[tax_id] = []
                    bucket = rank_map[tax_id]
                    for genome_record in genome_records:
                        bucket.append(genome_record)
                        genome_bucket_map[id(genome_record)].append(bucket)
        except ValueError as e:
#           log.warning(e)
            print(e) # ToDo
    return per_rank_map, genome_bucket_map

"""
Sorts the otus in the profile by abundance
"""
def sort_by_abundance(profile):
    sorted_keys = []
    for otu in profile:
        lineage, abundances = profile[otu]
        avg_abundance = sum(abundances)/len(abundances)
        sorted_keys.append((avg_abundance, otu)) # average abundance
    sorted_keys = sorted(sorted_keys, reverse=True)
    return [key for ab,key in sorted_keys]


"""
Given a BIOM lineage, create a NCBI tax id lineage
"""
def transform_lineage(lineage, ncbi):
    return transform_lineage_with_cache(lineage, ncbi, {}, {})


def transform_lineage_with_cache(lineage, ncbi, name_taxid_cache, rank_cache):
    new_lineage = []
    for member in lineage:
        name = member.split("__")[-1] # name is on the right hand side
        if len(name) == 0:
            continue
        taxid = resolve_name_to_allowed_taxid(name, ncbi, name_taxid_cache, rank_cache)
        if taxid is not None:
            new_lineage.append(taxid) # should contain only one element
    return new_lineage[::-1] # invert list, so lowest rank appears first (last in BIOM)


def truncated_geometric(p, l, u, max_tries=100000):
    if not (0 < p <= 1) or l > u:
        raise ValueError("Invalid parameters: p must be in (0,1) and l <= u.")

    for _ in range(max_tries):
        k = np_rand.geometric(p)
        if l <= k <= u:
            return k

    raise RuntimeError(
        f"Failed to draw a sample in the range [{l}, {u}] within {max_tries} tries. "
        f"The probability of a draw falling in this range may be too low for p={p}."
    )


def validate_strain_bounds(min_strains, max_strains):
    if min_strains < 1 or max_strains < 1:
        raise ValueError(
            f"Invalid strain bounds: min_strains={min_strains}, max_strains={max_strains}. "
            "Both values must be at least 1."
        )
    if min_strains > max_strains:
        raise ValueError(
            f"Invalid strain bounds: min_strains={min_strains}, max_strains={max_strains}. "
            "min_strains must be less than or equal to max_strains."
        )


def draw_num_strains(min_strains, max_strains):
    validate_strain_bounds(min_strains, max_strains)

    if min_strains == max_strains:
        return min_strains

    # The historical geometric draw becomes pathological for max_strains <= 2:
    # with max=2 it can only return 1, and with max=1 it is invalid.
    # Use a bounded uniform draw for these small ranges and preserve the
    # previous geometric behaviour for larger ranges.
    if max_strains <= 2:
        return int(np_rand.randint(min_strains, max_strains + 1))

    return truncated_geometric(2. / max_strains, min_strains, max_strains)


def set_abundances(otu_genome_map, abundances, used_genomes, otu, tax_id, mu, sigma):
    log_normal_vals = np_rand.lognormal(mu, sigma, len(used_genomes))
    sum_log_normal = sum(log_normal_vals)
    for i, genome in enumerate(used_genomes):
        otu_id = otu + "." + str(i)
        otu_genome_map[otu_id] = (tax_id, genome["genome_id"], genome["path"], genome["original_path"], genome["source_file"], genome["source_line"], genome["novelty"], genome["label"], [])  # taxid, genomeid, copy path, original path, source file, source line, novelty, label, abundances per sample
        relative_abundance = log_normal_vals[i] / sum_log_normal
        for abundance in abundances:  # calculate abundance per sample
            current_abundance = relative_abundance * abundance
            otu_genome_map[otu_id][-1].append(current_abundance)


def sample_genomes(genomes, amount):
    if amount <= 0 or not genomes:
        return []
    draw_count = min(len(genomes), amount)
    used_indices = np_rand.choice(len(genomes), draw_count, replace=False)
    return [genomes[i] for i in used_indices]


def get_candidate_groups(available_genomes, prioritize_additional_genomes):
    high_quality_additional = [genome for genome in available_genomes if genome["source"] == "additional" and genome["quality_pass"]]
    reference_genomes = [genome for genome in available_genomes if genome["source"] == "reference"]
    low_quality_additional = [genome for genome in available_genomes if genome["source"] == "additional" and not genome["quality_pass"]]

    if prioritize_additional_genomes:
        candidate_groups = [high_quality_additional, reference_genomes]
    else:
        candidate_groups = [high_quality_additional + reference_genomes]

    candidate_groups.append(low_quality_additional)
    return candidate_groups


def remove_genomes_from_rank_map(used_genomes, genome_bucket_map):
    for genome in used_genomes:
        bucket_refs = genome_bucket_map.pop(id(genome), None)
        if not bucket_refs:
            continue
        for bucket in bucket_refs:
            try:
                bucket.remove(genome)
            except ValueError:
                continue


def pick_genomes(tax_id, max_rank, min_strains, max_strains, per_rank_map, ncbi, no_replace,
                 prioritize_additional_genomes, rank_cache, genome_bucket_map):
    rank = get_cached_rank(tax_id, ncbi, rank_cache)
    if RANKS.index(rank) > RANKS.index(max_rank):
        return []
    available_genomes = per_rank_map[rank].get(tax_id)
    if not available_genomes:
        return []
    strains_to_draw = draw_num_strains(min_strains, max_strains)
    used_genomes = []
    for candidate_group in get_candidate_groups(available_genomes, prioritize_additional_genomes):
        remaining = strains_to_draw - len(used_genomes)
        if remaining <= 0:
            break
        used_genomes.extend(sample_genomes(candidate_group, remaining))

    if no_replace: # sampling without replacement:
        remove_genomes_from_rank_map(used_genomes, genome_bucket_map)

    return used_genomes


"""
Given the OTU to lineage/abundances map and the genomes to lineage map, create map otu: taxid, genome, abundances
"""
def map_otus_to_genomes(tax_profile, per_rank_map, mu, sigma, min_strains, max_strains, no_replace,
                        max_genomes, max_rank, prioritize_additional_genomes, ncbi, rank_cache,
                        genome_bucket_map):
    unmatched_otus = []
    matched_otus = set()
    otu_genome_map = {}
    sorted_otus = sort_by_abundance(tax_profile)
    genome_set_size = 0

    otu_to_lineage = {}
    otus_list = []
    max_lineage_len = 0
    name_taxid_cache = {}
    for otu in sorted_otus:
        lin, _ = tax_profile[otu]
        otu_to_lineage[otu] = transform_lineage_with_cache(lin, ncbi, name_taxid_cache, rank_cache)
        if len(otu_to_lineage[otu]) == 0:
            unmatched_otus.append(otu)
            continue
        otus_list.append(otu)
        max_lineage_len = max(max_lineage_len, len(otu_to_lineage[otu]))

    for idx_lineage in range(max_lineage_len + 1):
        for otu in otus_list:
            if genome_set_size >= max_genomes and no_replace:  # cancel if no genomes are available anymore
                break
            if otu in matched_otus:
                continue

            try:
                tax_id = otu_to_lineage[otu][idx_lineage]
            except IndexError:
                continue

            used_genomes = pick_genomes(
                tax_id, max_rank, min_strains, max_strains, per_rank_map, ncbi, no_replace,
                prioritize_additional_genomes, rank_cache, genome_bucket_map
            )
            if not used_genomes:
                continue

            genome_set_size += len(used_genomes)  # how many genomes are used

            matched_otus.add(otu)
            set_abundances(otu_genome_map, tax_profile[otu][1], used_genomes, otu, tax_id, mu, sigma)

    for otu in otus_list:
        if otu not in matched_otus:
            unmatched_otus.append(otu)

    return otu_genome_map, unmatched_otus

def fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile, debug, prioritize_additional_genomes):
    genomes = {}
    for rank in per_rank_map:
        for taxid in per_rank_map[rank]:
            for genome in per_rank_map[rank][taxid]:
                genome_key = (genome["path"], genome["genome_id"])
                if genome_key not in genomes:
                    genomes[genome_key] = {
                        "tax_id": taxid,
                        "genome_id": genome["genome_id"],
                        "path": genome["path"],
                        "original_path": genome["original_path"],
                        "source_file": genome["source_file"],
                        "source_line": genome["source_line"],
                        "novelty": genome["novelty"],
                        "source": genome["source"],
                        "label": genome["label"],
                        "quality_pass": genome["quality_pass"]
                    }

    ordered_genomes = []
    for candidate_group in get_candidate_groups(list(genomes.values()), prioritize_additional_genomes):
        ordered_genomes.extend(sample_genomes(candidate_group, len(candidate_group)))

    if not ordered_genomes:
        return otu_genome_map

    otu_indices = np_rand.choice(len(unmatched_otus),len(unmatched_otus),replace=False)
    for i, genome in enumerate(ordered_genomes[:len(unmatched_otus)]):
        curr_otu = unmatched_otus[otu_indices[i]] #so we choose a random genome
        lineage, abundances = tax_profile[curr_otu]
        otu_genome_map[curr_otu] = (genome["tax_id"], genome["genome_id"], genome["path"], genome["original_path"], genome["source_file"], genome["source_line"], genome["novelty"], genome["label"], abundances)
#        if debug:
#            _log.warning("Filling up OTU %s (mapped tax id: %s) to genome with tax id %s" % (curr_otu, lin[0], tax_id))
    return otu_genome_map

"""
Take fasta input file and split by any N occurence (and remove Ns)
"""
def split_by_N(fasta_path, out_path, script):
    #os.system("scripts/split_fasta.pl %s %s" % (fasta_path, out_path))
    os.system(script+" %s %s" % (fasta_path, out_path))
    os.remove(fasta_path)

"""
Downloads the given genome and returns the out path
"""
def download_genome(genome, output_path, script):
    genome_path = os.path.dirname(output_path)
    if genome_path and not os.path.exists(genome_path):
        os.makedirs(genome_path)
    out_name = genome.rstrip().split('/')[-1]
    http_address = os.path.join(genome, out_name + "_genomic.fna.gz")
    opened = urllib.request.urlopen(http_address)
    out = output_path
    tmp_out = output_path + ".tmp"
    out_gz = out + ".gz"
    with open(out_gz,'wb') as outF:
        outF.write(opened.read())
    gf = gzip.open(out_gz)
    new_out = open(tmp_out,'wb')
    new_out.write(gf.read())
    gf.close()
    os.remove(out_gz)
    new_out.close()
    split_by_N(tmp_out, out, script)
    return out


def format_genome_transfer_error(error, otu, taxid, genome_id, genome_label, source_file,
                                 source_line, path, original_path, target_genome_path):
    source_entry = {
        "source_file": source_file,
        "source_line": source_line
    }
    return (
        "Could not copy/download selected genome after 10 tries.\n"
        "OTU: {otu}\n"
        "TaxID: {taxid}\n"
        "Genome ID: {genome_id}\n"
        "Genome label: {genome_label}\n"
        "Input entry: {input_entry}\n"
        "Original path: {original_path}\n"
        "Resolved path: {path}\n"
        "Target path: {target_path}\n"
        "Last error: {error}"
    ).format(
        otu=otu,
        taxid=taxid,
        genome_id=genome_id,
        genome_label=genome_label,
        input_entry=format_genome_record_source(source_entry),
        original_path=original_path,
        path=path,
        target_path=target_genome_path,
        error=repr(error)
    )


"""
Given the created maps and the old config files, creates the required files and new config
"""
def write_config(otu_genome_map, no_samples, script, out_path_genomes):
    out_path = "./"
    genome_to_id = os.path.join(out_path, "genome_to_id.tsv")
    #config.set('community0','id_to_genome_file', genome_to_id)
    metadata = os.path.join(out_path, "metadata.tsv")
    #config.set('community0','metadata',metadata)
    genome_filename_mapping = os.path.join(out_path, "genome_filename_mapping.tsv")
    #no_samples = int(config.get("Main","number_of_samples"))
    abundances = [os.path.join(out_path,"abundance_%s.tsv" % i) for i in range(no_samples)]
    #log.info("Downloading %s genomes" % len(otu_genome_map))
    
    #create_path = os.path.join(out_path_genomes,"genomes")
    create_path = out_path_genomes

    if not os.path.exists(create_path):
        os.makedirs(create_path)
    used_output_names = {}
    selected_genomes = []
    for otu in otu_genome_map:
        taxid, genome_id, path, original_path, source_file, source_line, novelty, genome_label, curr_abundances = otu_genome_map[otu]
        output_entry = {
            "source": "selected genome",
            "file": source_file,
            "line": source_line,
            "label": genome_label,
            "path": original_path
        }
        target_genome_path = build_output_genome_path(create_path, output_entry, used_output_names)
        selected_genomes.append(
            (
                otu,
                taxid,
                genome_id,
                path,
                original_path,
                source_file,
                source_line,
                novelty,
                genome_label,
                curr_abundances,
                target_genome_path
            )
        )

    abundance_streams = [open(abundance, 'w') for abundance in abundances]
    try:
        with open(genome_to_id, 'w') as gid, open(metadata,'w') as md, open(genome_filename_mapping, 'w') as mapping:
            md.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n") # write header
            mapping.write("original_path\tworkdir_path\n")
            for (otu, taxid, genome_id, path, original_path, source_file, source_line,
                 novelty, genome_label, curr_abundances, target_genome_path) in selected_genomes:
                counter = 0
                last_error = None
                while counter < 10:
                    try:
                        if is_remote_path(path):
                            genome_path = download_genome(path, target_genome_path, script)
                        else:
                            genome_path = target_genome_path
                            shutil.copy2(path, genome_path)
                        break
                    except Exception as e:
                        last_error = e
                        counter += 1
                        #log.error("Caught exception %s while moving/downloading genomes" % repr(e))
                        if counter >= 10:
                            raise ValueError(
                                format_genome_transfer_error(
                                    e,
                                    otu,
                                    taxid,
                                    genome_id,
                                    genome_label,
                                    source_file,
                                    source_line,
                                    path,
                                    original_path,
                                    target_genome_path
                                )
                            )
                if counter == 10:
    #            log.error("Genome %s (from %s, path %s) could not be downloaded after 10 tries, check your connection settings" % (otu, genome_id, path))
                    raise ValueError(
                        format_genome_transfer_error(
                            last_error,
                            otu,
                            taxid,
                            genome_id,
                            genome_label,
                            source_file,
                            source_line,
                            path,
                            original_path,
                            target_genome_path
                        )
                    )
                gid.write("%s\t%s\n" % (otu, genome_path))
                md.write("%s\t%s\t%s\t%s\n" % (otu,taxid,genome_id,novelty))
                mapping.write("%s\t%s\n" % (normalize_mapping_path(original_path), normalize_mapping_path(genome_path)))
                for abundance_stream, abundance_value in zip(abundance_streams, curr_abundances):
                    abundance_stream.write("%s\t%s\n" % (otu, abundance_value))
    finally:
        for abundance_stream in abundance_streams:
            abundance_stream.close()


def str2bool(x):
    return str(x).lower() in ("1", "true", "yes", "y")


def set_ranks(ncbi):
    global RANKS

    # Get all possible taxonomic ranks
    conn = sqlite3.connect(ncbi.dbfile)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT rank FROM species WHERE rank != ''")
    all_ranks = {row[0] for row in cursor.fetchall()}
    conn.close()

    if 'acellular root' in all_ranks:
        RANKS = NEW_RANKS
    else:
        RANKS = LEGACY_RANKS


def main(biom_profile, no_samples, reference_genomes, seed, mu, sigma,
         min_strains,  max_strains, debug, no_replace, fill_up, script,
         genomes_out_dir, additional_references, additional_genomes_quality_file, ncbi_taxdump_file,
         max_rank, prioritize_additional_genomes, additional_genomes_max_contamination,
         additional_genomes_min_completeness, additional_genomes_max_num_contigs):

    if additional_references == "None":
        additional_references = None
    if additional_genomes_quality_file == "None":
        additional_genomes_quality_file = None

    validate_strain_bounds(min_strains, max_strains)

    np_rand.seed(seed)

    genomes_map, total_genomes = read_genomes_list(reference_genomes, additional_references,
                                                   additional_genomes_quality_file,
                                                   additional_genomes_max_contamination,
                                                   additional_genomes_min_completeness,
                                                   additional_genomes_max_num_contigs)
    validate_local_genome_paths(genomes_map)

    ncbi = NCBITaxa(taxdump_file=ncbi_taxdump_file)
    set_ranks(ncbi)
    rank_cache = {}

    tax_profile = read_taxonomic_profile(biom_profile, no_samples)
    per_rank_map, genome_bucket_map = get_genomes_per_rank(genomes_map, ncbi, rank_cache)
    otu_genome_map, unmatched_otus = map_otus_to_genomes(tax_profile, per_rank_map, mu, sigma,
                                                         min_strains, max_strains, no_replace,
                                                         total_genomes, max_rank,
                                                         prioritize_additional_genomes, ncbi, rank_cache,
                                                         genome_bucket_map)

    if fill_up and len(unmatched_otus) > 0:
        otu_genome_map = fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile,
                                         debug, prioritize_additional_genomes)

    write_config(otu_genome_map, no_samples, script, genomes_out_dir)


if __name__ == "__main__":
    try:
        main(biom_profile = sys.argv[1],
            no_samples = int(sys.argv[2]),
            reference_genomes = sys.argv[3],
            seed = int(sys.argv[4]),
            mu = int(sys.argv[5]),
            sigma = int(sys.argv[6]),
            min_strains = int(sys.argv[7]),
            max_strains = int(sys.argv[8]),
            debug = str2bool(sys.argv[9]),
            no_replace = str2bool(sys.argv[10]),
            fill_up = str2bool(sys.argv[11]),
            script = sys.argv[12],
            genomes_out_dir = sys.argv[13],
            additional_references = sys.argv[14],
            additional_genomes_quality_file = sys.argv[15],
            ncbi_taxdump_file = sys.argv[16],
            max_rank = sys.argv[17],
            prioritize_additional_genomes = str2bool(sys.argv[18]),
            additional_genomes_max_contamination = float(sys.argv[19]),
            additional_genomes_min_completeness = float(sys.argv[20]),
            additional_genomes_max_num_contigs = int(float(sys.argv[21])))
    except ValueError as e:
        print("ERROR: {}".format(e), file=sys.stderr)
        sys.exit(1)
