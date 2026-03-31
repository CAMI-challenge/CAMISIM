"""Tests for convert_config.py — v1 INI to v2 Nextflow config conversion."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from convert_config import convert_to_nextflow_config


class TestConvertToNextflowConfig:
    def test_output_directory_renamed(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\noutput_directory = /my/output\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        content = out.read_text()
        assert 'outdir = "/my/output"' in content

    def test_gsa_true(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\ngsa = True\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "gsa = true" in out.read_text()

    def test_gsa_false(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\ngsa = False\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "gsa = false" in out.read_text()

    def test_anonymous_renamed_to_anonymization(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nanonymous = True\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "anonymization = true" in out.read_text()

    def test_type_quoted(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\ntype = nanosim3\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert 'type = "nanosim3"' in out.read_text()

    def test_skipped_keys(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nphase = 0\nmax_processors = 8\ncompress = True\nsamtools = /usr/bin/samtools\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        content = out.read_text()
        assert "phase" not in content.split("biom_profile")[0]  # before defaults
        assert "max_processors" not in content.split("biom_profile")[0]
        assert "compress" not in content.split("biom_profile")[0]

    def test_metadata_renamed(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nmetadata = /path/to/meta.tsv\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert 'metadata_file = "/path/to/meta.tsv"' in out.read_text()

    def test_id_to_genome_file_renamed(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nid_to_genome_file = /path/to/genomes.tsv\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert 'genome_locations_file = "/path/to/genomes.tsv"' in out.read_text()

    def test_error_profiles_renamed(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nerror_profiles = /path/to/profiles\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert 'base_profile_name = "/path/to/profiles"' in out.read_text()

    def test_fragment_size_sd_renamed(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nfragment_size_standard_deviation = 27\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "fragment_size_sd = 27" in out.read_text()

    def test_pooled_gsa_true(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\npooled_gsa = True\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "pooled_gsa = true" in out.read_text()

    def test_pooled_gsa_false(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\npooled_gsa = False\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert 'pooled_gsa = "[]"' in out.read_text()

    def test_view_renamed_to_verbose(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nview = True\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "verbose = true" in out.read_text()

    def test_defaults_appended(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nsize = 0.1\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        content = out.read_text()
        assert 'biom_profile=""' in content
        assert "conda.enabled = true" in content

    def test_passthrough_key(self, tmp_path):
        ini = tmp_path / "test.ini"
        ini.write_text("[Main]\nnumber_of_samples = 5\n")
        out = tmp_path / "out.config"
        convert_to_nextflow_config(str(ini), str(out))
        assert "number_of_samples = 5" in out.read_text()
