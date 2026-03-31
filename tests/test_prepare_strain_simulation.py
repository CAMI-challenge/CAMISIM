"""Tests for pipelines/metagenomic/scripts/prepare_strain_simulation.py — strain distribution."""

import random
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metagenomic", "scripts"))
from prepare_strain_simulation import PrepareStrainSimulation


@pytest.fixture
def pss():
    return PrepareStrainSimulation()


class TestInit:
    def test_init_without_seed(self, pss):
        pss.init(seed=None)  # should not raise

    def test_init_with_seed(self, pss):
        pss.init(seed=42)  # should not raise

    def test_seed_reproducibility(self, pss):
        pss.init(seed=123)
        val1 = random.random()
        pss.init(seed=123)
        val2 = random.random()
        assert val1 == val2


class TestGetGenomeAmountsGeometric:
    def test_sum_equals_max(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric(0.5, 20)
        assert sum(amounts) == 20

    def test_all_positive(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric(0.5, 10)
        assert all(a > 0 for a in amounts)

    def test_high_probability_many_singletons(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric(0.99, 100)
        singletons = sum(1 for a in amounts if a == 1)
        assert singletons > len(amounts) * 0.5

    def test_small_max(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric(0.5, 1)
        assert sum(amounts) == 1


class TestGetGenomeAmountsGeometricFix:
    def test_sum_equals_max(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric_fix(5, 20)
        assert sum(amounts) == 20

    def test_correct_number_of_entries(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric_fix(5, 20)
        assert len(amounts) == 5

    def test_all_at_least_one(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric_fix(3, 10)
        assert all(a >= 1 for a in amounts)

    def test_equal_when_max_equals_real(self, pss):
        pss.init(seed=42)
        amounts = pss.get_genome_amounts_geometric_fix(5, 5)
        assert amounts == [1, 1, 1, 1, 1]
