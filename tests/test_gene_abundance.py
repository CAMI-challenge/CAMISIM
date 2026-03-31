"""Tests for pipelines/metatranscriptomic/scripts/get_gene_abundance.py."""

import random
import pytest

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metatranscriptomic", "scripts"))
from get_gene_abundance import GeneAbundace


@pytest.fixture
def ga():
    return GeneAbundace()


class TestLtZero:
    def test_positive(self, ga):
        assert ga.lt_zero(5.0) == 5.0

    def test_zero(self, ga):
        assert ga.lt_zero(0) == 0.001

    def test_negative(self, ga):
        assert ga.lt_zero(-1.0) == 0.001


class TestAddInitialLogDistribution:
    def test_fills_specified_sample(self, ga):
        abundances = {"g1": [0.0, 0.0], "g2": [0.0, 0.0]}
        random.seed(42)
        ga.add_initial_log_distribution(abundances, 1.0, 1.0, 0)
        for val in abundances.values():
            assert val[0] > 0
            assert val[1] == 0.0

    def test_fills_second_sample(self, ga):
        abundances = {"g1": [0.0, 0.0]}
        random.seed(42)
        ga.add_initial_log_distribution(abundances, 1.0, 1.0, 1)
        assert abundances["g1"][0] == 0.0
        assert abundances["g1"][1] > 0

    def test_reproducible(self, ga):
        a1 = {"g1": [0.0], "g2": [0.0]}
        random.seed(99)
        ga.add_initial_log_distribution(a1, 2.0, 0.5, 0)

        a2 = {"g1": [0.0], "g2": [0.0]}
        random.seed(99)
        ga.add_initial_log_distribution(a2, 2.0, 0.5, 0)
        assert a1 == a2


class TestAddTimeseriesGauss:
    def test_extinction_propagates(self, ga):
        abundances = {"g1": [0.0, 0.0, 0.0]}
        ga.add_timeseries_gauss(abundances, 0, 1.0)
        assert abundances["g1"][1] == 0.0
        assert abundances["g1"][2] == 0.0

    def test_positive_chain(self, ga):
        abundances = {"g1": [5.0, 0.0, 0.0]}
        random.seed(42)
        ga.add_timeseries_gauss(abundances, 0, 0.001)
        assert all(v > 0 for v in abundances["g1"])


class TestAddGenewiseLogDistribution:
    def test_fills_from_initial(self, ga):
        abundances = {"g1": [2.0, 0.0, 0.0]}
        random.seed(42)
        ga.add_genewise_log_distribution(abundances, 0.5)
        assert all(v > 0 for v in abundances["g1"])


class TestRandomDistributionToRelativeAbundance:
    def test_sums_to_one(self, ga):
        abundances = {"g1": [3.0, 6.0], "g2": [7.0, 4.0]}
        ga.random_distribution_to_relative_abundance(abundances, 2)
        for s in range(2):
            total = sum(v[s] for v in abundances.values())
            assert abs(total - 1.0) < 1e-9

    def test_proportions(self, ga):
        abundances = {"g1": [2.0], "g2": [8.0]}
        ga.random_distribution_to_relative_abundance(abundances, 1)
        assert abs(abundances["g1"][0] - 0.2) < 1e-9
        assert abs(abundances["g2"][0] - 0.8) < 1e-9

    def test_single_gene(self, ga):
        abundances = {"g1": [5.0, 3.0]}
        ga.random_distribution_to_relative_abundance(abundances, 2)
        assert abundances["g1"][0] == 1.0
        assert abundances["g1"][1] == 1.0
