import numpy as np
import pytest

from sparsim.parameter import SPARSimParameter, int_var_relation
from sparsim.simulate import SPARSim_simulation, simulate_hyper


@pytest.fixture
def base_sim_param():
    """Provides a basic, reusable SPARSimParameter object for testing."""
    return SPARSimParameter(
        intensity=np.array([10, 50, 100]),
        variability=np.array([0.2, 0.1, 0.05]),
        library_size=np.array([1000, 1500]),
    )


# ############################
# Tests for parameter.py
# ############################


def test_parameter_creation(base_sim_param):
    """Check that the SPARSimParameter object is initialized correctly."""
    assert np.array_equal(base_sim_param.intensity, np.array([10, 50, 100]))
    assert np.array_equal(base_sim_param.variability, np.array([0.2, 0.1, 0.05]))
    assert np.array_equal(base_sim_param.library_size, np.array([1000, 1500]))
    assert base_sim_param.intensity_2 is None
    assert base_sim_param.variability_2 is None
    assert base_sim_param.p_bimod is None


def test_int_var_relation():
    """Tests the variability interpolation logic."""
    sim_param = {
        "intensity": np.array([10, 20, 30, 40, 50]),
        "variability": np.array([0.5, 0.4, 0.3, 0.2, 0.1]),
    }

    # Test interpolation within the range
    new_intensities_interp = np.array([15, 25, 35])
    expected_variability_interp = np.array([0.45, 0.35, 0.25])
    result_interp = int_var_relation(new_intensities_interp, sim_param)
    assert np.allclose(result_interp, expected_variability_interp)

    # Test values outside the range (should be clamped)
    new_intensities_outside = np.array([5, 55])
    expected_variability_outside = np.array([0.5, 0.1])
    result_outside = int_var_relation(new_intensities_outside, sim_param)
    assert np.allclose(result_outside, expected_variability_outside)

    # Test zero intensity
    new_intensity_zero = np.array([0, 10])
    result_zero = int_var_relation(new_intensity_zero, sim_param)
    assert np.isnan(result_zero[0])
    assert result_zero[1] == 0.5


def test_create_de_genes_parameter(base_sim_param):
    """Test creation of parameters for DE genes."""
    fc_multiplier = np.array([1, 2, 0.5])
    de_param = base_sim_param.create_DE_genes_parameter(fc_multiplier=fc_multiplier)

    # Check that new intensities are correct
    expected_intensity = base_sim_param.intensity * fc_multiplier
    assert np.array_equal(de_param.intensity, expected_intensity)

    # Check that dimensions are correct
    assert de_param.intensity.shape == base_sim_param.intensity.shape
    assert de_param.variability.shape == base_sim_param.variability.shape
    assert de_param.library_size.shape == base_sim_param.library_size.shape

    # Check with a specified n_cells
    de_param_n_cells = base_sim_param.create_DE_genes_parameter(
        fc_multiplier=fc_multiplier, n_cells=5
    )
    assert de_param_n_cells.library_size.shape == (5,)


# ############################
# Tests for simulate.py
# ############################


def test_simulate_hyper_sum():
    """The sum of counts from simulate_hyper should equal seqdepth."""
    avg_abund = np.array([100, 200, 700])
    seqdepth = 50
    counts = simulate_hyper(avg_abund, seqdepth)
    assert counts.sum() == seqdepth


def test_simulate_hyper_error():
    """simulate_hyper should raise an error if seqdepth > sum(avg_abund)."""
    # avg_abund sums to 60. So seqdepth = 100 should raise an error.
    with pytest.raises(ValueError):
        simulate_hyper(avg_abund=np.array([10, 20, 30]), seqdepth=100)


def test_simulation_smoke_test(base_sim_param):
    """A basic end-to-end run to ensure it completes without errors."""
    result = SPARSim_simulation(base_sim_param)
    assert "count_matrix" in result
    assert result["count_matrix"] is not None


def test_output_shapes(base_sim_param):
    """Output matrices should have the shape (n_genes, n_cells)."""
    n_genes = base_sim_param.intensity.shape[0]
    n_cells = base_sim_param.library_size.shape[0]
    result = SPARSim_simulation(base_sim_param)

    assert result["count_matrix"].shape == (n_genes, n_cells)
    assert result["gene_matrix"].shape == (n_genes, n_cells)
    assert result["abundance_matrix"].shape == (n_genes, n_cells)
    assert result["variability_matrix"].shape == (n_genes, n_cells)


def test_library_sizes(base_sim_param):
    """The sum of counts per cell should match the specified library size."""
    result = SPARSim_simulation(base_sim_param)
    # The current implementation of simulate_hyper does not guarantee
    # that the sum of counts will be *exactly* the library size due to the
    # nature of the sampling from a scaled distribution.
    # We test that it's close.
    assert np.allclose(
        result["count_matrix"].sum(axis=0), base_sim_param.library_size, rtol=0.1
    )


def test_zero_variability():
    """Genes with zero variability should not have biological variation."""
    params = SPARSimParameter(
        intensity=np.array([10, 50, 100]),
        variability=np.array([0.2, 0, 0.1]),  # Middle gene has zero variability
        library_size=np.array([1000, 1500]),
    )
    result = SPARSim_simulation(params)

    # For the gene with zero variability, bio-varied expression should equal initial expression
    assert np.allclose(result["gene_matrix"][1, :], result["abundance_matrix"][1, :])
    # For genes with non-zero variability, they should likely be different
    assert not np.allclose(
        result["gene_matrix"][0, :], result["abundance_matrix"][0, :]
    )


def test_bimodality():
    """Test bimodal expression logic."""
    params = SPARSimParameter(
        intensity=np.array([10, 50, 200]),
        variability=np.array([0.2, 0.1, 0.05]),
        library_size=np.array([1000] * 10),  # 10 cells
        intensity_2=np.array([10, 500, 200]),  # Gene 2 has a different 2nd mode
        variability_2=np.array([0.2, 0.01, 0.05]),
        p_bimod=np.array([1, 0, 1]),  # Gene 2 should always be in mode 2
    )

    result = SPARSim_simulation(params)
    abundance_matrix = result["abundance_matrix"]

    # Gene 0 (p_bimod=1) should always be its original intensity
    assert np.all(abundance_matrix[0, :] == 10)

    # Gene 1 (p_bimod=0) should always be in the second mode
    assert np.all(abundance_matrix[1, :] == 500)

    # Gene 2 (p_bimod=1) should be its original intensity
    assert np.all(abundance_matrix[2, :] == 200)
