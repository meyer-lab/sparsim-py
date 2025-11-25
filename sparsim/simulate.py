import numpy as np

from .parameter import SPARSimParameter


def simulate_hyper(avg_abund, seqdepth, max_val, digits=None):
    """
    Simulate technical variability (multivariate hypergeometric) using NumPy.

    Args:
        avg_abund (np.array): Array containing intensity values for each feature.
        seqdepth (int): Sequencing depth (sample size).
        max_val (int): Max value for random number generation space.
        digits (int, optional): Used for sanity check against seqdepth.
                                Defaults to None (skips check).

    Returns:
        np.array: Array of count values matching the length of avg_abund.
    """

    # 1. Sanity Check
    # In Python, we check if seqdepth exceeds the population size (max_val)
    # because we cannot sample unique integers if sample_size > population.
    if seqdepth > max_val:
        raise ValueError(
            f"seqdepth ({seqdepth}) cannot be greater than max_val ({max_val})"
        )

    if digits is not None and seqdepth > 10**digits:
        print("seqdepth must be <= than 10^digits")
        return None

    # 2. Generate Random Indices (The "While" Loop replacement)
    # The R code generates 'seqdepth' unique integers between 0 and max_val.
    # NumPy does this in a single highly-optimized C-level call.
    rng = np.random.default_rng()

    # unique indices sampled from range [0, max_val)
    random_indices = rng.choice(max_val, size=seqdepth, replace=False)

    # 3. Map Indices to Features (The "Cumulative" Step)
    cumulative_orig = np.cumsum(avg_abund)

    # R's findInterval(..., left.open=T) is equivalent to searchsorted(..., side='right')
    # This finds which 'bin' (gene/feature) each random index falls into.
    bin_indices = np.searchsorted(cumulative_orig, random_indices, side="right")

    # 4. Count Occurrences (The "Tabulate" Step)
    # Count how many times each bin index appears.
    # minlength ensures we get a count for every gene, even if 0.
    counts = np.bincount(bin_indices, minlength=len(avg_abund))

    # 5. Handle "Overflow"
    # If max_val > sum(avg_abund), some random indices might fall beyond the last gene.
    # searchsorted assigns these to index == len(avg_abund).
    # We slice the result to discard these 'missed' hits, strictly matching the R output shape.
    return counts[: len(avg_abund)]


def SPARSim_simulation(
    params: SPARSimParameter,
):
    """
    Simulate a raw count table.
    """

    # Initialize gene expression and variability matrices
    # genes x cells
    gene_expression_matrix = np.tile(
        params.intensity[:, np.newaxis],
        len(params.library_size),
    )

    variability_values = params.variability.copy()
    variability_values[np.isnan(variability_values)] = 0
    variability_values[variability_values < 0] = np.min(
        variability_values[variability_values > 0]
    )
    gene_expression_var_matrix = np.tile(
        variability_values[:, np.newaxis], len(params.library_size)
    )

    # Sample cells with alternative expression
    if params.p_bimod is not None:
        bimodal_genes_id = (~np.isnan(params.intensity_2)) & ((1 - params.p_bimod) > 0)

        for cell in range(params.library_size.size):
            mode2_ind = np.random.binomial(
                1, 1 - params.p_bimod[bimodal_genes_id]
            ).astype(bool)
            gene_expression_matrix[mode2_ind, cell] = params.intensity_2[mode2_ind]
            gene_expression_var_matrix[mode2_ind, cell] = params.variability_2[
                mode2_ind
            ]

    # Generate biological variability
    gene_expression_matrix_bio_var = gene_expression_matrix.copy()

    with np.errstate(divide="ignore", invalid="ignore"):
        shape = 1 / gene_expression_var_matrix
        scale = gene_expression_var_matrix * gene_expression_matrix

    gene_expression_matrix_bio_var = np.random.gamma(shape=shape, scale=scale)

    zero_var_index = gene_expression_var_matrix == 0
    gene_expression_matrix_bio_var[zero_var_index] = gene_expression_matrix[
        zero_var_index
    ]

    max_lib_size = params.library_size.max()
    input_fragment_lib_size = max_lib_size * 100
    digits = np.floor(np.log10(input_fragment_lib_size)) + 1
    new_libsize = 10**digits

    gene_expression_matrix_bio_var_scaled = np.round(
        (gene_expression_matrix_bio_var.T * new_libsize)
        / gene_expression_matrix_bio_var.sum(axis=0)
    ).T
    num_fragment = gene_expression_matrix_bio_var_scaled.sum(axis=0)

    print("Simulating technical variability ...")
    sim_count_matrix = np.zeros_like(gene_expression_matrix_bio_var_scaled)

    for cell_idx, sample_lib_size in enumerate(params.library_size):
        sim_count_matrix[:, cell_idx] = simulate_hyper(
            avg_abund=gene_expression_matrix_bio_var_scaled,
            seqdepth=sample_lib_size,
            max_val=int(num_fragment),
            digits=digits,
        )

    return {
        "count_matrix": sim_count_matrix,
        "gene_matrix": gene_expression_matrix_bio_var,
        "abundance_matrix": gene_expression_matrix,
        "variability_matrix": gene_expression_var_matrix,
    }
