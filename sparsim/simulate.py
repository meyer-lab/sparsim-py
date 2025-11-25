import numpy as np

from .parameter import SPARSimParameter


def simulate_hyper(avg_abund: np.ndarray, seqdepth: int, digits=None) -> np.ndarray:
    """
    Simulate technical variability (multivariate hypergeometric) using NumPy.

    Args:
        avg_abund (np.array): Intensity values for each feature (for a single cell).
        seqdepth (int): Sequencing depth (sample size).
        digits (int, optional): Used for sanity check against seqdepth.
                                Defaults to None (skips check).

    Returns:
        np.array: Array of count values matching the length of avg_abund.
    """
    # max_val should be the sum of avg_abund to represent the total pool size
    max_val = int(np.sum(avg_abund))

    # 1. Sanity Check
    # In Python, we check if seqdepth exceeds the population size (max_val)
    # because we cannot sample unique integers if sample_size > population.
    if seqdepth > max_val:
        raise ValueError(
            f"seqdepth ({seqdepth}) cannot be greater than max_val (sum of avg_abund: {max_val})"
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
    # We slice the result to discard these 'missed' hits.
    return counts[: len(avg_abund)]


def SPARSim_simulation(
    params: SPARSimParameter, rng: np.random.Generator | None | int = None
):
    """
    Simulate a raw count table.
    """
    if rng is None:
        rng = np.random.default_rng()
    elif rng is int:
        rng = np.random.default_rng(rng)

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
        # Identify genes having bimodal parameters defined and are eligible for mode 2
        # This is a boolean array of length num_genes
        eligible_for_mode2 = (~np.isnan(params.intensity_2)) & (
            (1 - params.p_bimod) > 0
        )

        for cell in range(params.library_size.size):
            # For the eligible genes, determine which switch to mode 2 for this cell
            # Results in bool arr of length len(eligible_for_mode2[eligible_for_mode2])
            # which is the number of eligible genes.
            switch_to_mode2_subset = rng.binomial(
                1, 1 - params.p_bimod[eligible_for_mode2]
            ).astype(bool)

            # Now, we need to create a boolean index for the full gene_expression_matrix
            # Only apply the switch to mode 2 to those genes that were eligible.
            actual_switch_indices = np.zeros(params.intensity.shape, dtype=bool)
            actual_switch_indices[eligible_for_mode2] = switch_to_mode2_subset

            # Apply the update using the full-length boolean index
            gene_expression_matrix[actual_switch_indices, cell] = params.intensity_2[
                actual_switch_indices
            ]
            gene_expression_var_matrix[actual_switch_indices, cell] = (
                params.variability_2[actual_switch_indices]
            )

    # Generate biological variability
    gene_expression_matrix_bio_var = gene_expression_matrix.copy()

    with np.errstate(divide="ignore", invalid="ignore"):
        shape = 1 / gene_expression_var_matrix
        scale = gene_expression_var_matrix * gene_expression_matrix

    gene_expression_matrix_bio_var = rng.gamma(shape=shape, scale=scale)

    zero_var_index = gene_expression_var_matrix == 0
    gene_expression_matrix_bio_var[zero_var_index] = gene_expression_matrix[
        zero_var_index
    ]

    max_lib_size = params.library_size.max()
    input_fragment_lib_size = max_lib_size * 100
    digits = np.floor(np.log10(input_fragment_lib_size)) + 1
    new_libsize = 10**digits

    cell_sums = gene_expression_matrix_bio_var.sum(axis=0)
    cell_sums_reshaped = cell_sums[
        np.newaxis, :
    ]  # Reshape to (1, n_cells) for broadcasting with (n_genes, n_cells)

    gene_expression_matrix_bio_var_scaled = np.round(
        (gene_expression_matrix_bio_var * new_libsize)  # (n_genes, n_cells) * scalar
        / cell_sums_reshaped  # (1, n_cells)
    )

    print("Simulating technical variability ...")
    sim_count_matrix = np.zeros_like(gene_expression_matrix_bio_var_scaled)

    for cell_idx, sample_lib_size in enumerate(params.library_size):
        sim_count_matrix[:, cell_idx] = simulate_hyper(
            avg_abund=gene_expression_matrix_bio_var_scaled[:, cell_idx],
            seqdepth=sample_lib_size,
            digits=digits,
        )

    return {
        "count_matrix": sim_count_matrix,
        "gene_matrix": gene_expression_matrix_bio_var,
        "abundance_matrix": gene_expression_matrix,
        "variability_matrix": gene_expression_var_matrix,
    }
