import numpy as np
from scipy.interpolate import interp1d


def int_var_relation(intensity_values: np.ndarray, sim_param: dict) -> np.ndarray:
    """
    Associate gene intensity and gene variability values

    This function associates variability values to a set of intensity values (obtained,
    for example, during DE genes simulation) `intensity_values`.
    The association between gene intensity values and gene variability values is learned
    from the one described in an existing simulation parameter `sim_param`,
    using a linear interpolation.

    Parameters
    ----------
    intensity_values : numpy.ndarray
        The (new) intensity values to which we want to associate variability values.
    sim_param : dict
        Simulation parameter from which learn gene intensity vs gene var association.

    Returns
    -------
    numpy.ndarray
        A vector size of `intensity_values` containing the computed var values.
    """
    not_na_var_ind = ~np.isnan(sim_param["variability"])
    int_ord = np.sort(sim_param["intensity"][not_na_var_ind])
    var_not_na = sim_param["variability"][not_na_var_ind]
    var_ord = var_not_na[np.argsort(sim_param["intensity"][not_na_var_ind])]

    interp_variability = np.empty_like(intensity_values, dtype=float)
    interp_variability.fill(np.nan)

    # handle extreme variability values
    interp_variability[intensity_values >= np.max(int_ord)] = sim_param["variability"][
        np.argmax(int_ord)
    ]
    interp_variability[intensity_values <= np.min(int_ord)] = sim_param["variability"][
        np.argmin(int_ord)
    ]

    # do interpolation
    to_interp = np.where(np.isnan(interp_variability))[0]

    if len(to_interp) > 0:
        f = interp1d(int_ord, var_ord)
        interp_variability[to_interp] = f(intensity_values[to_interp])

    interp_variability[intensity_values == 0] = np.nan

    return interp_variability


class SPARSimParameter:
    """
    SPARSim simulation parameter class.
    """

    def __init__(
        self,
        intensity: np.ndarray,
        variability: np.ndarray,
        library_size: np.ndarray,
        intensity_2: None | np.ndarray = None,
        variability_2: None | np.ndarray = None,
        p_bimod: None | np.ndarray = None,
    ):
        """
        Create a SPARSim simulation parameter.

        Parameters
        ----------
        intensity : numpy.ndarray
            Array of gene expression intensity values.
        variability : numpy.ndarray
            Array of gene expression variability values.
        library_size : numpy.ndarray
            Array of library size values.
        intensity_2 : numpy.ndarray, optional
            Array of gene expression intensity values for the second expression mode.
        variability_2 : numpy.ndarray, optional
            Array of gene expression variability values for the second expression mode.
        p_bimod : numpy.ndarray, optional
            Array of bimodal gene expression probabilities.
        """
        self.intensity = intensity
        self.variability = variability
        self.library_size = library_size

        if (
            intensity_2 is not None
            and variability_2 is not None
            and p_bimod is not None
        ):
            self.intensity_2 = intensity_2
            self.variability_2 = variability_2
            self.p_bimod = p_bimod
        else:
            self.intensity_2 = None
            self.variability_2 = None
            self.p_bimod = None

    def create_DE_genes_parameter(
        self,
        fc_multiplier: np.ndarray,
        n_cells: int | None = None,
        lib_size_DE=None,
    ):
        """
        Create a SPARSim simulation parameter with DE genes.

        Parameters
        ----------
        fc_multiplier : numpy.ndarray
            Array of fold change multipliers.
        n_cells : int, optional
            Number of cells to simulate.
        lib_size_DE : numpy.ndarray, optional
            Array of library size values.

        Returns
        -------
        SPARSimParameter
            A new SPARSimParameter object with DE genes.
        """
        if lib_size_DE is None:
            if n_cells is None:
                lib_size_DE = self.library_size
            else:
                lib_size_DE = np.random.choice(
                    self.library_size, size=n_cells, replace=True
                )

        sim_param_DE_intensity = fc_multiplier * self.intensity

        sim_param_DE_variability = int_var_relation(
            sim_param_DE_intensity,
            {
                "intensity": self.intensity,
                "variability": self.variability,
            },
        )

        return SPARSimParameter(
            intensity=sim_param_DE_intensity,
            variability=sim_param_DE_variability,
            library_size=lib_size_DE,
        )
