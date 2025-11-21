#' Create SPARSim simulation parameter
#'
#' Function to create a SPARSim simulation parameter.
#'
#' To simulate N feature (e.g. genes), user must specify N values of gene expression level and gene expression variability in the function input parameters \code{intensity} and \code{variability}, respectively.
#' To simulate M samples (i.e. cells), user must specify M values of sample library size in the function input parameter \code{library_size}.
#'
#' User can optionally specify the names to assign at the single feature and sample to simulate (function input parameters \code{feature_names} and \code{sample_names}, respectively,
#' as well as the name of the experimental condition (function input parameter \code{condition_name}). If the user does not specify such information, the function will set some default values.
#'
#' To simulate T different experimental conditions in a single count table, then T different simulation parameters must be created.
#'
#' @param intensity Array of gene expression intensity values
#' @param variability Array of gene expression variability values
#' @param library_size Array of library size values
#' @param feature_names Array of feature names. It must be of the same length of \code{intensity} array. If NA (default), feature will be automatically named "gene_1", "gene_2", ... "gene_<N>", where N = length(intensity)
#' @param sample_names Array of sample names. It must be of the same length of \code{library_size} array. If NA (defatul), sample will be automatically named "<condition_name>_cell1", "<condition_name>_cell2", ..., "<condition_name>_cell<M>", where M = length(library_size)
#' @param condition_name Name associated to the current experimental condition. If NA (default), it will be set to "cond<l1><l2>", where l1 and l2 are two random letters.
#' @param intensity_2 Array of gene expression intensity values for the second expression mode, if simulating genes with bimodal gene expression. Entries containing \code{NAs} will be ignored. If NULL (default), no bimodal gene expression is simulated.
#' @param variability_2 Array of gene expression variability values for the second expression mode, if simulating genes with bimodal gene expression. If NULL (default), no bimodal gene expression is simulated.
#' @param p_bimod Array of bimodal gene expression probabilities; the i-th value indicates the probability \code{p} of the i-th gene to be expressed in the first mode
#' (i.e. the one specified in the i-th entries of parameters \code{intensity} and \code{variability}); with probability \code{1-p} the i-th gene will be expressed in the second mode (i.e. the one specified in the i-th entries of parameters \code{intensity_2} and \code{variability_2})
#' @return SPARSim simulation parameter describing one experimental condition
#' @export
SPARSim_create_simulation_parameter <- function(intensity, variability, library_size, feature_names = NA, sample_names = NA, condition_name = NA,
                                                intensity_2 = NULL, variability_2 = NULL, p_bimod = NULL){

  # assign feature name to "intensity" and "variability"
  if(length(feature_names)==1){ # potential NA
    if(is.na(feature_names)){
      feature_names <- paste0("gene_", c(1:length(intensity)) )
    }
  }
  names(intensity) <- feature_names
  names(variability) <- feature_names

  # assign condition name
  if(is.na(condition_name)){
    condition_name <- paste0("cond_", paste0(sample(LETTERS, size = 2), collapse = "") )
  }

  # assign sample names to "library_size"
  if(length(sample_names)==1){
    if(is.na(sample_names)){
      sample_names <- paste0(condition_name,"_cell",c(1:length(library_size)))
    }
  }
  names(library_size) <- sample_names

  cond_param <- list()
  cond_param$intensity <- intensity
  cond_param$variability <- variability
  cond_param$lib_size <- library_size
  cond_param$name <- condition_name

  # if bimodal gene expression is required, then set the related fields
  if(!is.null(intensity_2) & !is.null(variability_2) & !is.null(p_bimod)){
    names(intensity_2) <- feature_names
    names(variability_2) <- feature_names
    cond_param$intensity_2 <- intensity_2
    cond_param$variability_2 <- variability_2
    cond_param$p_bimod <- p_bimod
  }

  return(cond_param)
}




#' Create SPARSim simulation parameter with DE genes
#'
#' Function to create a SPARSim simulation parameter with DE genes respect to another SPARSim simulation parameter provided as input.
#'
#' To simulate N feature (e.g. genes), user must specify N values of gene expression level and gene expression variability in the function input parameters \code{intensity} and \code{variability}, respectively.
#' To simulate M samples (i.e. cells), user must specify M values of sample library size in the function input parameter \code{library_size}.
#'
#' User can optionally specify the names to assign at the single feature and sample to simulate (function input parameters \code{feature_names} and \code{sample_names}, respectively,
#' as well as the name of the experimental condition (function input parameter \code{condition_name}). If the user does not specify such information, the function will set some default values.
#'
#' To simulate T different experimental conditions in a single count table, then T different simulation parameters must be created.
#'
#' @param sim_param SPARSim simulation parameter used as prototype
#' @param fc_multiplier Array of fold change multipliers
#' @param N_cells Number of cells to simulate.
#' If NULL (default), the number of cells will be set to the length of \code{lib_size_DE} or to the number of cells in \code{sim_param}
#' @param lib_size_DE Array of library size values. It must be of the same length of \code{sample_names} and equal to \code{N_cells}, if they are specified.
#' If NULL (default), and \code{N_cells} is NULL, the library size values from sim_param will be used.
#' If NULL (default), and \code{N_cells} is specified, then N_cell randomly chosen library size values from sim_param will be used.
#' @param sample_names Array of sample names. It must be of the same length of \code{lib_size_DE} and equal to \code{N_cells}, if they are specified.
#' If NULL (defatul), sample will be automatically named "<condition_name>_cell1", "<condition_name>_cell2", ...
#' @param condition_name Name associated to the experimental condition.
#' If NULL (default), the name is set to the concatenation of the name in \code{sim_param} with the suffix \code{"_DE"}
#'
#' @return SPARSim simulation parameter with DE genes respect to the input SPARSim simulation parameter \code{sim_param}
#'
#' @export
SPARSim_create_DE_genes_parameter <- function(sim_param, fc_multiplier, N_cells = NULL, lib_size_DE = NULL, sample_names = NULL, condition_name = NULL){

  if(!is.null(lib_size_DE)|!is.null(sample_names)|!is.null(N_cells)){ # at least one of N_cells, lib_size_DE and sample_names was specified by user
    if(!is.null(lib_size_DE)&!is.null(sample_names)&!is.null(N_cells)){
      if( (N_cells!=length(lib_size_DE)) || (N_cells!=length(sample_names)) ){
        stop("The number of cells is not consistent across N_cells, lib_size_DE and sample_names")
      }
    }
  }


  if(is.null(lib_size_DE)){ # if no library size values are specified,
    if(is.null(N_cells)){ # if no N_cell is specified
      lib_size_DE <- sim_param$lib_size # then use the library size values from sim_param
    }else{
      lib_size_DE <- sample(sim_param$lib_size, size = N_cells, replace = TRUE) # else sample N_cells values from library size values in sim_param
    }
  }

  # if no name is specified, then concatenate the name in sim_param with the suffix "_DE"
  if(is.null(condition_name)){
    condition_name <- paste0(sim_param_DE$name, "_DE")
  }


  if(is.null(sample_names)){# if no sample names are specified, then concatenate the condition name with the suffices "_cell1", "_cell2", ...
    sample_names <- paste0(condition_name,"_cell",c(1:length(lib_size_DE)))
  }else{ # if sample names are specified, then check that they are consistent with the number of cells
    if(length(sample_names)!=lib_size_DE){
      stop("The number of samples names is not compatible with the number of cells to simulate")
    }
  }

  sim_param_DE_intensity <- fc_multiplier * sim_param$intensity
  sim_param_DE <- SPARSim_create_simulation_parameter(intensity = sim_param_DE_intensity,
                                                      variability = int_var_relation(sim_param_DE_intensity, sim_param),
                                                      library_size = lib_size_DE,
                                                      feature_names = names(sim_param$intensity),
                                                      sample_names = sample_names,
                                                      condition_name = condition_name)

  return(sim_param_DE)
}
