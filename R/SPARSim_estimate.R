
#' Estimate SPARSIm "intensity" parameter
#'
#' Function to estimate the intensity values from the genes in \code{data}. The intensity is computed as mean of normalized counts for each gene.
#'
#' This function is used in \code{SPARSim_estimate_parameter_from_data} to compute SPARSim "intensity" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data normalized count data matrix (gene on rows, samples on columns). \code{rownames(data)} must contain gene names.
#' @return An array of intensity values having \code{N_genes} elements (\code{N_genes = nrow(data)}). Array entries are named with gene names.
#' @export
SPARSim_estimate_intensity <- function(data){

  N_genes <- nrow(data)
  intensity_val <- array(0, dim = N_genes)
  names(intensity_val) <- rownames(data)

  # intesity as mean of normalized counts
  intensity_val <- apply(data, 1, mean)
  rm(data); gc()
  return (intensity_val)

}


#' Estimate SPARSim "variability" parameter
#'
#' Function to estimate the variability values from the genes in \code{data}.
#'
#' This function is used in \code{SPARSim_estimate_parameter_from_data} to compute SPARSim "variability" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data raw count data matrix (gene on rows, samples on columns)
#' @return An array of variability values having \code{N_genes} elements (\code{N_genes = nrow(data)})
#' @export
SPARSim_estimate_variability <- function (data){

  # compute percentage of entry >=0 for each gene
  perc_not_zeros_per_gene <- rowSums(data>0)/ncol(data)

  perc_not_zeros = 0
  ind_pass_filter <- (perc_not_zeros_per_gene > perc_not_zeros)
  N_gene_pass_filter <- sum(ind_pass_filter)

  variability <- array(NA, dim = nrow(data))

  # variability as edgeR dispersion of RAW counts in data
  if(N_gene_pass_filter > 1){

    #library(edgeR)
    f <- edgeR::DGEList(as.matrix(data[ind_pass_filter,]))
    gc()

    f <- edgeR::estimateGLMCommonDisp(f)
    gc()

    f <- edgeR::estimateGLMTrendedDisp(f)
    gc()

    if(any(is.na(f$trended.dispersion))){
      f$trended.dispersion <- NULL
      print("Trended dispresion is NA")
    }

    f <- edgeR::estimateGLMTagwiseDisp(f)
    gc()

    variability[ind_pass_filter] <- f$tagwise.dispersion
  }
  rm(data); gc()
  return(variability)

}



#' Estimate SPARSim "library size" parameter
#'
#' Function to estimate the library sizes from the samples in \code{data}.
#'
#' This function is used in \code{SPARSim_estimate_parameter_from_data} to compute SPARSim "library size" parameter, given a real count table as input.
#' If the count table contains more than one experimental condition, then the function is applied to each experimental conditions.
#'
#' @param data raw count data matrix (gene on rows, samples on columns)
#' @return An array of library size values having \code{N_samples} elements (\code{N_samples = ncol(data)})
#' @export
SPARSim_estimate_library_size <- function (data){
  return(colSums(data))
}



#' Estimate SPARSim simulation parameter from a given count table
#'
#' Function to estimate SPARSim simulation parameters (intensity, variability and library sizes) from a real count table.
#' If the real count table contains more than one experimental condition, it is possible to estimate the parameters for each experimental condition.
#'
#' @param raw_data count matrix (gene on rows, samples on columns) containing raw count data
#' @param norm_data count matrix (gene on rows, samples on columns) containing normalized count data
#' @param conditions list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @return A SPARSim simulation parameters
#' @export
SPARSim_estimate_parameter_from_data <- function (raw_data, norm_data, conditions){

  # Number of experimental conditions
  N_cond <- length(conditions)

  # Number of genes
  N_feature <- nrow(raw_data)


  dataset_parameter <- list()

  cond_i <- 1
  for (cond in conditions){# for each experimental condition

    print(paste0("Experimental condition ", cond_i))

    gene_names <- rownames(raw_data)
    sample_names <- colnames(raw_data)

    print("...estimating gene intensity")
    # estimate intesity
    #estim_intensity <- estimate_intensity (data = norm_data[,cond])
    N_genes <- length(gene_names)
    estim_intensity <- array(0, dim = N_genes)
    names(estim_intensity) <- gene_names
    for(i in 1:N_genes){
      estim_intensity[i] <- mean(norm_data[i,cond])
    }
    #estim_intensity <- apply(norm_data[,cond], 1, mean)# intesity as mean of normalized counts
    gc()

    print("...estimating gene variability")
    # estimate variability
    estim_variability <- SPARSim_estimate_variability (data = raw_data[,cond])

    print("...estimating library size")
    # estimate library size
    estim_lib_size <- SPARSim_estimate_library_size(data = raw_data[,cond])

    # get feature names, if present
    estim_feature_names <- NA
    if(!is.null(gene_names)){
      estim_feature_names <- gene_names
    }

    # get sample names, if present
    estim_sample_names <- NA
    if(!is.null(sample_names)){
      estim_sample_names <- sample_names[cond]
    }

    # get condition name, if present
    estim_cond_name <- paste0("cond_",cond_i)
    if(!is.null(names(conditions))){
      estim_cond_name <- names(conditions)[cond_i]
    }

    print("...creating SPARSim simulation parameter")
    cond_parameter <- SPARSim_create_simulation_parameter (intensity = estim_intensity,
                                                           variability = estim_variability,
                                                           library_size = estim_lib_size,
                                                           feature_names = estim_feature_names,
                                                           sample_names = estim_sample_names,
                                                           condition_name = estim_cond_name)

    dataset_parameter[[cond_i]] <- cond_parameter
    cond_i <- cond_i + 1

  }

  #names(dataset_parameter) <- names(conditions)

  return(dataset_parameter)
}
