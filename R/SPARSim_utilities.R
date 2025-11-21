#' scran normalization
#'
#' Function to compute the count normalization proposed by Lun et. al in the scran package
#'
#' @param data a raw count matrix (genes on rows, samples on columns) to normalize
#' @param sizes parameter \code{size} of function \code{scran::computeSumFactors}
#' @param positive parameter \code{positive} of function \code{scran::computeSumFactors}
#' @return The normalized count matrix using the scran method
#' @export
scran_normalization <- function (data, sizes = seq(20, 100, 5), positive = FALSE){

  # load scran
  #library(scran)

  # create SingleCellExperiment object
  sce_scran <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(data)))

  # do scran normalization
  sce_scran <- scran::computeSumFactors(sce_scran, sizes=sizes, positive=positive) # compute normalization factors
  sce_scran_norm <- scater::normalizeCounts(sce_scran, log = FALSE) # apply normalization factors

  return(sce_scran_norm)
}




#' Associate gene intensity and gene variability values
#'
#' This function associates variability values to a set of intensity values (obtained, for example, during DE genes simulation)
#' \code{intensity_values}.
#' The association between gene intensity values and gene variability values is learned from the one
#' described in an existing simulation parameter \code{sim_param}, using a linear interpolation.
#'
#' @param intensity_values the (new) intensity values to which we want to associate variability values
#' @param sim_param simulation parameter from which learn the gene intensity vs gene variability association
#'
#' @return A vector of the same size of \code{intensity_values} containing the computed variability values
int_var_relation<-function(intensity_values, sim_param){

  # to learn the association between gene intensity and gene variability in sim_param
  # first exclude intensity values associate to NA variability values (i.e. flat genes)
  # second sort intensities values
  # last re-order variability values using the order of the intensity values
  not_NA_var_ind <- which(!is.na(sim_param$variability))
  int_ord <- sort(sim_param$intensity[not_NA_var_ind])
  var_not_NA <- sim_param$variability[not_NA_var_ind]
  var_ord <- var_not_NA[names(int_ord)]


  interp_variability<-vector()

  # handle extreme variability values
  interp_variability[intensity_values >= max(int_ord)] <- sim_param$variability[which.max(int_ord)]
  interp_variability[intensity_values <= min(int_ord)] <- sim_param$variability[which.min(int_ord)]


  # do interpolation
  to_interp <- which(is.na(interp_variability))
  for (k in to_interp){

    if(intensity_values[k] == 0){
      interp_variability[k] <- NA
    }else{
      interpolation_ind <- which.max(int_ord>=intensity_values[k])
      interp_variability[k] <- approx(c(int_ord[interpolation_ind-1],int_ord[interpolation_ind]),
                                      c(var_ord[interpolation_ind-1],var_ord[interpolation_ind]),
                                      xout=intensity_values[k])$y
    }
  }
  return(interp_variability)
}

#' Function to retrieve the current state of the random number generation
#'
#' This function return the \code{.Random.seed} object in the global environment.
#'
#' @return The \code{.Random.seed} object if it exists, NULL otherwise.
get_rng_state <- function(){
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

#' Function to restore the a (previously saved) state of the random number generation
#'
#' This function restore the \code{.Random.seed} object in the global environment with the one passed as input argument.
#' If the input argument is NULL, then no action is performed.
#'
#' @param rng_state A (previously saved) state of the random number generation
set_rng_state <- function(rng_state) {
  if (!is.null(rng_state)) {
    assign(".Random.seed", rng_state, envir = .GlobalEnv, inherits = FALSE)
  }
}
