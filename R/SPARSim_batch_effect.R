#' Create a SPARSim batch object
#'
#' Function to create a SPARSim batch object.
#'
#' A SPARSim batch object describe the characteristics of the batch effect to apply during SPARSim simulation.
#' Input parameter \code{name} contains a batch identifier (e.g. "Lane_1", "Lane_2", etc.).
#' Batch effect is simulated multiplying gene expression levels by batch specific multiplicative factors.
#' Batch specific multiplicative factors will be sampled from the distribution specified in the
#' input parameter \code{distribution}, with parameter \code{param_A} and \code{param_B}.
#' Currently, distribution could be one between \code{"normal"} (i.e. \code{rnorm(..., mean = param_A, sd = param_B)})
#' and \code{"gamma"} (i.e. \code{rgamma(..., scale = param_A, shape = param_B)})
#'
#' @param name name of the batch
#' @param distribution batch effect factors distribution. current available distributions:
#' \itemize{
#'   \item "normal": normal distribution (mean = param_A, sd = param_B)
#'   \item "gamma": gamma distribution (scale = param_A, shape = param_B)
#' }
#' @param param_A parameter to describe the distribution
#' @param param_B parameter to describe the distribution
#' @return A SPARSim batch object
#'
#' @examples
#' # create 2 batches, called "Lane1" and "Lane2"
#' batch_1 <- SPARSim_create_batch(name = "Lane1", distribution = "normal", param_A = 1, param_B = 1)
#' batch_2 <- SPARSim_create_batch(name = "Lane2", distribution = "normal", param_A = 3, param_B = 1)
#'
#' @export
SPARSim_create_batch <- function(name, distribution = "normal", param_A = 1, param_B = 1){
  batch <- list()
  batch[["name"]] <- name
  batch[["param_A"]] <- param_A
  batch[["param_B"]] <- param_B
  batch[["distribution"]] <- distribution
  return(batch)
}


#' Create a SPARSim batch set
#'
#' Function to create a set of SPARSim batches
#'
#' @param batch_list a list of SPARSim batch objects. Each element of the list must be the result of a call to function \code{\link{SPARSim_create_batch}}
#'
#' @return A SPARSim batch set
#'
#' @examples
#' # create 2 batches, called "Lane1" and "Lane2"
#' batch_1 <- SPARSim_create_batch(name = "Lane1", distribution = "normal", param_A = 1, param_B = 1)
#' batch_2 <- SPARSim_create_batch(name = "Lane2", distribution = "normal", param_A = 3, param_B = 1)
#'
#' # create batch set
#' batches <- list(); batches[[1]] <- batch_1; batches[[2]] <- batch_2
#' batch_set <- SPARSim_create_batch_set(batch_list = batches)
#'
#' @export
SPARSim_create_batch_set <- function(batch_list){

  if(is.null(batch_list)){
    stop("ERROR: the input parameter of function SPARSim_create_batch_set() is NULL. It should be a list of batches; each element in the list should be created using SPARSim_create_batch()")
  }

  N_batch_type <- length(batch_list)

  batch_names <- batch_list[[1]][["name"]]

  if(N_batch_type>1){
    for(i in 2:N_batch_type){
      batch_names <- c(batch_names, batch_list[[i]][["name"]])
    }
  }

  if(length(unique(batch_names)) < N_batch_type){
    stop("ERROR: in the input parameter of function SPARSim_create_batch_set(), there are batches having the same name (field batch_list[[i]][\"name\"]. Batch name must be unique.")
  }

  names(batch_list) <- batch_names

  return(batch_list)
}


#' Create SPARSim batch effect simulation parameter
#'
#' Function to create the SPARSim input batch parameter.
#'
#' @param batch_set result of function \code{\link{SPARSim_create_batch_set}}
#' @param batch_sample a character array of length N, where N is the number of cells/samples to simulate.
#' The i-th element in the array contains the name of the batch type to simulate in the i-th sample;
#' if the i-th element is set to NULL, then the i-th sample is affected by no batch effect
#'
#' @return A SPARSim batch effect simulation parameter
#' @examples
#' # create 2 batches, called "Lane1" and "Lane2"
#' batch_1 <- SPARSim_create_batch(name = "Lane1", distribution = "normal", param_A = 1, param_B = 1)
#' batch_2 <- SPARSim_create_batch(name = "Lane2", distribution = "normal", param_A = 3, param_B = 1)
#'
#' # create batch set
#' batches <- list(); batches[[1]] <- batch_1; batches[[2]] <- batch_2
#' batch_set <- SPARSim_create_batch_set(batch_list = batches)
#'
#' # create the association between batches and samples/cells
#' # as example, consider simulating 40 cells:
#' # the first 30 samples are affected by batch "Lane1",
#' # the second 10 samples are affected by batch "Lane2"
#' batch_sample <- c(rep("Lane1", 30),rep("Lane2", 10))
#'
#' # having both the batch set (i.e. batch_set) and
#' # the association between batches and samples/cells (i.e. batch_sample)
#' # it is now possible to create the SPARSim batch effect simulation parameter
#' SPARSim_batch_parameter <- SPARSim_create_batch_parameter(batch_set = batch_set,
#'                                                           batch_sample = batch_sample)
#'
#' @export
SPARSim_create_batch_parameter <- function(batch_set, batch_sample){

  SPARSim_batch_parameter <- list()
  SPARSim_batch_parameter[["batch_set"]] <- batch_set
  SPARSim_batch_parameter[["batch_sample"]] <- batch_sample

  return(SPARSim_batch_parameter)
}



#' Generate batch effect multiplicative factors for each batch
#'
#' Function to generate the batch effect multiplicative factors for each batch, following the batch description
#' in \code{SPARSim_batch_parameter}.
#'
#' @param N_batch_factors number of features that will be affected by batch effect multiplicative factors, i.e. number of features in the count matrix
#' @param batch_factors_id IDs of features that will be affected by batch effect multiplicative factors, i.e. names of the features in the count matrix
#' @param SPARSim_batch_parameter SPARSim batch parameter containing the information about the batch effect
#'
#' @return An updated version of the input parameter \code{SPARSim_batch_parameter} containing the values of the multiplicative factors for each batch
#'
SPARSim_compute_batch_effect_factors <- function(N_batch_factors, batch_factors_id, SPARSim_batch_parameter){

  # get the set of batches
  batch_set <- SPARSim_batch_parameter[["batch_set"]]

  batch_effect_factors <- list()

  N_batch <- length(batch_set)

  for(i in 1:N_batch){

    batch_tmp <- batch_set[[i]]

    if(batch_tmp[["distribution"]] == "normal"){
      batch_factors <- rnorm(N_batch_factors, mean = batch_tmp[["param_A"]], sd = batch_tmp[["param_B"]])
      batch_factors[batch_factors<=0] <- batch_tmp[["param_A"]] # to avoid values <= 0, set them to the mean value
      names(batch_factors) <- batch_factors_id
      batch_tmp[["factors"]] <- batch_factors
    }

    if(batch_tmp[["distribution"]] == "gamma"){
      batch_factors <- rgamma(N_batch_factors, scale = batch_tmp[["param_A"]], shape = batch_tmp[["param_B"]])
      batch_factors[batch_factors<=0] <- batch_tmp[["param_A"]]*batch_tmp[["param_B"]] # to avoid values <= 0, set them to the mean value (scale*shape)
      names(batch_factors) <- batch_factors_id
      batch_tmp[["factors"]] <- batch_factors
    }

    batch_set[[i]] <- batch_tmp
  }

  SPARSim_batch_parameter[["batch_set"]] <- batch_set

  return(SPARSim_batch_parameter)
}




#' Create batch effect multiplicative factors matrix
#'
#' Function to create batch effects multiplicative factors matrix. This function is used internally by \code{\link{SPARSim_simulation}}.
#'
#' This function return a matrix containing the batch effect multiplicative factors.
#' The matrix has the same dimensions and names of \code{batch_factor_matrix}.
#' The element \code{[i,j]} in the matrix contains the batch effect factor for i-th feature in j-th sample.
#'
#' @param batch_factor_matrix matrix to fill with batch factor values (features on rows; samples on columns)
#' @param SPARSim_batch_parameter SPARSim batch parameter containing the information about the batch effect
#'
#' @return A numeric matrix, having the same dimensions and names of \code{batch_factor_matrix}.
#' The element \code{[i,j]} in the matrix contains the batch effect factors for i-th feature in j-th sample.
SPARSim_create_batch_effect_matrix <- function(batch_factor_matrix, SPARSim_batch_parameter){

  # get the set of batches
  batch_set <- SPARSim_batch_parameter[["batch_set"]]

  # get association between sample and batches
  batch_sample <- SPARSim_batch_parameter[["batch_sample"]]

  for(batch in batch_set){

    batch_name <- batch[["name"]] # get batch name
    batch_factors <- batch[["factors"]] # get batch factor values

    cell_id <- batch_sample == batch_name # get cell IDs related to the current batch
    cell_id[is.na(cell_id)] <- FALSE

    batch_factor_matrix[, cell_id] <- rep(batch_factors, sum(cell_id)) # fill matrix by column

  }

  return(batch_factor_matrix)

}







