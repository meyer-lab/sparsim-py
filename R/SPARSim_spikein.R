#' Create a spike-in mix
#'
#' Function to create a spike-in mix
#'
#' A spike-in mix is a set of spike-ins. Each spike-in mix is identified by a name (input parameter \code{name}).
#'
#' The abundances of spike-ins in the mix are described in the parameter \code{abundance}, a numeric vector containing
#' the abundance of each spike-in in the spike-in mix.
#' The length of \code{abundance} defines the number of spike-in in the spike-in mix.
#'
#' By default, spike-ins are named "spikein_1", "spikein_2", ..., . The user could specified different spike-in IDs using
#' the parameter \code{spike_in_IDS}, which should contain an array of spike-in IDs having the same number of element of \code{abundance}.
#'
#' Spike-ins should be affected only by technical variability, since their are synthetic sequences added at kwown concentration.
#' However, in some studies an extra variability was detected.
#' Such extra variability could be optionally added using the \code{extra_variability} input parameter, having the same meaning of the \code{variability}
#' parameter used in function \code{\link{SPARSim_create_simulation_parameter}}.
#' By default, no extra variability is simulated.
#'
#' @param mix_name name of the spike-in mix
#' @param abundance numeric array, the i-th element contains the abundance of the i-th spike-in.
#' @param spike_in_IDS optional array of spike-in IDs; if NA (default) spike-ins will be named "spikein_1", "spikein_2", ...,
#' @param extra_variability optional array of extra variability values; if NA (default), no extra variability is added during spike-ins simulation
#'
#' @return A spike-in mix object
#'
#' @examples
#' # create a first spike-in mix made by 100 spike-ins. Let's call it "spikein_M1"
#' # and set spike-in IDs to "spikein_M1_1", "spikein_M1_2", ..., "spikein_M1_100"
#' spikein_mix1_abundance <- runif(n = 100, 0.01, 1000)
#' spikein_mix1_IDS <- paste0("spikein_M1_", c(1:100))
#' spikein_mix1 <- SPARSim_create_spikein_mix(mix_name= "spikein_M1",
#'                                            abundance = spikein_mix1_abundance,
#'                                            spike_in_IDS = spikein_mix1_IDS)
#'
#' # create a second spike-in mix made by 90 spike-ins
#' # let's call it "spikein_M2" and add some extra variability to it
#' spikein_mix2_abundance <- runif(n = 90, 0.001, 10000)
#' spikein_mix2_extra_var <- runif(n = 90, 0.01, 0.02)
#' spikein_mix2 <- SPARSim_create_spikein_mix(mix_name= "spikein_M2",
#'                                            abundance = spikein_mix2_abundance,
#'                                            extra_variability = spikein_mix2_extra_var)
#'
#' @export
SPARSim_create_spikein_mix <- function(mix_name, abundance, spike_in_IDS = NA, extra_variability = NA){

  # get the number of spike-in
  N_spikein <- length(abundance)

  # if no extra variability values are provided by the user, set the variability parameter to 0
  if(length(extra_variability)==1){
    if(is.na(extra_variability)){
      extra_variability <- rep(0, N_spikein)
    }
  }

  # if no spike-in IDs are provided, spike-in will be named "spikein_1", "spikein_2", ..., "spikein_<N_spikein>"
  if(length(spike_in_IDS)==1){
    if(is.na(spike_in_IDS)){
      spike_in_IDS <- paste0("spikein_", c(1:N_spikein))
    }
  }

  # assign spikein IDs to values
  names(abundance) <- names(extra_variability) <- spike_in_IDS

  spikein <- list()
  spikein[["mix_name"]] <- mix_name
  spikein[["abundance"]] <- abundance
  spikein[["variability"]] <- extra_variability
  spikein[["ids"]] <- spike_in_IDS

  return(spikein)
}



#' Create a set of spike-in mixes
#'
#' Function to create a pool of spike-in mixes
#'
#'
#' If different spike-in mixes have spike-ins with the same names, such spike-ins are considered to be the same spike-ins, but maintaining their original values.
#'
#' As a toy example, consider a spike-in mix A made by 5 spike-ins "spike_1", "spike_2", "spike_3", "spike_4", and "spike_5" with abundances 10, 20, 30, 40 and 50, respectively.
#' Consider a spike-in mix B made by 4 spike-ins "spike_4", "spike_5", "spike_6" and "spike_7", with abundances 400, 500, 600 and 700 respectively.
#'
#' The set of spike-in mixes is so made by 7 spike-ins: "spike_1", "spike_2", "spike_3", "spike_4", "spike_5", "spike_6" and "spike_7".
#' In this scenario, "spike_4" and "spike_5" are common to spike-in mix A and mix B.
#'
#' However, they will maintain their original values: "spike_4" and "spike_5" will have abundances of 40 and 50 when associated to spike-in mix A,
#' while they will have abundances of 400 and 500 when associated to spike-in mix B.
#'
#' @param spikein_mixes a list of existing spike-in mixes. Each element of the list must be the result of a call to function \code{\link{SPARSim_create_spikein_mix}}
#'
#' @return A set of spike-in mixes
#'
#' @examples
#' # create a first spike-in mix made by 100 spike-ins. Let's call it "spikein_M1"
#' # and set spike-in IDs to "spikein_M1_1", "spikein_M1_2", ..., "spikein_M1_100"
#' spikein_mix1_abundance <- runif(n = 100, 0.01, 1000)
#' spikein_mix1_IDS <- paste0("spikein_M1_", c(1:100))
#' spikein_mix1 <- SPARSim_create_spikein_mix(mix_name= "spikein_M1",
#'                                            abundance = spikein_mix1_abundance,
#'                                            spike_in_IDS = spikein_mix1_IDS)
#'
#' # create a second spike-in mix made by 90 spike-ins
#' # let's call it "spikein_M2" and add some extra variability to it
#' spikein_mix2_abundance <- runif(n = 90, 0.001, 10000)
#' spikein_mix2_extra_var <- runif(n = 90, 0.01, 0.02)
#' spikein_mix2 <- SPARSim_create_spikein_mix(mix_name= "spikein_M2",
#'                                            abundance = spikein_mix2_abundance,
#'                                            extra_variability = spikein_mix2_extra_var)
#'
#' # create set of spike-in mixes
#' spikein_mixes_list <- list(mix1 = spikein_mix1, mix2 = spikein_mix2)
#' spikein_set <- SPARSim_create_spikein_set(spikein_mixes = spikein_mixes_list)
#'
#' @export
SPARSim_create_spikein_set <- function(spikein_mixes){

  if(is.null(spikein_mixes)){
    stop("ERROR: the input parameter of function create_spikein_set() is NULL. It should be a list of spike-in mix; each element in the list should be created using create_spikein()")
  }

  N_spikein_mix <- length(spikein_mixes)

  spikein_mix_names <- spikein_mixes[[1]][["mix_name"]]

  if(N_spikein_mix>1){
    for(i in 2:N_spikein_mix){
      spikein_mix_names <- c(spikein_mix_names, spikein_mixes[[i]][["mix_name"]])
    }
  }

  if(length(unique(spikein_mix_names)) < N_spikein_mix){
    stop("ERROR: in the input parameter of function create_spikein_set(), there are spike-in mix having the same name (field spikein_list[[i]][\"name\"]. Spike-in mix name must be unique.")
  }

  names(spikein_mixes) <- spikein_mix_names

  return(spikein_mixes)

}



#' Create SPARSim spike-in simulation parameter
#'
#' Function to create a SPARSim spike-in simulation parameter
#'
#' @param spikein_set a set of spike-in mixes, as result of a call to function \code{\link{SPARSim_create_spikein_set}}
#' @param spikein_sample a character array of length equal to the number of samples.
#' The i-th element in the array contains the name of the spike-in mix to simulate in the i-th sample;
#' if the i-th element is set to NULL, then the i-th sample will contain no spike-in mix.
#' @param spikein_proportion a numeric array of length \code{length(spikein_set)}, the i-th element in the array is associated to the i-th spike-in mix in \code{spikein_set}.
#' Given a reference amount of endogenous material X (here we'll set X to the median of the amount of endogenous material among all the input samples in all the experimental conditions),
#' the i-th element in the array represents the proportion (i.e. a value in [0,1]) of the i-th spike-in mix material respect to the reference amount X.
#' As example, if \code{spikein_set} contains 2 spike-in mixes, e.g. "spikein_M1" and "spikein_M2", having spikein_proportion set to \code{c(0.03,0.05)}
#' and having the reference amount of endogenous material X = 100ng, then:
#' \itemize{
#'   \item calling Y1 the simulated amount of material from spike-in mix "spikein_M1", Y1 = 0.03*100ng = 3ng of spike-in mix "spikein_M1" will be add to samples which require spike-in mix "spikein_M1"
#'   \item calling Y2 the simulated amount of material from spike-in mix "spikein_M2", Y2 = 0.05*100ng = 5ng of spike-in mix "spikein_M2" will be add to samples which require spike-in mix "spikein_M2"
#' }
#'
#' @return A SPARSim spike-in simulation parameter
#'
#' @examples
#' # create a first spike-in mix made by 100 spike-ins. Let's call it "spikein_M1"
#' # and set spike-in IDs to "spikein_M1_1", "spikein_M1_2", ..., "spikein_M1_100"
#' spikein_mix1_abundance <- runif(n = 100, 0.01, 1000)
#' spikein_mix1_IDS <- paste0("spikein_M1_", c(1:100))
#' spikein_mix1 <- SPARSim_create_spikein_mix(mix_name= "spikein_M1",
#'                                            abundance = spikein_mix1_abundance,
#'                                            spike_in_IDS = spikein_mix1_IDS)
#'
#' # create a second spike-in mix made by 90 spike-ins
#' # let's call it "spikein_M2" and add some extra variability to it
#' spikein_mix2_abundance <- runif(n = 90, 0.001, 10000)
#' spikein_mix2_extra_var <- runif(n = 90, 0.01, 0.02)
#' spikein_mix2 <- SPARSim_create_spikein_mix(mix_name= "spikein_M2",
#'                                            abundance = spikein_mix2_abundance,
#'                                            extra_variability = spikein_mix2_extra_var)
#'
#' # create set of spike-in mixes
#' spikein_mixes_list <- list(mix1 = spikein_mix1, mix2 = spikein_mix2)
#' spikein_set <- SPARSim_create_spikein_set(spikein_mixes = spikein_mixes_list)
#'
#' # create the association between samples/cells and spike-in mixes
#' # as an example, consider having 40 samples
#' # the first 10 samples will contain spike-in mix "spikein_M1",
#' # the following 10 samples will contain no spike-in mix
#' # and the last 20 samples will contain spike-in mix "spikein_M2"
#' spikein_sample <- c(rep("spikein_M1", 10), rep(NA, 10), rep("spikein_M2", 20))
#'
#' # define the proportion of the spike-in mixes respect to the endogenous material
#' spikein_proportion <- c(0.03,0.05)
#'
#' # having the set of spike-in mixes, the association between samples
#' # and spike-in mixes and the proportion of spike-ins abundances respect to the endogenous material,
#' # it is now possible to create the SPARSim spike-in simulation parameter
#' SPARSim_spikein_param <- SPARSim_create_spikein_parameter(spikein_set = spikein_set,
#'                                                           spikein_sample = spikein_sample,
#'                                                           spikein_proportion = spikein_proportion)
#'
#' @export
SPARSim_create_spikein_parameter <- function(spikein_set, spikein_sample, spikein_proportion){

  SPARSim_spikein_parameter <- list()
  SPARSim_spikein_parameter[["spikein_set"]] <- spikein_set
  SPARSim_spikein_parameter[["spikein_sample"]] <- spikein_sample
  SPARSim_spikein_parameter[["spikein_proportion"]] <- spikein_proportion

  N_spikein_mix <- length(spikein_set)

  # spike-in IDs
  spikein_id <- spikein_set[[1]][["ids"]]
  if(N_spikein_mix > 1){
    for(i in 2:N_spikein_mix){
      spikein_id <- union(spikein_id, spikein_set[[i]][["ids"]])
    }
  }
  SPARSim_spikein_parameter[["spikein_ids"]] <- spikein_id

  return(SPARSim_spikein_parameter)
}






### Function to compute the "absolute" abundance of spike-in, given a reference amount of endogenous material (endogenous_abundance)
# param:
#   - SPARSim_spikein_parameter: SPARSim spike-in simulation parameter (created using create_SPARSim_spikein_parameter())
#   - endogenous_abundance: numeric value representing the reference amount of endogenous material
# return: a new SPARSim_spikein_parameter containing the new spike-in "absolute" abundance
compute_spikein_abundance <- function(SPARSim_spikein_parameter, endogenous_abundance){

  # Create the new SPARSim_spikein_parameter starting from the input one
  SPARSim_spikein_parameter_new <- SPARSim_spikein_parameter

  # Extract the spikein mixes and their proportion related to the reference amount of endogenous material
  spikein_set <- SPARSim_spikein_parameter_new[["spikein_set"]]
  spikein_proportion <- SPARSim_spikein_parameter_new[["spikein_proportion"]]

  # Number of spike-in mixes
  N_spikein_mix <- length(spikein_set)

  # Rescale the spike-in abundance to the specified value
  for (i in 1:N_spikein_mix){
    i_th_spikein_mix_values <- spikein_set[[i]][["abundance"]]

    scaled_values <- (i_th_spikein_mix_values / sum(i_th_spikein_mix_values)) * (spikein_proportion[i] * endogenous_abundance)

    spikein_set[[i]][["abundance"]] <- scaled_values
  }

  SPARSim_spikein_parameter_new[["spikein_set"]] <- spikein_set

  return(SPARSim_spikein_parameter_new)

}



### Function to add spike-in to a single experimental condition
# param:
#   - spikein_matrix: spikein matrix to fill with abundance or variability values (features on rows; samples on columns)
#   - type: "abundance" or "variability". it define if the matrix to fill is the one containing the spike-in abundance or variability
#   - SPARSim_spikein_parameter: SPARSim spike-in simulation parameter (created using create_SPARSim_spikein_parameter())
# return: spike-in abundance or variability matrix
create_spikein_matrix <- function(spikein_matrix, type, SPARSim_spikein_parameter){

  spikein_set <- SPARSim_spikein_parameter[["spikein_set"]]
  spikein_sample <- SPARSim_spikein_parameter[["spikein_sample"]]

  # spike-in IDs
  spikein_id <- SPARSim_spikein_parameter[["spikein_ids"]]

  # Number of spike-in mixes
  N_spikein_mix <- length(spikein_set)


  for(spikein_mix in spikein_set){# for each spike-in mix

    spikein_mix_name <- spikein_mix[["mix_name"]] # get spike-in mix name

    spikein_mix_value <- spikein_mix[[type]] # get spike-in abundance or variability values (depending on "type")

    cell_id <- spikein_sample == spikein_mix_name # get cell IDs related to the current spike-in mix
    cell_id[is.na(cell_id)] <- FALSE

    spikein_value_id <- spikein_mix[["ids"]] # get spike-in IDs (to avoid ID sorting problem)

    spikein_matrix[spikein_value_id, cell_id] <- rep(spikein_mix_value[spikein_value_id], sum(cell_id)) # fill matrix by column

  }

  return(spikein_matrix)

}



### Function to add spike-in mixes to a single experimental condition
# param:
#   - fragment_matrix: matrix containing the number of fragment for each feature in each sample (features on rows; samples on columns)
#   - spikein_matrix: matrix containing the spike-in abundances as result of function create_spikein_abundances()
# return: matrix containing the number of fragment for each feature in each sample (features on rows; samples on columns). feature now includes the added spike-in fragments
add_spikein <- function(fragment_matrix, spikein_matrix){
  # add spike-in to the fragment matrix
  final_count_matrix <- rbind (fragment_matrix, spikein_matrix)

  return(final_count_matrix)
}

