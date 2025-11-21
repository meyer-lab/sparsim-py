#' Simulate technical variability
#'
#' Function to simulate the technical variability (i.e. a multivariate hypergeometric on a gamma expression value array)
#'
#' @param avgAbund array containing the intensity values for each feature. It describes the intensity of a single sample
#' @param seqdepth sequencing depth (i.e. sample size of the MH)
#' @param digits number of digits for random number generation
#' @param max_val max value for random number generation
#' @return An array of \code{length(avgAbund)} elements representing the count values for the current sample
simulate_hyper <- function(avgAbund, seqdepth = NULL, digits, max_val) {

  if(seqdepth > 10^digits){
    print("seqdepth must be <= than 10^digits")
    return(NA)
  }

  index <- random_unif_interval(round(seqdepth*1.1),max_val)

  while_i <- 1
  index <- unique(index)
  while (length(index)<seqdepth){
    new_ind <- random_unif_interval( (seqdepth-length(index))*(while_i*10),max_val)
    index<-c(index,new_ind)
    while_i <- while_i +1
    index<-unique(index)
  }
  index <- index[1:seqdepth]

  cumulative_orig<-cumsum(avgAbund)
  check_orig<-findInterval(index,cumulative_orig,left.open = T)+1
  tmp_tabulate <- tabulate(check_orig, nbins = length(avgAbund))
  names(tmp_tabulate)<-names(avgAbund)

  return(tmp_tabulate)
}



#' Simulate technical variability (with control on random number generation through seed)
#'
#' Function to simulate the technical variability (i.e. a multivariate hypergeometric on a gamma expression value array)
#' with control on random number generation through seed.
#'
#' @param avgAbund array containing the intensity values for each feature. It describes the intensity of a single sample
#' @param seqdepth sequencing depth (i.e. sample size of the MH)
#' @param digits number of digits for random number generation
#' @param max_val max value for random number generation
#' @param seed numeric seed
#' @return An array of \code{length(avgAbund)} elements representing the count values for the current sample
simulate_hyper_seed <- function(avgAbund, seqdepth = NULL, digits, max_val, seed) {

  if(seqdepth > 10^digits){
    print("seqdepth must be <= than 10^digits")
    return(NA)
  }

  # index <- random_unif_interval(round(seqdepth*1.1),max_val)
  #
  # while_i <- 1
  # index <- unique(index)
  # while (length(index)<seqdepth){
  #   new_ind <- random_unif_interval( (seqdepth-length(index))*(while_i*10),max_val)
  #   index<-c(index,new_ind)
  #   while_i <- while_i +1
  #   index<-unique(index)
  # }
  # index <- index[1:seqdepth]

  index <- runif_int_no_duplicates(seqdepth, max_val, seed)

  cumulative_orig<-cumsum(avgAbund)
  check_orig<-findInterval(index,cumulative_orig,left.open = T)+1
  tmp_tabulate <- tabulate(check_orig, nbins = length(avgAbund))
  names(tmp_tabulate)<-names(avgAbund)

  return(tmp_tabulate)
}



#' Generate integer number from an uniform distribution with no duplicates
#'
#' Function to generate \code{num} integer values from an uniform distribution \code{[1, max_val]} with no duplicates
#' with control on random number generation through seed.
#'
#' @param num sample size of the MH
#' @param max_val max value for random number generation in the range \code{[1, max_val]}
#' @param seed numeric seed
#' @return An array of \code{length(avgAbund)} elements representing the count values for the current sample
runif_int_no_duplicates <- function(num, max_val, seed){
  return(random_unif_interval_seed(num, max_val, seed)[1:num])
}



#' Function to simulate a raw count table
#'
#' @param dataset_parameter list containing, the intensity, variability and lib sizes of each experimental condition. It is the return value of "estimate_parameter_from_data" or could be created by the users
#' @param batch_parameter batch effect simulation parameter
#' @param spikein_parameter spike-in simulation parameter
#' @param output_sim_param_matrices boolean flag. If TRUE, the function will output two additional matrices, called abundance_matrix and variability_matrix, containing the gene intensities and gene variabilities used as simulation input. (Default: FALSE)
#' @param output_batch_matrix boolean flag. If TRUE, the function will output an additional matrix, called batch_factors_matrix, containing the multiplicative factors used in batch effect simulation. (Default: FALSE)
#' @param gene_expr_simulation_seed integer value used as seed for random number generation during gene expression level simulation. If NULL (default) the current R seed is used.
#' @param count_data_simulation_seed integer value used as seed for random number generation during count data simulation. If NULL (default) the current R and C++ seed is used.
#' @param bimodal_expr_simulation_seed integer value used as seed for random number generation during bimodal gene expression level simulation. If NULL (default) the current R seed is used.
#' @param batch_effect_simulation_seed integer value used as seed for random number generation during batch effect simulation. If NULL (default) the current R seed is used.
#' @param preserve_global_rng boolean flag. If TRUE, the function will not alter the state of the random number generator in the global R environment. If FALSE (default), the state of the random number generator may be altered during the execution of the function.
#' @return A list of 5 elements:
#'
#'   - \code{count_matrix}: the simulated count matrix (genes on rows, samples on columns)
#'
#'   - \code{gene_matrix}: the simulated gene expression levels (genes on rows, samples on columns)
#'
#'   - \code{abundance_matrix}: the input gene intensity values provided as input (genes on rows, samples on columns), if \code{output_sim_param_matrices} = TRUE. NULL otherwise.
#'
#'   - \code{variability_matrix}: the input gene variability values provided as input (genes on rows, samples on columns), if \code{output_sim_param_matrices} = TRUE. NULL otherwise.
#'
#'   - \code{batch_factors_matrix}: the multiplicative factor used in batch generation (genes on rows, samples on columns), if \code{output_batch_matrix} = TRUE. NULL otherwise.
#'
#' @export
SPARSim_simulation <- function(dataset_parameter,
                               batch_parameter = NULL,
                               spikein_parameter = NULL,
                               output_sim_param_matrices = FALSE,
                               output_batch_matrix = FALSE,
                               gene_expr_simulation_seed = NULL,
                               count_data_simulation_seed = NULL,
                               bimodal_expr_simulation_seed = NULL,
                               batch_effect_simulation_seed = NULL,
                               preserve_global_rng = FALSE){

  # should the state of the random number generator in the global R environment be preserved?
  if (preserve_global_rng == TRUE){
    # retrieve the current random number generation state
    current_rng_state <- get_rng_state()
    # restore the saved state before the current function exits
    on.exit(set_rng_state(current_rng_state))
  }



  # number of experimental condition
  N_cond <- length(dataset_parameter)

  # number of genes
  N_genes <- length(dataset_parameter[[1]]$intensity)

  # total number of cells
  N_cell <- sum(unlist(lapply( dataset_parameter, function(x){return(length(x$lib_size))})))

  cat("Number of experimental conditions: ", N_cond, "\n")
  cat("Number of genes: ", N_genes, "\n")
  cat("Number of cells: ", N_cell, "\n")


  cat("Setting gene expression intensity... ", "\n")
  # initialize gene expression level matrix
  gene_expression_matrix <- matrix(0, nrow = N_genes, ncol = N_cell)
  rownames(gene_expression_matrix) <- names(dataset_parameter[[1]]$intensity)
  column_index <- 1
  new_col_names <- ""
  for(cond in 1:N_cond){

    cell_cond_name <- names(dataset_parameter[[cond]]$lib_size)
    new_col_names <- c( new_col_names ,  cell_cond_name)

  }
  new_col_names <- new_col_names[-1]
  colnames(gene_expression_matrix) <- new_col_names


  # fill gene expression level matrix with initial (still without biological variability) gene expression values
  for(cond in 1:N_cond){
     tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
     tmp_N_cell <- length(tmp_cell_ids)

     N_cell_mem <- 5000
     if(tmp_N_cell < N_cell_mem){
       gene_expression_matrix[, tmp_cell_ids] <- rep(dataset_parameter[[cond]]$intensity, tmp_N_cell)
     }else{
       tmp_data <- rep(dataset_parameter[[cond]]$intensity, N_cell_mem)
       N_mult <- floor(tmp_N_cell/N_cell_mem)

       for(ind in c(1:N_mult)){
         ind_range <- c( ((ind-1)*N_cell_mem+1) : (ind*N_cell_mem) )
         gene_expression_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data
       }
       if(N_mult*N_cell_mem < tmp_N_cell){
         ind_range <- c( (N_mult*N_cell_mem+1) : tmp_N_cell)
         gene_expression_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data[c(1: (length(ind_range)*N_genes) )]
       }

       rm(tmp_data); gc(verbose = FALSE)
     }

  }


  cat("Setting gene expression variability ... ", "\n")
  # initialize gene expression variability matrix
  gene_expression_var_matrix <- matrix(0, nrow = N_genes, ncol = N_cell)
  rownames(gene_expression_var_matrix) <- rownames(gene_expression_matrix)
  colnames(gene_expression_var_matrix) <- colnames(gene_expression_matrix)

  # fill gene expression variability matrix
  for(cond in 1:N_cond){
    tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
    tmp_N_cell <- length(tmp_cell_ids)

    # get variability values and fix possible NA or negative values
    variability_values <- dataset_parameter[[cond]]$variability
    variability_values[is.na(variability_values)]<-0
    variability_values[variability_values<0]<-min(variability_values[variability_values>0])

    N_cell_mem <- 5000
    if(tmp_N_cell < N_cell_mem){
      gene_expression_var_matrix[, tmp_cell_ids] <- rep(variability_values, tmp_N_cell)
    }else{
      tmp_data <- rep(variability_values, N_cell_mem)
      N_mult <- floor(tmp_N_cell/N_cell_mem)

      for(ind in c(1:N_mult)){
        ind_range <- c( ((ind-1)*N_cell_mem+1) : (ind*N_cell_mem) )
        gene_expression_var_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data
      }
      if(N_mult*N_cell_mem < tmp_N_cell){
        ind_range <- c( (N_mult*N_cell_mem+1) : tmp_N_cell)
        gene_expression_var_matrix[, tmp_cell_ids[ind_range]  ] <- tmp_data[c(1: (length(ind_range)*N_genes) )]
      }

      rm(tmp_data); gc(verbose = FALSE)
    }

  }



  # if user asks for bimodal gene expression
  if(!is.null(bimodal_expr_simulation_seed)){ # if user set a seed for bimodal gene expression
    # retrieve the current random number generation state
    pre_bimodal_rng_state <- get_rng_state()
    # set the user specified seed
    set.seed(bimodal_expr_simulation_seed)
  }
  for(cond in 1:N_cond){

    if(!is.null(dataset_parameter[[cond]]$intensity_2)){ # check if bimodal gene expression is present in this experimental condition
      tmp_cell_ids <- names(dataset_parameter[[cond]]$lib_size)
      tmp_N_cell <- length(tmp_cell_ids)

      # get intensity values for mode 2
      intensity_values <- dataset_parameter[[cond]]$intensity_2

      # get variability values for mode 2 and fix possible NA or negative values
      variability_values <- dataset_parameter[[cond]]$variability_2
      variability_values[is.na(variability_values)]<-0
      variability_values[variability_values<0]<-min(variability_values[variability_values>0])

      # get probability for mode 2 (as 1 - p_mode)
      p_mode_2 <- 1 - dataset_parameter[[cond]]$p_bimod

      # identify bimodal genes (i.e. having intensities != NA and probability of mode 2 beign greater than 0)
      bimodal_genes_id <- ( ( !is.na(intensity_values) ) & (p_mode_2>0) )
      bimodal_genes_name <- rownames(gene_expression_matrix)[bimodal_genes_id]

      intensity_values_tmp <- intensity_values[bimodal_genes_id]
      variability_values_tmp <- variability_values[bimodal_genes_id]
      p_mode_2_tmp <- p_mode_2[bimodal_genes_id]

      N_bimod_gene <- length(intensity_values_tmp)

      # simulate bimodal genes in the cells
      for(cell in tmp_cell_ids){

        # use the probability to be in mode 2 to select which genes are in the second expression mode in the current cell
        mode2_ind <- as.logical(rbinom(N_bimod_gene, 1, p_mode_2_tmp))

        # only for the selected genes, change the gene expression intensity matrix and the gene expression variability matrix
        # with the values describing the second expression mode
        gene_expression_matrix[bimodal_genes_name[mode2_ind],cell] <- intensity_values_tmp[mode2_ind]
        gene_expression_var_matrix[bimodal_genes_name[mode2_ind],cell] <- variability_values_tmp[mode2_ind]
      }

    }

  }
  if(!is.null(bimodal_expr_simulation_seed)){ # if user set a seed for bimodal gene expression
    # restore the saved rng state before the bimodal gene expression simulation
    set_rng_state(pre_bimodal_rng_state)
  }



  # get samples library sizes
  sample_lib_size <- unlist(lapply( dataset_parameter, function(x){return(x$lib_size)} ), use.names = FALSE)
  names(sample_lib_size) <- colnames(gene_expression_matrix)



  # if user asks for spike-in presence
  if(! is.null(spikein_parameter)){

    cat("Simulating spike-ins presence ... ", "\n")

    # compute the abundances of endogenous material in each experimental condition as the sum of gene expression levels
    end_abund <- unlist( lapply( dataset_parameter, function(x){return(sum(x$intensity))}) )

    # use the median of all the endogenous material abundances as reference for spike-in addition
    reference_abundance_value <- median(end_abund)

    # compute spike-in abundance given the reference value of endogenous material
    spikein_parameter_new <- compute_spikein_abundance(SPARSim_spikein_parameter = spikein_parameter, endogenous_abundance = reference_abundance_value)

    # get IDs and number of spike-in
    spikein_ids <- spikein_parameter_new[["spikein_ids"]]
    N_spikein <- length(spikein_ids)


    # initialize spike-in abundance matrix
    spikein_abund_matrix <- matrix(0, nrow = N_spikein, ncol = N_cell)
    rownames(spikein_abund_matrix) <- spikein_ids
    colnames(spikein_abund_matrix) <- colnames(gene_expression_matrix)


    # fill spike-in abundance matrix
    spikein_abund_matrix <- create_spikein_matrix (spikein_matrix = spikein_abund_matrix, type = "abundance", SPARSim_spikein_parameter = spikein_parameter_new)

    # initialize spike-in variability matrix
    spikein_var_matrix  <- matrix(0, nrow = N_spikein, ncol = N_cell)
    rownames(spikein_var_matrix) <- spikein_ids
    colnames(spikein_var_matrix) <- colnames(gene_expression_matrix)


    # fill spike-in variability matrix
    spikein_var_matrix <- create_spikein_matrix (spikein_matrix = spikein_var_matrix, type = "variability", SPARSim_spikein_parameter = spikein_parameter_new)


    # append spike-in matrices to gene matrices
    gene_expression_matrix <- rbind(gene_expression_matrix[,colnames(gene_expression_matrix)], spikein_abund_matrix[, colnames(gene_expression_matrix)])

    gene_expression_var_matrix <- rbind(gene_expression_var_matrix[, colnames(gene_expression_var_matrix)], spikein_var_matrix[, colnames(gene_expression_var_matrix)])

  }



  # if user asks for batch effects
  if(!is.null(batch_effect_simulation_seed)){ # if user set a seed for batch effect simulation
    # retrieve the current random number generation state
    pre_batch_effect_rng_state <- get_rng_state()
    # set the user specified seed
    set.seed(batch_effect_simulation_seed)
  }
  batch_factor_matrix <- NULL
  if(!is.null(batch_parameter)){

    cat("Simulating batch effects ... ", "\n")

    # compute batch effect factors
    batch_parameter <- SPARSim_compute_batch_effect_factors(N_batch_factors = nrow(gene_expression_matrix),
                                                            batch_factors_id = rownames(gene_expression_matrix),
                                                            SPARSim_batch_parameter = batch_parameter)


    # initialize batch factor matrix with no batch effect
    batch_factor_matrix <- matrix(1, nrow = nrow(gene_expression_matrix), ncol = N_cell)
    rownames(batch_factor_matrix) <- rownames(gene_expression_matrix)
    colnames(batch_factor_matrix) <- colnames(gene_expression_matrix)

    # fill batch factor matrix with the required batch effect factors
    batch_factor_matrix <- SPARSim_create_batch_effect_matrix (batch_factor_matrix = batch_factor_matrix, SPARSim_batch_parameter = batch_parameter)

    # apply batch effect factors
    feature_id <- rownames(gene_expression_matrix)
    cell_id <- colnames(gene_expression_matrix)

    gene_expression_matrix[feature_id,cell_id] <- gene_expression_matrix[feature_id,cell_id] * batch_factor_matrix[feature_id,cell_id]

  }
  if(!is.null(batch_effect_simulation_seed)){ # if user set a seed for batch effect simulation
    # restore the saved rng state before the batch effect simulation
    set_rng_state(pre_batch_effect_rng_state)
  }



  # if user does not require the batch factor matrix as output, free the memory and return a NULL matrix
  if(output_batch_matrix == FALSE){
    rm(batch_factor_matrix); gc(verbose = FALSE);
    batch_factor_matrix <- NULL
  }



  ### Simulate biological variability
  cat("Simulating biological variability ... ", "\n")

  # if user set a seed for gene expression simulation
  if(!is.null(gene_expr_simulation_seed)){
    # retrieve the current random number generation state
    pre_gene_expr_rng_state <- get_rng_state()
    # set the user specified seed
    set.seed(gene_expr_simulation_seed)
  }

  gene_expression_matrix_bio_var <- gene_expression_matrix

  # simulate biological variability using a gamma
  gene_expression_matrix_bio_var <- matrix(
                                          rgamma(n = nrow(gene_expression_matrix_bio_var)*ncol(gene_expression_matrix_bio_var),
                                                shape = 1/gene_expression_var_matrix,
                                                scale = gene_expression_var_matrix*gene_expression_matrix),
                                          ncol = ncol(gene_expression_matrix)
                                          )

  # feature having null (i.e. zero) variability have no variability, so "undo" the gamma
  zero_var_index <- (gene_expression_var_matrix == 0)
  gene_expression_matrix_bio_var [zero_var_index] <- gene_expression_matrix [zero_var_index]

  rownames(gene_expression_matrix_bio_var) <- rownames(gene_expression_matrix)
  colnames(gene_expression_matrix_bio_var) <- colnames(gene_expression_matrix)

  # if user does not require the gene intensity and gene variability matrices as output, free the memory and return NULL matrices
  if(output_sim_param_matrices == FALSE){
    rm(gene_expression_matrix, gene_expression_var_matrix); gc(verbose = FALSE)
    gene_expression_matrix <- NULL
    gene_expression_var_matrix <- NULL
  }

  # identify the maximum library size
  max_lib_size <- max(sample_lib_size)

  # set the input fragment library size (the population size parameter for MH)
  input_fragment_lib_size <- max_lib_size * 10^2
  digits<-floor(log10(input_fragment_lib_size))+1
  new_libsize<-10^digits

  gene_expression_matrix_bio_var_scaled <- round(t(t(gene_expression_matrix_bio_var)*new_libsize/colSums(gene_expression_matrix_bio_var)))
  num_fragment <- colSums(gene_expression_matrix_bio_var_scaled)

  # if user set a seed for gene expression simulation, restore the saved rng state before the gene expression simulation
  if(!is.null(gene_expr_simulation_seed)){
    set_rng_state(pre_gene_expr_rng_state)
  }


  ### Simulate technical variability
  cat("Simulating technical variability ... ", "\n")

  sim_count_matrix <- matrix(0, ncol = ncol(gene_expression_matrix_bio_var_scaled), nrow = nrow(gene_expression_matrix_bio_var_scaled) )
  rownames(sim_count_matrix) <- rownames(gene_expression_matrix_bio_var_scaled)
  colnames(sim_count_matrix) <- colnames(gene_expression_matrix_bio_var_scaled)

  # if the user set NO seed for count data simulation
  if(is.null(count_data_simulation_seed)){
    for(sample in colnames(sim_count_matrix)){
      sim_count_matrix[, sample] <- simulate_hyper(avgAbund = gene_expression_matrix_bio_var_scaled[, sample],
                                                   seqdepth = sample_lib_size[sample],
                                                   digits = digits,
                                                   max_val = num_fragment[sample])
    }
  }
  # if the user set a seed for count data simulation
  else{

    #retrieve the current random number generation state
    pre_count_data_rng_state <- get_rng_state()
    #set the user specified seed
    set.seed(count_data_simulation_seed)

    # create a seed for each cell
    count_data_cell_seed <- stats::runif(ncol(sim_count_matrix), min = 0, max = 2147483647)
    names(count_data_cell_seed) = colnames(sim_count_matrix)

    for(sample in colnames(sim_count_matrix)){
      sim_count_matrix[, sample] <- simulate_hyper_seed(avgAbund = gene_expression_matrix_bio_var_scaled[, sample],
                                                        seqdepth = sample_lib_size[sample],
                                                        digits = digits,
                                                        max_val = num_fragment[sample],
                                                        seed = count_data_cell_seed[sample])
    }

    #restore the saved rng state before the count data simulation step
    set_rng_state(pre_count_data_rng_state)

  }

  ### collect and return the results
  return(list(count_matrix = sim_count_matrix,
              gene_matrix = gene_expression_matrix_bio_var,
              abundance_matrix = gene_expression_matrix,
              variability_matrix = gene_expression_var_matrix,
              batch_factors_matrix = batch_factor_matrix )
  )
}


