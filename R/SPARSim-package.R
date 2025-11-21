#' SPARSim single cell: a count data simulator for scRNA-seq data
#'
#' This package implements the SPARSim simulator, a tool to simulate single cell RNA-seq count table.
#' For additional details about the simulation model, please refer to
#' \emph{Baruzzo, G., Patuzzi, I., & Di Camillo, B. (2020). SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, 36(5), 1468-1475.}
#'
#' @docType package
#'
#' @author Giacomo Baruzzo \email{giacomo.baruzzo@@unipd.it}
#' @author Ilaria Patuzzi \email{IPatuzzi@@izsvenezie.it}
#' @author Barbara Di Camillo \email{barbara.dicamillo@@unipd.it}
#'
#' @references Baruzzo, G., Patuzzi, I., & Di Camillo, B. (2020). SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, 36(5), 1468-1475.
#'
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats median rbinom rgamma rnorm approx
#' @useDynLib SPARSim
#' @name SPARSim
#'
NULL
