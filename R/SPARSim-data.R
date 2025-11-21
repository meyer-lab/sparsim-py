#' Example count matrix
#'
#' Example count matrix having 5000 genes and 150 cells (across 3 different experimental conditions).
#' First 50 cells correspond to condition A, second 50 cells correspond to condition B and last 50 cells correspond to conditions C.
#'
#' @docType data
#'
#' @usage data(Example_count_matrix)
#'
#' @format A dataframe having genes on rows and cells on columns.
#'
#' @examples
#' # Load example count matrix
#' data(Example_count_matrix)
#'
#' # Normalize the count matrix
#' Example_count_matrix_norm <- scran_normalization(Example_count_matrix)
#'
#' # Get column index for each experimental condition
#' cond_A_column_index <- c(1:50) # Condition A column indices: from column 1 to column 50
#' cond_B_column_index <- c(51:100) # Condition B column indices: from column 51 to column 100
#' cond_C_column_index <- c(101:150) # Condition C column indices: from column 101 to column 150
#'
#' # Create conditions parameter
#' Example_count_matrix_conditions <- list(cond_A = cond_A_column_index,
#'                                         cond_B = cond_B_column_index,
#'                                         cond_C = cond_C_column_index)
#'
#' # Create SPARSim simulation parameter through the estimation from the example count matrix
#' SPARSim_example_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = Example_count_matrix,
#'                                                        norm_data = Example_count_matrix_norm,
#'                                                        conditions = Example_count_matrix_conditions)
#'
#' # Run SPARSim simulation using the just created simulation parameter
#' example_sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_example_sim_param)
#'
"Example_count_matrix"




#' SPARSim simulation parameter preset from Tung dataset
#'
#' SPARSim simulation parameter preset estimated from Tung (\emph{Tung et al., 2017; GSE77288}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Tung count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 18464 genes
#'   \item 564 cells (across 8 different experimental conditions)
#'   \item average read/UMI count per cell of 70473
#'   \item sparsity ~53\%
#' }
#' Information about the original Tung dataset:
#' \itemize{
#'   \item Number of genes: 18464
#'   \item Number of cells: 564
#'   \item Number of different experimental conditions: 8 (Tung_C1, Tung_C2, Tung_C3, Tung_C4, Tung_C5, Tung_C6, Tung_C7, Tung_C8)
#'   \item Sparsity: ~53.16\%
#'   \item Average read/UMI count per cell: 70473
#'   \item Species: Human
#'   \item Cell types: Induced pluripotent stem cells
#'   \item Protocol and platform: Fluidigm C1, modified SMARTer
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the Tung dataset:
#'
#' scRNA-seq data were produced from induced pluripotent stem cell (iPSC) lines of three Yoruba individuals (NA19098, NA19101 and NA19239).
#' Three C1 96 well-integrated fluidic circuit replicates were collected from each individual and ERCC spike-in control were added.
#' Sample were sequenced using the SMARTer protocol, using 5-bp UMIs. A total of 864 cells were sequenced.
#' Cells from second replicate of individual NA19098 were removed due to a pipetting error during in ERCC addition.
#' Additional cells were removed using an ad-hoc quality control procedure (see (\emph{Tung et al., 2017})).
#' Finally, only 564 cells were marked as high quality samples, resulting in:
#' \itemize{
#'   \item Tung_C1: 85 cells from first replicate of NA19098
#'   \item Tung_C2: 57 cells from third replicate of NA19098
#'   \item Tung_C3: 80 cells from first replicate of NA19101
#'   \item Tung_C4: 70 cells from second replicate of NA19101
#'   \item Tung_C5: 51 cells from third replicate of NA19101
#'   \item Tung_C6: 74 cells from first replicate of NA19239
#'   \item Tung_C7: 68 cells from second replicate of NA19239
#'   \item Tung_C8: 79 cells from third replicate of NA19239
#' }
#' @docType data
#'
#' @usage data(Tung_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Tung parameter preset
#' data(Tung_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Tung_like_dataset <- SPARSim_simulation(Tung_param_preset)
#'
#' # simulate a count table resembling only the last three experimental conditions
#' # i.e. Tung_C6, Tung_C7 and Tung_C8
#' Tung_like_dataset_2 <- SPARSim_simulation(Tung_param_preset[c("Tung_C6", "Tung_C7", "Tung_C8")])
#'
#' @references Tung, P.Y. et al. Batch effects and the effective design of single-cell gene expression studies. Sci Rep. 2017;7:39921.
"Tung_param_preset"






#' SPARSim simulation parameter preset from Camp dataset
#'
#' SPARSim simulation parameter preset estimated from Camp (\emph{Camp et al., 2015; GSE75140}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Camp count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 39913 genes
#'   \item 434 cells (across 6 different experimental conditions)
#'   \item average read/UMI count per cell of 981747
#'   \item sparsity ~84\%
#' }
#' Information about the original Camp dataset:
#' \itemize{
#'   \item Number of genes: 39913
#'   \item Number of cells: 434
#'   \item Number of different experimental conditions: 6 (Camp_C1, Camp_C2, Camp_C3, Camp_C4, Camp_C5, Camp_C6)
#'   \item Sparsity: ~84.00\%
#'   \item Average read/UMI count per cell: 981747
#'   \item Species: Human
#'   \item Cell types: Whole brain organoids, Microdissected Cerebral Organoids
#'   \item Protocol and platform: Fluidigm C1, SMARTer
#'   \item UMIs: No
#' }
#'
#' Additional details about the Camp dataset:
#'
#' scRNA-seq data were produced from human fetal neocortex and human cerebral organoids.
#' A total of 734 cells were sequenced. All samples were processed on the microfluidic Fluidigm C1 and contain 92 external RNA spike-ins.
#' Samples were sequenced using the SMARTer protocol.
#' Cerebral organoid data were generated from dissociated whole organoids derived from induced pluripotent stem cell line 409B2 (iPSC 409B2) at 33 days (40 cells), 35 days (68 cells), 37 days (71 cells), 41 days (74 cells), and 65 days (80 cells) after the start of embryoid body culture.
#' Cerebral organoid data were also generated from microdissected cortical-like regions from H9 embryonic stem cell derived organoids at 53 days (region 1, 48 cells; region 2, 48 cells) or from iPSC 409B2 organoids at 58 days (region 3, 43 cells; region 4, 36 cells).
#' Samples from Whole brain organoids at day 41 were removed due to an error during Splat parameter estimation, resulting in 434 cells:
#' \itemize{
#'   \item Camp_C1: 96 cells from Microdissected Cerebral Organoids at day 53
#'   \item Camp_C2: 79 cells from Microdissected Cerebral Organoids at day 58
#'   \item Camp_C3: 40 cells from Whole brain organoids at day 33
#'   \item Camp_C4: 68 cells from Whole brain organoids at day 35
#'   \item Camp_C5: 71 cells from Whole brain organoids at day 37
#'   \item Camp_C6: 80 cells from Whole brain organoids at day 65
#' }
#' @docType data
#'
#' @usage data(Camp_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Camp parameter preset
#' data(Camp_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Camp_like_dataset <- SPARSim_simulation(Camp_param_preset)
#'
#' # simulate a count table resembling only the Whole brain organoids
#' # i.e. Camp_C3, Camp_C4, Camp_C5 and Camp_C6
#' Camp_like_dataset_2 <- SPARSim_simulation(Camp_param_preset[c("Camp_C3","Camp_C4",
#'                                                              "Camp_C5","Camp_C6")])
#'
#' @references Camp, J.G., et al. Human cerebral organoids recapitulate gene expression programs of fetal neocortex development. Proc Natl Acad Sci U S A. 2015;112:15672–7.
"Camp_param_preset"


#' Camp dataset
#'
#' Count matrix from Camp study (\emph{Camp et al., 2015; GSE75140}):
#' \itemize{
#'   \item Number of genes: 39913
#'   \item Number of cells: 434
#'   \item Number of different experimental conditions: 6 (Camp_C1, Camp_C2, Camp_C3, Camp_C4, Camp_C5, Camp_C6)
#'   \item Sparsity: ~84.00\%
#'   \item Average read/UMI count per cell: 981747
#'   \item Species: Human
#'   \item Cell types: Whole brain organoids, Microdissected Cerebral Organoids
#'   \item Protocol and platform: Fluidigm C1, SMARTer
#'   \item UMIs: No
#' }
#'
#'
#' Additional details about the Camp dataset:
#'
#' scRNA-seq data were produced from human fetal neocortex and human cerebral organoids.
#' A total of 734 cells were sequenced. All samples were processed on the microfluidic Fluidigm C1 and contain 92 external RNA spike-ins.
#' Samples were sequenced using the SMARTer protocol.
#' Cerebral organoid data were generated from dissociated whole organoids derived from induced pluripotent stem cell line 409B2 (iPSC 409B2) at 33 days (40 cells), 35 days (68 cells), 37 days (71 cells), 41 days (74 cells), and 65 days (80 cells) after the start of embryoid body culture.
#' Cerebral organoid data were also generated from microdissected cortical-like regions from H9 embryonic stem cell derived organoids at 53 days (region 1, 48 cells; region 2, 48 cells) or from iPSC 409B2 organoids at 58 days (region 3, 43 cells; region 4, 36 cells).
#' Samples from Whole brain organoids at day 41 were removed due to an error during Splat parameter estimation, resulting in 434 cells:
#' \itemize{
#'   \item Camp_C1: 96 cells from Microdissected Cerebral Organoids at day 53
#'   \item Camp_C2: 79 cells from Microdissected Cerebral Organoids at day 58
#'   \item Camp_C3: 40 cells from Whole brain organoids at day 33
#'   \item Camp_C4: 68 cells from Whole brain organoids at day 35
#'   \item Camp_C5: 71 cells from Whole brain organoids at day 37
#'   \item Camp_C6: 80 cells from Whole brain organoids at day 65
#' }
#' @docType data
#'
#' @usage data(Camp_count_matrix)
#'
#' @format A data frame having genes on rows and cells on columns.
#'
#' @examples
#' # load Camp data count matrix
#' data(Camp_count_matrix)
#'
#' # estimate library size from Camp count table
#' Camp_lib_size <- SPARSim_estimate_library_size(Camp_count_matrix)
#'
#' @references Camp, J.G., et al. Human cerebral organoids recapitulate gene expression programs of fetal neocortex development. Proc Natl Acad Sci U S A. 2015;112:15672–7.
"Camp_count_matrix"




#' SPARSim simulation parameter preset from Engel dataset
#'
#' SPARSim simulation parameter preset estimated from Engel (\emph{Engel et al., 2016; GSE74596}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Engel count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 25865 genes
#'   \item 203 cells (across 4 different experimental conditions)
#'   \item average read/UMI count per cell of 4215547
#'   \item sparsity ~69\%
#' }
#' Information about the original Engel dataset:
#' \itemize{
#'   \item Number of genes: 25865
#'   \item Number of cells: 203
#'   \item Number of different experimental conditions: 4 (Engel_C1, Engel_C2, Engel_C3, Engel_C4)
#'   \item Sparsity: ~68.59\%
#'   \item Average read/UMI count per cell: 4215547
#'   \item Species: Mouse
#'   \item Cell types: Natural killer T cells
#'   \item Protocol and platform: Flow cytometry, Modified Smart-seq2
#'   \item UMIs: No
#' }
#'
#' Additional details about the Engel dataset:
#'
#' scRNA-seq data were produced from mouse Natural killer T cells.
#' NKT cells subtypes (NKT0, NKT1, NKT17 and NKT2) were isolated from thymuses and directly sorted by flow cytometry into lysis buffer (96 well plate single cell sort).
#' A modified Smart-seq2 protocol was used for sequencing.
#' A total of 203 cells were sequenced, resulting in
#' \itemize{
#'   \item Engel_C1: 45 NTK0 cells
#'   \item Engel_C2: 46 NTK1 cells
#'   \item Engel_C3: 44 NTK17 cells
#'   \item Engel_C4: 68 NTK2 cells
#' }
#' @docType data
#'
#' @usage data(Engel_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Engel parameter preset
#' data(Engel_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Engel_like_dataset <- SPARSim_simulation(Engel_param_preset)
#'
#' # simulate a count table resembling only NTK1 and NTK2 cells, i.e. Engel_C2 and Engel_C4
#' Engel_like_dataset_2 <- SPARSim_simulation(Engel_param_preset[c("Engel_C2", "Engel_C4")])
#'
#' @references Engel, I., et al., Innate-like functions of natural killer T cell subsets result from highly divergent gene programs. Nat Immunol. 2016;17:728–39.
"Engel_param_preset"






#' SPARSim simulation parameter preset from Chu dataset
#'
#' SPARSim simulation parameter preset estimated from Chu (\emph{Chu et al., 2016; GSE75748)}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Chu count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 17782 genes
#'   \item 758 cells (across 6 different experimental conditions)
#'   \item average read/UMI count per cell of 1335606
#'   \item sparsity ~51\%
#' }
#' Information about the original Chu dataset:
#' \itemize{
#'   \item Number of genes: 17782
#'   \item Number of cells: 758
#'   \item Number of different experimental conditions: 6 (Chu_C1, Chu_C2, Chu_C3, Chu_C4, Chu_C5, Chu_C6)
#'   \item Sparsity: ~51.30\%
#'   \item Average read/UMI count per cell: 1335606
#'   \item Species: Human
#'   \item Cell types: iPSC to Endoderm
#'   \item Protocol and platform: Fluidigm C1, SMARTer
#'   \item UMIs: No
#' }
#'
#' Additional details about the Chu dataset:
#'
#' scRNA-seq data were produced from a 6 step time series experiment of Human pluripotent stem cells (hPSCs) differentiation into definitive endoderm.
#' The six time points correspond to 0 h, 12 h, 24 h, 36 h, 72 h, and 96 h along the differentiation protocol from embryonic stem cell to definitive endoderm.
#' Single-cell capture and library preparations were performed using the Fluidigm C1 and a SMARTer protocol.
#' A total of 758 cells were sequenced, resulting in:
#' \itemize{
#'   \item Chu_C1: 92 cells at 0 h
#'   \item Chu_C2: 102 cells at 12 h
#'   \item Chu_C3: 66 cells at 24 h
#'   \item Chu_C4: 172 cells at 36 h
#'   \item Chu_C5: 138 cells at 72 h
#'   \item Chu_C6: 188 cells at 96 h
#' }
#' @docType data
#'
#' @usage data(Chu_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Chu parameter preset
#' data(Chu_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Chu_like_dataset <- SPARSim_simulation(Chu_param_preset)
#'
#' # simulate a count table resembling only cells at time points 0h and 96h, i.e. Chu_C1 and Chu_C6
#' Chu_like_dataset_2 <- SPARSim_simulation(Chu_param_preset[c("Chu_C1", "Chu_C6")])
#'
#' @references Chu, L. F., et al. Single-cell RNA-seq reveals novel regulators of human embryonic stem cell differentiation to definitive endoderm. Genome Biol, 2016, 17.1: 173.
"Chu_param_preset"






#' SPARSim simulation parameter preset from Horning dataset
#'
#' SPARSim simulation parameter preset estimated from Horning (\emph{Horning et al., 2018; GSE99795)}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Horning count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 37347 genes
#'   \item 144 cells (across 3 different experimental conditions)
#'   \item average read/UMI count per cell of 2572504
#'   \item sparsity ~39\%
#' }
#' Information about the original Horning dataset:
#' \itemize{
#'   \item Number of genes: 37347
#'   \item Number of cells: 144
#'   \item Number of different experimental conditions: 3 (Horning_C1, Horning_C2, Horning_C3)
#'   \item Sparsity: ~39.25\%
#'   \item Average read/UMI count per cell: 2572504
#'   \item Species: Human
#'   \item Cell types: LNCap (prostate cancer)
#'   \item Protocol and platform: Smart-seq2
#'   \item UMIs: No
#' }
#'
#' Additional details about the Horning dataset:
#'
#' scRNA-seq data were produced from LNCaP prostate cancer cells treated and untreated with androgen after cell cycle synchronization.
#' The three treatment groups (Group 1: cells 0h after synchronization and androgen deprivation; Group 2: cells 12h after synchronization and androgen deprivation; Group 3: cells 12h after synchronization and androgen addition) were processed using a SMART-seq2 protocol.
#' A total of 144 cells were sequenced, resulting in:
#' \itemize{
#'   \item Horning_C1: 48 cells of Group 1
#'   \item Horning_C2: 48 cells of Group 2
#'   \item Horning_C3: 48 cells of Group 3
#' }
#' @docType data
#'
#' @usage data(Horning_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Horning parameter preset
#' data(Horning_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Horning_like_dataset <- SPARSim_simulation(Horning_param_preset)
#'
#' # simulate a count table resembling only cells at time points 12h, i.e. Horning_C2 and Horning_C3
#' Horning_like_dataset_2 <- SPARSim_simulation(Horning_param_preset[c("Horning_C2", "Horning_C3")])
#'
#' @references Horning, A.M., et al. Single-Cell RNA-seq Reveals a Subpopulation of Prostate Cancer Cells with Enhanced Cell-Cycle-Related Transcription and Attenuated Androgen Response. Cancer Res 2018 Feb 15;78(4):853-864.
"Horning_param_preset"







#' SPARSim simulation parameter preset from Bacher dataset
#'
#' SPARSim simulation parameter preset estimated from Bacher (\emph{Bacher et al., 2017; GSE85917}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Bacher count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 17128 genes
#'   \item 366 cells (across 4 different experimental conditions)
#'   \item average read/UMI count per cell of 2504380
#'   \item sparsity ~39\%
#' }
#' Information about the original Bacher dataset:
#' \itemize{
#'   \item Number of genes: 17128
#'   \item Number of cells: 366
#'   \item Number of different experimental conditions: 4 (Bacher_C1, Bacher_C2, Bacher_C3, Bacher_C4)
#'   \item Sparsity: ~39.08\%
#'   \item Average read/UMI count per cell: 2504380
#'   \item Species: Human
#'   \item Cell types: Embryonic stem cells
#'   \item Protocol and platform: Fluidigm C1, SMARTer
#'   \item UMIs: No
#' }
#'
#' Additional details about the Bacher dataset:
#'
#' scRNA-seq data were produced from H9 and H1 hESCs. 91 hESCs were captured using a Fluidigm C1 system and processed using the SMARTer protocol. The identical single-cell indexed and fragmented cDNA were pooled at 24 cells per lane or at 96 cells per lane, resulting in approximately 4 million and 1 million mapped reads per cell in the two groups, respectively. The same process was then repeated for H1 hESCs. A total of 366 cells were sequenced, resulting in:
#' \itemize{
#'   \item Bacher_C1: 91 H9 hESCs (at 24 single cell cDNA libraries per lane)
#'   \item Bacher_C2: 92 H1 hESCs (at 24 single cell cDNA libraries per lane)
#'   \item Bacher_C3: 91 H9 hESCs (at 96 single cell cDNA libraries per lane)
#'   \item Bacher_C4: 92 H1 hESCs (at 96 single cell cDNA libraries per lane)
#' }
#' @docType data
#'
#' @usage data(Bacher_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Bacher parameter preset
#' data(Bacher_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Bacher_like_dataset <- SPARSim_simulation(Bacher_param_preset)
#'
#' # simulate a count table resembling only the first two experimental conditions
#' # i.e. Bacher_C1 and Bacher_C2
#' Bacher_like_dataset_2 <- SPARSim_simulation(Bacher_param_preset[c("Bacher_C1", "Bacher_C2")])
#'
#' @references Bacher, R., et al. SCnorm: robust normalization of single-cell RNA-seq data. Nat Methods, 2017, 14.6: 584.
"Bacher_param_preset"



#' SPARSim simulation parameter preset from Zheng dataset
#'
#' SPARSim simulation parameter preset estimated from Zheng (\emph{Zheng et al., 2017}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Zheng count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 19536 genes
#'   \item 3388 cells (across 4 different experimental conditions)
#'   \item average read/UMI count per cell of 15041
#'   \item sparsity ~83\%
#' }
#' Information about the original Zheng dataset:
#' \itemize{
#'   \item Number of genes: 19536
#'   \item Number of cells: 3388
#'   \item Number of different experimental conditions: 4 (Zheng_C1, Zheng_C2, Zheng_C3, Zheng_C4)
#'   \item Sparsity: ~82.56\%
#'   \item Average read/UMI count per cell: 15041
#'   \item Species: Human
#'   \item Cell types: Jurkat and 293T cells
#'   \item Protocol and platform: 10X Genomics
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the Zheng dataset:
#'
#' scRNA-seq data were produced from a sample where an equal number of 293T and Jurkat cells was mixed.
#' Cells were processed using the GemCode Single-Cell 3′ platform.
#' A total of 3388 cells were sequenced, resulting in:
#' \itemize{
#'   \item Zheng_C1: 1440 Jurkat cells
#'   \item Zheng_C2: 1718 293T cells
#'   \item Zheng_C3: 184 mixed cells, i.e. samples containing both Jurkat and 293T cells
#'   \item Zheng_C4: 46 unclassified cells
#' }
#' @docType data
#'
#' @usage data(Zheng_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Zheng parameter preset
#' data(Zheng_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Zheng_like_dataset <- SPARSim_simulation(Zheng_param_preset)
#'
#' # simulate a count table resembling only Jurkat and 293T cells, i.e. Zheng_C1 and Zheng_C2
#' Zheng_like_dataset_2 <- SPARSim_simulation(Zheng_param_preset[c("Zheng_C1", "Zheng_C2")])
#'
#' @references Zheng, G.H.Y., et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun. 2017; 8:14049
#' @references https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/jurkat:293t_50:50
"Zheng_param_preset"



#' SPARSim simulation parameter preset from Macosko dataset
#'
#' SPARSim simulation parameter preset estimated from Macosko (\emph{Macosko et al., 2015; GSE63472-GSM1626794}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Macosko count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 20230 genes
#'   \item 9000 cells
#'   \item average read/UMI count per cell of 1129
#'   \item sparsity ~97\%
#' }
#' Information about the original Macosko dataset:
#' \itemize{
#'   \item Number of genes: 20230
#'   \item Number of cells: 9000
#'   \item Number of different experimental conditions: -
#'   \item Sparsity: ~96.54\%
#'   \item Average read/UMI count per cell: 1129
#'   \item Species: Mouse
#'   \item Cell types: Retinal cells
#'   \item Protocol and platform: Drop-Seq
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the Macosko dataset:
#'
#' scRNA-seq data were produced from wild-type P14 mouse retinal cells.
#' Cells were processed using the Drop-Seq platform.
#' A total of 9000 cells were sequenced.
#'
#' @docType data
#'
#' @usage data(Macosko_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Macosko parameter preset
#' data(Macosko_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Macosko_like_dataset <- SPARSim_simulation(Macosko_param_preset)
#'
#' @references Makosco, E.Z., et al. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets. Cell. 2015:161(5):1202-1214
"Macosko_param_preset"





#' SPARSim simulation parameter preset from Saunders dataset
#'
#' SPARSim simulation parameter preset estimated from Saunders (\emph{Saunders et al., 2018; GSE116470}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Saunders count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 19940 genes
#'   \item 5688 cells
#'   \item average read/UMI count per cell of 2183
#'   \item sparsity ~94\%
#' }
#' Information about the original Saunders dataset:
#' \itemize{
#'   \item Number of genes: 19940
#'   \item Number of cells: 5688
#'   \item Number of different experimental conditions: -
#'   \item Sparsity: ~94.10\%
#'   \item Average read/UMI count per cell: 2183
#'   \item Species: Mouse
#'   \item Cell types: Polydendrocytes
#'   \item Protocol and platform: Drop-Seq
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the Saunders dataset:
#'
#' scRNA-seq data were produced from P60-70 C57BL/6 male mice.
#' Cells were processed using the Drop-Seq platform.
#' A total of 5688 cells were sequenced.
#'
#' @docType data
#'
#' @usage data(Saunders_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Saunders parameter preset
#' data(Saunders_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Saunders_like_dataset <- SPARSim_simulation(Saunders_param_preset)
#'
#' @references Saunders, A., et al. Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain. Cell. 2018:174(4) P1015-1030.E16
#' @references http://dropviz.org/
"Saunders_param_preset"





#' SPARSim simulation parameter preset from 10X Brain dataset
#'
#' SPARSim simulation parameter preset estimated from 10X Brain (\emph{https://support.10xgenomics.com/single-cell-gene-expression/datasets/}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to Brain_10X count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 21080 genes
#'   \item 11843 cells
#'   \item average read/UMI count per cell of 7601
#'   \item sparsity ~87\%
#' }
#' Information about the original 10X Brain dataset:
#' \itemize{
#'   \item Number of genes: 21080
#'   \item Number of cells: 11843
#'   \item Number of different experimental conditions: -
#'   \item Sparsity: ~87.37\%
#'   \item Average read/UMI count per cell: 7601
#'   \item Species: Mouse
#'   \item Cell types: Brain cells
#'   \item Protocol and platform: 10X Genomics
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the 10X Brain dataset:
#'
#' scRNA-seq data were produced from brain of an E18 mouse.
#' Cells were processed using the Chromium Single-Cell 3′ v3 platform.
#' A total of 11843 cells were sequenced.
#'
#' @docType data
#'
#' @usage data(Brain_10X_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load Brain_10X parameter preset
#' data(Brain_10X_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' Brain_10X_like_dataset <- SPARSim_simulation(Brain_10X_param_preset)
#'
#' @references https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/neuron_10k_v3
"Brain_10X_param_preset"





#' SPARSim simulation parameter preset from 10X T dataset
#'
#' SPARSim simulation parameter preset estimated from 10X T (\emph{https://support.10xgenomics.com/single-cell-gene-expression/datasets/}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to T_10X count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 20267 genes
#'   \item 8093 cells
#'   \item average read/UMI count per cell of 3880
#'   \item sparsity ~95\%
#' }
#' Information about the original 10X T dataset:
#' \itemize{
#'   \item Number of genes: 20267
#'   \item Number of cells: 8093
#'   \item Number of different experimental conditions: -
#'   \item Sparsity: ~94.59\%
#'   \item Average read/UMI count per cell: 3880
#'   \item Species: Human
#'   \item Cell types: Pan T cells
#'   \item Protocol and platform: 10X Genomics
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the 10X T dataset:
#'
#' scRNA-seq data were produced from pan T cells isolated from mononuclear cells of a healthy donor.
#' Cells were processed using the Chromium Single-Cell 3′ v2 platform.
#' A total of 8093 cells were sequenced.
#'
#' @docType data
#'
#' @usage data(T_10X_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load T_10X parameter preset
#' data(T_10X_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' T_10X_like_dataset <- SPARSim_simulation(T_10X_param_preset)
#'
#' @references https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_3k_4k_aggregate
"T_10X_param_preset"






#' SPARSim simulation parameter preset from 10X PBMC dataset
#'
#' SPARSim simulation parameter preset estimated from 10X PBMC (\emph{https://support.10xgenomics.com/single-cell-gene-expression/datasets/}) count table.
#' It allows simulating a count table having characteristics (i.e. number of genes, number of cells, gene intensity, gene variability, sparsity and library sizes) similar to PBMC_10X count matrix.
#' The simulation parameter will generate a count matrix having:
#' \itemize{
#'   \item 17220 genes
#'   \item 5419 cells
#'   \item average read/UMI count per cell of 3880
#'   \item sparsity ~96\%
#' }
#' Information about the original 10X PBMC dataset:
#' \itemize{
#'   \item Number of genes: 17220
#'   \item Number of cells: 5419
#'   \item Number of different experimental conditions: -
#'   \item Sparsity: ~95.66\%
#'   \item Average read/UMI count per cell: 3880
#'   \item Species: Human
#'   \item Cell types: Peripheral blood mononuclear cells
#'   \item Protocol and platform: 10X Genomics
#'   \item UMIs: Yes
#' }
#'
#' Additional details about the 10X PBMC dataset:
#'
#' scRNA-seq data were produced from peripheral blood mononuclear cells (PBMCs) from a healthy donor.
#' Cells were processed using the Chromium Single-Cell 3′ v1 platform.
#' A total of 5419 cells were sequenced.
#'
#' @docType data
#'
#' @usage data(PBMC_10X_param_preset)
#'
#' @format A list of one or more SPARSim simulation parameters (one of each experimental condition) as output of function \code{\link{SPARSim_estimate_parameter_from_data}}.
#' The single simulation parameter has the format described in \code{\link{SPARSim_create_simulation_parameter}}.
#'
#' @examples
#' # load PBMC_10X parameter preset
#' data(PBMC_10X_param_preset)
#'
#' # simulate a count table using the just loaded parameter preset
#' PBMC_10X_like_dataset <- SPARSim_simulation(PBMC_10X_param_preset)
#'
#' @references https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc6k
"PBMC_10X_param_preset"




