# SPARSim

SPARSim is an R tool for the simulation of single cell RNA-seq (scRNA-seq) count table.

To cite this package please use  ***"SPARSim single cell: a count data simulator for scRNA-seq data", Giacomo Baruzzo, Ilaria Patuzzi, Barbara Di Camillo, Bioinformatics, Volume 36, Issue 5, March 2020, Pages 1468–1475,*** https://doi.org/10.1093/bioinformatics/btz752



## Applications of SPARSim in scientific research

SPARSim has been ranked among the top single cell transcriptomics data simulators:
- **SPARSim** was extended using _simAdaptor_, a tool that adds spatial variables to single-cell simulations. The combined SPARSim+simAdaptor approach ranked among the top spatial transcriptomics simulators in a comprehensive benchmark study involving 13 simulation methods. [_“Multi-task benchmarking of spatially resolved gene expression simulation models”_, **Genome Biology**, 2025](https://doi.org/10.1186/s13059-025-03505-w)
- **SPARSim** was ranked among the top-performing simulators in an extensive benchmark of 12 scRNA-seq simulation tools. [_“A benchmark study of simulation methods for single-cell RNA sequencing data”_, **Nature Communications**, 2021](https://doi.org/10.1038/s41467-021-27130-w)


Also, SPARSim has been widely adopted by the research community for benchmarking and evaluating single cell transcriptomics data analysis tools. Below is a non-exhaustive list of peer-reviewed publications that have utilized SPARSim in various contexts:

- **SPARSim** data were used to test _scEVE_, a novel ensemble clustering algorithm for scRNA-seq data. [_“scEVE: a single-cell RNA-seq ensemble clustering algorithm capitalizing on the differences of predictions between multiple clustering methods”_, **NAR Genomics & Bioinformatics**, 2025](https://doi.org/10.1093/nargab/lqaf073)
- **SPARSim**-simulated data were used to evaluate _SciGeneX_, a novel R package to generate co-expression gene modules in scRNA-seq and spatial transcriptomics. [_“SciGeneX: enhancing transcriptional analysis through gene module detection in single-cell and spatial transcriptomics data”_, **NAR Genomics & Bioinformatics**, 2025](https://doi.org/10.1093/nargab/lqaf043)
- **SPARSim**-simulated data were used to evaluate _SPARROW_, a computational framework for integrating scRNA-seq with spatial transcriptomics. [_“SPARROW reveals microenvironment-zone-specific cell states in healthy and diseased tissues”_, **Cell Systems**, 2025](https://doi.org/10.1016/j.cels.2025.101235)
- **SPARSim** is employed in _MOSim_ as the simulator for scRNA-seq data within a broader multilayer regulatory network framework, supporting both bulk (RNA-seq, ATAC-seq, miRNA-seq, ChIP-seq, and Methyl-seq) and single-cell datasets (scRNA-seq and scATAC-seq). [_“MOSim: bulk and single-cell multilayer regulatory network simulator”_, **Briefings in Bioinformatics**, 2025](https://doi.org/10.1093/bib/bbaf110)
- **SPARSim** serves as the scRNA-seq count data simulator in _AsaruSim_, a tool designed to generate synthetic single-cell and spatial RNA-seq Nanopore long reads datasets. [_“AsaruSim: a single-cell and spatial RNA-Seq Nanopore long-reads simulation workflow”_, **Bioinformatics**, 2025](https://doi.org/10.1093/bioinformatics/btaf087)
- **SPARSim**-simulated data supported the development of _Mcadet_, a novel feature selection framework for scRNA-seq data. [_“Mcadet: A feature selection method for fine-resolution single-cell RNA-seq data based on multiple correspondence analysis and community detection”_, **PLOS Computational Biology**, 2024](https://doi.org/10.1371/journal.pcbi.1012560)
- **SPARSim** data were used to test _Piccolo_, a feature selection methods for scRNA-seq data. [_“Feature selection followed by a novel residuals-based normalization that includes variance stabilization simplifies and improves single-cell gene expression analysis”_, **BMC Bioinformatics**, 2024](https://doi.org/10.1186/s12859-024-05872-w)
- **SPARSim**-simulated data were used to evaluate _jrSiCKLSNMF_, a clustering technique for single-cell multimodal omics data. [_“Clustering single-cell multimodal omics data with jrSiCKLSNMF”_, **Frontiers in Genetics**, 2023](https://doi.org/10.3389/fgene.2023.1179439)
- **SPARSim** was used to test _scGMM-VGAE_, a hybrid statistical and deep learning clustering methods for scRNA-seq data. [_“scGMM-VGAE: a Gaussian mixture model-based variational graph autoencoder algorithm for clustering single-cell RNA-seq data”_, **Machine Learning: Science and Technology**, 2023](https://doi.org/10.1088/2632-2153/acd7c3)
- **SPARSim** was used to validate _CONGAS+_, a computational method for inferring copy number variations from scRNA-seq and scATAC-seq. [_“A Bayesian method to infer copy number clones from single-cell RNA and ATAC sequencing”_, **PLOS Computational Biology**, 2023](https://doi.org/10.1371/journal.pcbi.1011557)
- **SPARSim**-simulated data were used to evaluate _UIPBC_, a novel clustering method for scRNA-seq data that requires no user input. [_“UIPBC: An effective clustering for scRNA-seq data analysis without user input”_, **Knowledge-Based Systems**, 2022](https://doi.org/10.1016/j.knosys.2022.108767)
- **SPARSim** data were used to test the _Polar Gini Curve_ method for identifying spatial gene expression patterns from scRNA-seq data. [_“Polar Gini Curve: A Technique to Discover Gene Expression Spatial Patterns from Single-Cell RNA-Seq Data”_, **Genomics, Proteomics & Bioinformatics**, 2021](https://doi.org/10.1016/j.gpb.2020.09.006)
- **SPARSim** was used to test _UICPC_, a clustering method for scRNA-seq data based on graph centrality measures. [_“UICPC: Centrality-based clustering for scRNA-seq data analysis without user input”_, **Computers in Biology and Medicine**, 2021](https://doi.org/10.1016/j.compbiomed.2021.104820)


## Installation

#### Install from GitLab

SPARSim is available on GitLab at https://gitlab.com/sysbiobig/sparsim

To install SPARSim from GitLab, please use the following commands:
```r
library(devtools)
install_gitlab("sysbiobig/sparsim", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```
The above commands would install SPARSim, the required dependencies and SPARSim vignette. 

To install SPARSim without its vignette, please use the following commands:
```r
library(devtools)
install_gitlab("sysbiobig/sparsim")
```



#### Install from source package

SPARSim R package can be downloaded at http://sysbiobig.dei.unipd.it/?q=SPARSim

It requires packages *RCpp*, *Matrix*, *scran* and *edgeR* to work.

To install from source, please use the following command:
```r
install.packages("SPARSim_0.9.5.tar.gz", repos = NULL, type = "source")
```

## Getting started

Please browse SPARSim vignettes using the following commands:

```r
library(SPARSim)
browseVignettes("SPARSim")
```

