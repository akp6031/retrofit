Welcome to RETROFIT, an R package for learning cell-type composition in spatial transcriptomic (ST) data.

Lack of single-cell resolution in spatial transcriptomic (ST) technologies produces gene expression from a mixture of potentially heterogeneous cells at each spot. We have developed a reference-free deconvolution method (RETROFIT), to decompose cell-type mixtures in ST data in an unsupervised manner, without relying on an external single-cell expression reference. 

This package can be used to obtain cell-type compositions for the spots in ST data generated by any platform such as 10x Genomics Visium data or Slide-seq data. Our model takes the input of an ST data matrix for deconvolution and a reference of either a list of known markers or single-cell RNAseq data for annotation. 

##
![Figure1](https://user-images.githubusercontent.com/90921267/220766755-daea9d4b-4ac0-4dd3-978c-7e71b31bc36e.png) <br />

## Installation

Please follow these steps to install `retrofit` from github:

``` r
install.packages("devtools") 
devtools::install_github("qunhualilab/retrofit")
```
Install `retrofit` from [Bioconductor](https://bioconductor.org/packages/devel/bioc/html/retrofit.html).

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')
BiocManager::install("retrofit")
```

## Example

Use the following code to apply RETROFIT to your data.

``` r
library(retrofit)
## load built in data
utils::data(testSimulationData)
x           <- testSimulationData$extra5_x
sc_ref      <- testSimulationData$sc_ref
marker_ref  <- testSimulationData$marker_ref

## decompose 
res         <- retrofit::decompose(x, L=16, iterations=100, verbose=TRUE)
W           <- res$w
H           <- res$h
TH          <- res$th

## annotate with correlations
res         <- retrofit::annotateWithCorrelations(sc_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells

## annotate with markers
res         <- retrofit::annotateWithMarkers(marker_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells	  
```

More details on how to use RETROFIT can be found in the vignettes below:
- [Simulation Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/SimulationVignette.Rmd) : This vignette uses sample data simulated by us to demonstrate the use of the package and how to interpret the generated results. This vignette can be used to get started with RETROFIT and understand its usage.

- [Colon Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/ColonVignette.Rmd) : This vignette utilizes real data from Human Colon tissue generated using the 10x Genomics Visium platform in this [paper](https://www.sciencedirect.com/science/article/pii/S009286742031686X). Herein, we demonstrate that our method is effective in identifying biologically relevant spatial patterns using ST tissues. This is a slightly more advanced vignette as it not only uses RETROFIT to deconvolve ST data but also analyzes its results.
