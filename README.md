# GUSrelate

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Genotyping Uncertainty with Sequencing data and RELATEdness.

An R package for constructing genomic relationship matrices (GRMs) for autopolyploid species using high-throughput sequencing data.

### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages("devtools")
devtools::install_github("tpbilton/GUSbase")
devtools::install_github("tpbilton/GUSrelate")
```

Note: Some of the functions are coded in C and therefore an appropriate C compiler is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a compiler.

### Example:

We give a simple example to illustrate how to use GUSrelate using the simulated dataset 
found in the `GUSbase` package. The data can be loaded and an RA object created using the following code
```
library(GUSrelate)
ra <- VCFtoRA(simDS()$vcf)
radata <- readRA(ra)
```

A genomic relationship matrix (GRM) can be constructed using the `makeGRM` function
```
grm <- makeGRM(RAobj=radata, name="VR", ploid=2, method="VanRaden", ep=0,
               filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL, BIN=NULL), est=list(MAF=F, HWE=F))
```
The arguments of the `makeGRM` function are:

* `RAobj`: An RA object created from the `makeRA` function
* `ploid`: An integer value giving the ploid level of individuals in the population
* `name`: A character string giving the name of the GRM analysis
* `method`: A character string specifying whether the VanRaden (`'VanRaden'`) based estimator or  the Weir-Goudet `'WG'` estimator is used to construct the GRM
* `ep`: A numeric vector specifying the sequencing error rate for each SNP or a numeric value specifying the overall sequencing error rate
* `filter`: 

Another GRM can be constructed with different parameters using the `makeGRM` function. We construct the GRM using the Weir-Goudet based estimator as follows:
```
grm$computeGRM(name="WG", method="WG")
```

A principal component analysis (PCA) plot of the GRM can be produced using the function `$plotPCA`:
```
grm$PCA(name="test")
```


### Funding:
The initial development of this package was partially funded by the Ministry of Business, Innovation and Employment via its funding of the “Genomics for Production & Security in a Biological Economy” programme (Contract ID C10X1306).
