# GUSrelate

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html) [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftpbilton%2FGUSrelate&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

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

The next step is to create a metadata file containing ploidy information of each sample (in case the ploidy level is not consistent in the population) and any other phenotype/grouping information that might be used later to explore the GRM. An example file is provide within the `GUSbase` package and the file location can be obtained using
```
(sim_metaInfo = simDS()$meta)
```

The next step is to create an GRM object.
```
grm <- makeGRM(RAobj=radata, samfile=sim_metaInfo, filter=list(MAF=NULL, MISS=NULL, BIN=100, MAXDEPTH=500))
```
The arguments of the `makeGRM` function are:

* `RAobj`: An RA object created from the `readRA` function.
* `samfile`: A character string of the path to the file containing the sample IDs and ploidy levels as described above.
* `filter`: list containing the filtering criteria to remove SNPs:
  * Minor allele frequency (`MAF`): SNPs are discarded if their MAF is less than the threshold (default is `NULL`) 
  * Proportion of missing data (`MISS`): SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype) is greater than the threshold value (default is `NULL`).
  * Bin size for SNP selection (`BIN`): SNPs are binned together if the distance (in base pairs) between them is less than the threshold value (default is 100). One SNP is then randomly selected from each bin and retained for final analysis. This filtering is to ensure that there is only one SNP on each sequence read.
  - Maximum average SNP depth (`MAXDEPTH`): SNPs with an average read depth above the threshold value are discarded (default is 500).

An additional filtering criteria based on Hardy Weinberg Disequilibrium (HWE) can be performed but requires performing the HWE test for each SNP as follows:
```
grm$HWEtest(nThreads = 3)
```
Comments on the `$HWEtest` function:

* Computing the HWE test can take some time. The `nThreads` argument specifies the number of threads to use in the calculation and can be used to speed up the processing.
* The HWE test has yet to be implemented for populations with multiple ploidy levels.

A GRM can be computed using the `$computeGRM` function as follows:
```
grm$computeGRM(name = "GRM_VR", method="VanRaden", ep=0, snpsubset=NULL, filter=list(MAF=NULL, MISS=NULL, PVALUE=0.01))
```
The arguments of the `$computeGRM` function are:

* `name`: A character string giving the name of the GRM.
* `method`: The method used to construct the GRM. Currently only `"VanRaden"` is implemented.
* `ep`: An sequencing error rate (based on a per read basis) to use when computing the GRM (default is 0). 
* `snpsubset`: An vector of indices of the SNPs that are to be retained in the analysis. Note: the SNPs must also pass the filtering criteria set in the `filter` argument. If `snpsubset=NULL`, then all the SNPs are retained prior to filtering.
* `filter`: A named list of thresholds for various filtering criteria. The GRM is constructed based on this filtering criteria and when a threshold is set to `NULL`, the filtering criterion is not applied. The different criteria are
  * Minor allele frequency (`MAF`): SNPs are discarded if their MAF is less than the threshold (default is `NULL`) 
  * Proportion of missing data (`MISS`): SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype) is greater than the threshold value (default is `NULL`).
  * Hardy-Weinberg equilibrium (HWE) test P-value (`PVALUE`): SNPs are discarded if the p-value from a HWE test is smaller than the threshold (default is `NULL`). This filters out badly behaving SNPs. Requires that `$HWEtest` function has been run first as described above. The HWE test implemented is the method described by Li (2009).

Another GRM can be constructed with different parameters using the `makeGRM` function. For example:
```
grm$computeGRM(name = "VR_filt", method="VanRaden", ep=0, snpsubset=NULL, filter=list(MAF=0.05, MISS=0.4, PVALUE=0.01))
```

A principal component analysis (PCA) plot of the GRM can be produced using the function `$plotPCA`:
```
## PCA of the GRM constructed using the WG method
grm$PCA(name = "VR_filt", colour=NULL, shape=NULL) 
```
The arguments of the `$PCA` function are:
* `name`: A character string giving the name of the GRM.
* `colour`: Column name of the file supplied to the `samfile` argument to use a grouping variable for colour on the PCA plot.
* `shape`: Column name of the file supplied to the `samfile` argument to use a grouping variable for shape of points on the PCA plot.

Lastly, the GRM created can be extracted from the GRM object using
```
grm$extractGRM(name = "VR_filt", IDvar=NULL)
```
or exported to a file using the function
```
grm$writeGRM(name = "VR_filt", filename = "VR_filt.csv", IDvar=NULL)
```
where the arguments are:

* `name`: A character string giving the name of the GRM.
* `filename`: The name of the file including the file extension (writes to csv format).
* `IDvar`: Which column in the file supplied to the `samfile` argument in the `makeGRM` function to use as IDs for the rows/columns of the GRM.

### Future development:

This package is still under development and additional methods and functions for analyzing GRMs is intended to be added at a later date.

### Funding:
The initial development of this package was partially funded by the Ministry of Business, Innovation and Employment via its funding of the “Genomics for Production & Security in a Biological Economy” programme (Contract ID C10X1306).
