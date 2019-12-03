##########################################################################
# Genotyping Uncertainty with Sequencing data and RELATEdness (GUSrelate)
# Copyright 2019 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################

#' Make an genomic relationship matix (GRM) object
#' 
#' Create a GRM object from an RA object, perform standard filtering and 
#' compute statistics specific required for constructing GRMs.
#' 
#' This function converts an RA object into a GRM object. A GRM object is
#' a R6 type obtain that contains RA data, various statistics related for 
#' to GRM analyses and functions (methods) for analysing GRMs.
#' 
#' The filtering criteria currently implemented are:
#' \iterize{
#' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is \code{NULL})}
#' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
#'  is greater than the threshold value (default is \code{NULL}).}
#' \item{Bin size for SNP selection (\code{BIN}): }{SNPs are binned together if the distance (in base pairs) between them is less than the threshold value (default is 100).
#' One SNP is then randomly selected from each bin and retained for final analysis. This filtering is to ensure that there is only one SNP on each sequence read.}
#' \item{Hardy Weinberg equilibrium (HWE) test P-value (\code{PVALUE}): }{SNPs are discarded if the p-value from a HWE test is smaller than the threshold (default is \code{NULL}).
#'  This filters out SNPs where the segregation type has been inferred wrong.}
#' \item{Maximum averge SNP depth (\code{MAXDEPTH}}{SNPs with an average read depth above the threshold value are discarded (default is 500).}
#' }
#' If a filtering criteria is set to \code{NULL}, then no filtering in regard to
#' that threshold is applied. 
#' 
#' @param RAobj Object of class RA created via the \code{\link[GUSbase]{readRA}} function.
#' @param name Name of GRM matrix to be constructed.
#' @param ploid Integer value giving the ploidy level of the individuals in the population.
#' @param method Character string specifying whether the VanRaden (\code{'VanRaden'}) based estimator or 
#' the Weir-Goudet (\code{'WG'}) estimator is used to construct the GRM.
#' @param indsubset
#' @param nThreads Integer vector specifying the number of threads to use in the OpenMP parallelization 
#' used in the estimation of allele frequencies when \code{est=list(mafEst=TRUE)} of in the 
#' estimation of the p-value from a HWE test when \code{est=list(HWE=TRUE)} 
#' @param ep   
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
#' @param est Named list specifying whether to estimate the minor frequency (\code{MAF=TRUE})
#' or to compute the pvalue for the Hardy Weinberg equilibrium (\code{HWE=TRUE}) test 
#' for each SNP using the method by \insertCite{li2011bioinform}{GUSrelate}
#' 
#' @references 
#' \insertRef{li2011bioinform}{GUSrelate}
#' 
#' @author Timothy P. Bilton
#' @return An R6 object of class GRM.
#' @export makeGRM

#### Make an unstructured population
makeGRM <- function(RAobj, name, ploid=2, method="VanRaden", indsubset=NULL, nThreads=1, ep=0,
                    filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL, MAXDEPTH=500, BIN=100),
                    est=list(MAF=TRUE, HWE=TRUE)){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(!is.vector(name) || !is.character(name) || length(name) != 1)
    stop("Argument `name` is invalid.")
  if(!is.vector(ploid) || !is.numeric(ploid) || length(ploid) != 1 || round(ploid/2) != ploid/2)
    stop("Argument for ploid level is invalid.")
  if(!is.numeric(nThreads) || length(nThreads) != 1 || nThreads < 0 || round(nThreads) != nThreads)
    stop("Argument for the number of cores for the MPI parallelization is invalid")
  if(is.null(indsubset))
    indsubset <- 1:RAobj$.__enclos_env__$private$nInd
  if(GUSbase::checkVector(indsubset, type="pos_integer", minv=1, maxv=RAobj$.__enclos_env__$private$nInd))
    stop("Argument for the indices of the individuals is invalid.")
  if(is.null(filter$MAF)) filter$MAF <- 0
  else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
    stop("Minor allele frequency filter is invalid")
  if(is.null(filter$MISS)) filter$MISS <- 1
  else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
    stop("Proportion of missing data filter is invalid")
  if(is.null(filter$PVALUE)) filter$PVALUE <- 0
  else if( length(filter$PVALUE) != 1 || !is.numeric(filter$PVALUE) || filter$PVALUE<0 || filter$PVALUE>1 )
    stop("P-value for Hardy-Weinberg equilibrium filter is invalid.")
  if(is.null(filter$PVALUE)) filter$MAXDEPTH <- 500
  else if( length(filter$PVALUE) != 1 || GUSbase::checkVector(filter$PVALUE, type="pos_integer"))
    stop("Maximum SNP depth filter is invalid.")

  ## initalize the GRM object
  GRMobj <- GRM$new(RAobj, ploid, indsubset)

  ## Compute proportion of missing data
  miss <- apply(GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt,
                  2, function(x) sum(x==0)/length(x))
  ## MAF based on observed reads
  pfreq <- colMeans(GRMobj$.__enclos_env__$private$ref/(GRMobj$.__enclos_env__$private$ref +
                                                        GRMobj$.__enclos_env__$private$alt), na.rm=T)
  maf <- pmin(pfreq,1-pfreq)

  ## Compute mean SNP depth
  snpdepth <- colMeans(GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt)
  samdepth <- rowMeans(GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt)

  indx <- which(miss < filter$MISS & maf > filter$MAF & snpdepth < filter$MAXDEPTH)

  chrom <- GRMobj$.__enclos_env__$private$chrom[indx]
  pos <- GRMobj$.__enclos_env__$private$pos[indx]
  ## determine the distance between adajcent SNPs
  if(filter$BIN > 0){
    oneSNP <- rep(FALSE,length(chrom))
    oneSNP[unlist(sapply(unique(chrom), function(x){
      ind <- which(chrom == x)
      if(length(ind > 1)){
        g1_diff <- diff(pos[ind])
        SNP_bin <- c(0,cumsum(g1_diff > filter$BIN)) + 1
        set.seed(58473+as.numeric(which(x==chrom))[1])
        keepPos <- sapply(unique(SNP_bin), function(y) {
          ind2 <- which(SNP_bin == y)
          if(length(ind2) > 1)
            return(sample(ind2,size=1))
          else if(length(ind2) == 1)
            return(ind2)
        })
        return(ind[keepPos])
      }
      else return(ind)
    },USE.NAMES = F ))] <- TRUE
    oneSNP <- which(oneSNP)
  } else
    oneSNP <- which(rep(TRUE, length(chrom)))
  GRMobj$.__enclos_env__$private$miss <- miss[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$maf <- maf[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$snpdepth <- snpdepth[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$samdepth <- samdepth
  GRMobj$.__enclos_env__$private$ref <- GRMobj$.__enclos_env__$private$ref[,indx[oneSNP]]
  GRMobj$.__enclos_env__$private$alt <- GRMobj$.__enclos_env__$private$alt[,indx[oneSNP]]
  GRMobj$.__enclos_env__$private$chrom <- GRMobj$.__enclos_env__$private$chrom[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$pos <- GRMobj$.__enclos_env__$private$pos[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$SNP_Names <- GRMobj$.__enclos_env__$private$SNP_Names[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$nSnps <- length(indx[oneSNP])
  GRMobj$.__enclos_env__$private$filter <- filter

  ## Compute allele frequencies if required
  if(isTRUE(est$MAF)){
    GRMobj$p_est(nThreads=nThreads)
  } else{
    GRMobj$.__enclos_env__$private$pfreq <- pfreq[indx[oneSNP]]
  }
  ## Compute p-value for HWE if required
  if(isTRUE(est$HWE)){
    GRMobj$HWE_est(nThreads=nThreads)
  }

  ## compute the GRM
  GRMobj$computeGRM(name=name, method=method, ep=ep, filter=filter)

  ## return the GRM object
  return(GRMobj)
}


