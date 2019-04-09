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
#' @return An R6 object of class GRM.
#' @export makeGRM

#### Make an unstructured population
makeGRM <- function(RAobj, ploid=2, method="VanRaden", indsubset=NULL, nThreads=1, ep=0,
                    filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL, MAXDEPTH=500, BIN=100),
                    est=list(MAF=TRUE, HWE=TRUE)){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
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

  ## Compute allele frequencies if required
  if(isTRUE(est$MAF)){
    GRMobj$p_est(nThreads=nThreads)
  } else{
    GRMobj$.__enclos_env__$private$pfreq <- pfreq
  }
  ## Compute p-value for HWE if required
  if(isTRUE(est$HWE)){
    GRMobj$HWE_est(nThreads=nThreads)
  }

  ## compute the GRM
  GRMobj$computeGRM(method=method, ep=ep, filter=filter)

  ## return the GRM object
  return(GRMobj)
}


