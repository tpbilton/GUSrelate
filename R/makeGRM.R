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
                    filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL),
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
  if(GUSbase::isValue(indsubset, type="pos_integer", minv=1, maxv=RAobj$.__enclos_env__$private$nInd))
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

  ## initalize the GRM object
  GRMobj <- GRM$new(RAobj, ploid, indsubset)

  ## Compute allele frequencies if required
  if(isTRUE(est$MAF)){
    GRMobj$.__enclos_env__$private$p_est(nThreads=nThreads)
  }
  ## Compute p-value for HWE if required
  if(isTRUE(est$HWE)){
    GRMobj$.__enclos_env__$private$HWE_est(nThreads=nThreads)
  }
  ## Compute proportion of missing data
  GRMobj$.__enclos_env__$private$miss <- apply(GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt,
                2, function(x) sum(x==0)/length(x))

  ## compute the GRM
  GRMobj$computeGRM(method=method, ep=ep, filter=filter)

  ## return the GRM object
  return(GRMobj)
}


