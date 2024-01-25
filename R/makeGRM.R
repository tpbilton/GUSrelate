##########################################################################
# Genotyping Uncertainty with Sequencing data and RELATEdness (GUSrelate)
# Copyright 2019-2021 Timothy P. Bilton
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

#' Make an GRM object
#' 
#' Create a genomic relationship matix (GRM) object from an RA object, perform standard filtering and 
#' compute statistics required for constructing GRMs.
#' 
#' This function converts an RA object into a GRM object. A GRM object is
#' a R6 type object that contains RA data, various statistics related 
#' to GRM analyses and functions (methods) for analysing GRMs.
#' 
#' The sample information as specified in the \code{samfile} argument should be a csv file with the first column
#' giving the ID of the sample (and must match the IDs in the RA object supplied in the \code{RAobj} argument) and
#' the second column giving the ploidy level of each individual. Additional columns can then be added to the that
#' gives more information about the sample (e.g., cultivar, location, population). For example,
#' \tabular{lll}{
#' ID \tab Ploidy \tab Group \cr
#' HE109 \tab 2 \tab Wild \cr
#' PE202 \tab 4 \tab Wild \cr
#' PE243 \tab 4 \tab Domesticated   
#' }
#' In this example, there are three individuals, the first (HE109) is a diploid and belongs to the Wild group,
#' the second individual (PE202) is a tetraploid that also belongs to the Wild group and the third individual
#' (PE243) is also a tetraploid but belongs to the Domesticated group. 
#' Note that the names for the first two columns must be "ID" and "Ploid" respectively, but any names can
#' be used for the remaining columns but it is recommended meaning name without spaced. Remeber that the first
#' two are required, any extra columns are optional but can be used later. 
#' 
#' The filtering criteria currently implemented are:
#' \itemize{
#' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is \code{NULL})}
#' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
#'  is greater than the threshold value (default is \code{NULL}).}
#' \item{Bin size for SNP selection (\code{BIN}): }{SNPs are binned together if the distance (in base pairs) between them is less than the threshold value (default is 100).
#' One SNP is then randomly selected from each bin and retained for final analysis. This filtering is to ensure that there is only one SNP on each sequence read.}
#' \item{Maximum averge SNP depth (\code{MAXDEPTH}): }{SNPs with an average read depth above the threshold value are discarded (default is 500).}
#' }
#' If a filtering criteria is set to \code{NULL}, then no filtering in regard to
#' that threshold is applied. 
#' 
#' @param RAobj Object of class RA created via the \code{\link[GUSbase]{readRA}} function.
#' @param samfile Character string giving the name of the file that contains the sample information of the population. 
#' See below for details.
#' @param method A character string specifying whether the VanRaden (\code{'VanRaden'}) based estimator or 
#' the Weir-Goudet (\code{'WG'}) estimator is used to construct the GRM.
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
#' 
#' @author Timothy P. Bilton
#' @return An R6 object of class GRM.
#' @export makeGRM

#### Make an unstructured population
makeGRM <- function(RAobj, samfile, filter=list(MAF=NULL, MISS=NULL, BIN=100, MAXDEPTH=500)){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(!is.vector(samfile) || !is.character(samfile) || length(samfile) != 1)
    stop("Argument `saminfo` is invalid. Must be a character of length 1.")
  else if(!file.exists(samfile))
    stop("File for sample information is not found")
  
  ## check filter values
  if(is.null(filter$MAF)) filter$MAF <- 0
  else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
    stop("Minor allele frequency filter is invalid")
  if(is.null(filter$MISS)) filter$MISS <- 1
  else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
    stop("Proportion of missing data filter is invalid")
  if(is.null(filter$MAXDETH)) filter$MAXDEPTH <- 500
  else if( length(filter$MAXDEPTH) != 1 || GUSbase::checkVector(filter$MAXDEPTH, type="pos_integer"))
    stop("Maximum SNP depth filter is invalid.")
  if(is.null(filter$BIN)) filter$BIN <- 0
  else if( length(filter$BIN) != 1 || GUSbase::checkVector(filter$BIN, type="pos_integer"))
    stop("Binning distance filter is invalid.")
  
  ## read in the same information file:
  saminfo = as.data.frame(data.table::fread(samfile, header=T))
  
  ## Checks on the sample info file
  if(nrow(saminfo)<2)
    stop("There is less than two samples in the sample file. Cannot proceed.")
  # Sample IDs and ploidy values are present sample info?
  if(all((names(saminfo) != "ID")))
    stop("No column for sample IDs in sample file")
  if(all((names(saminfo) != "Ploidy")))
    stop("No column for ploidy values in sample file")
  # Ploidy values are valid?
  if(any(sapply(saminfo$Ploidy, function(x) GUSbase::checkVector(x, type="pos_integer", minv=1))))
    stop("Ploidy values in sample file are invalid. These need to be non-missing positive integers")
  # Sample IDs in RA data?
  indIndex = saminfo$ID %in% RAobj$.__enclos_env__$private$indID
  if(any(!indIndex))
    stop(paste0("Sample ID(s) not present in RA dataset:\n",
                paste(saminfo$ID[!indIndex],collapse = "\n")))
  
  ## initalize the GRM object
  indIndex_ra = match(saminfo$ID,RAobj$.__enclos_env__$private$indID)
  GRMobj <- GRM$new(RAobj, saminfo$Ploidy, indIndex_ra, saminfo[indIndex,])
  
  ## compute depth matrix
  depth = GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt

  ## Compute proportion of missing data
  miss <- apply(depth, 2, function(x) sum(x==0)/length(x))
  ## MAF based on observed reads
  pfreq <- colMeans(GRMobj$.__enclos_env__$private$ref/depth, na.rm=T)
  maf <- pmin(pfreq,1-pfreq)

  ## Compute mean SNP depth
  snpdepth <- colMeans(depth)
  samdepth <- rowMeans(depth)

  ## filter SNPs
  indx <- which(miss < filter$MISS & maf > filter$MAF & snpdepth < filter$MAXDEPTH)
  chrom <- GRMobj$.__enclos_env__$private$chrom[indx]
  pos <- GRMobj$.__enclos_env__$private$pos[indx]
  
  ## determine the distance between adajcent SNPs
  if(filter$BIN > 0){
    oneSNP <- rep(FALSE,length(chrom))
    temp = sapply(unique(chrom), function(x){
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
    },USE.NAMES = F)
    oneSNP[unlist(temp)] <- TRUE
    oneSNP <- which(oneSNP)
  } else
    oneSNP <- which(rep(TRUE, length(chrom)))
  
  ## Update GRMobj
  GRMobj$.__enclos_env__$private$miss <- miss[indx[oneSNP]]
  GRMobj$.__enclos_env__$private$pfreq <- pfreq[indx[oneSNP]]
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

  ## update summary information
  temp <- GRMobj$.__enclos_env__$private$ref + GRMobj$.__enclos_env__$private$alt
  summaryInfo = list(
    header="Data Summary:\n",
    file=paste0("Data file:\t\t",GRMobj$.__enclos_env__$private$rafile,"\n"),
    meandepth=paste0("Mean Depth:\t\t", round(mean(temp),2),"\n"),
    callrate=paste0("Mean Call Rate:\t",round(sum(temp!=0)/length(temp),2),"\n"),
    num="Number of...\n",
    samples=paste0("  Samples:\t\t",GRMobj$.__enclos_env__$private$nInd,"\n"),
    snps=paste0("  SNPs:\t\t",GRMobj$.__enclos_env__$private$nSnps,"\n"),
    reads=paste0("  Reads:\t\t",sum(temp),"\n"))
  GRMobj$.__enclos_env__$private$summaryInfo = summaryInfo

  ## Print summary information
  GRMobj$print()
  
  ## return the GRM object
  return(GRMobj)
}


