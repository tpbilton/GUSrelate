##########################################################################
# Genotyping Uncertainty with Sequencing data for RELATEdness (GUSrelate)
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
#' GRM method: Construct a genomic relationship matrix (GRM)
#'
#' Method for constructing a genomic relationship matrix (GRM) for a diploid or autopolyploid population.
#'
#' @param name A character string giving the name of the GRM analysis
#' @param method A character string specifying whether the VanRaden (\code{'VanRaden'}) based estimator or
#' the Weir-Goudet (\code{'WG'}) estimator is used to construct the GRM.
#' @param ploid A positive integer vector specifying the ploidy level of each individual. If only a single value is given,
#' then the ploidy level is assumed equal in all the individuals.
#' @param ep Sequencing error value. Can be a single number (sequencing error rate same for all genotypes),
#' a vector equal to the number of SNPs (SNP specific sequencing error rate) or a matrix the same diminsion as the data (SNP and
#' individual specific sequencing error rate).
#' @param snpsubset Vector of indices of SNPs that should be retained in the analysis. Useful for subsetting the SNPs before
#' any filtering is applied.
#' @param phat Vector of allele frequencies to use in the constuction of the GRM.
#'
#' @seealso \code{\link{GRM}}
#' @author Timothy P. Bilton
#' @export computeGRM

computeGRM <- function(ref, alt, ploid, snpsubset=NULL, method="VanRaden", phat, ep=0, ...){

  if(!exists("thres")) thres = 0.001

  ## do some checks
  match.arg(method, choices = c("VanRaden","WG"))
  if(all(dim(ref) != dim(alt)))
    stop("The matrix of reference alleles has different dimensions to matrix of alternate alleles")
  else{
    nSnps = ncol(ref)
    nInd  = nrow(ref)
  }
  if(length(ploid) == 1)
    ploid = rep(ploid, nInd)
  else if(any(sapply(ploid, function(x) GUSbase::checkVector(x, type="pos_integer", minv=1))) || length(ploid) != nInd)
    stop("Input vector ploidy level is invalid.")

  depth <- ref + alt
  ratio <- ref/depth
  ## Check sequencing error argument
  if(is.vector(ep) || is.null(ep)){
    if(any(!(length(ep) %in% c(1,nSnps)))){
      warning("Vector for sequencing error parameter given but is not equal to 1 or the number of SNPs.\nSetting to zero.")
      ep = 0
    }
    if(length(ep) == 1) ep = matrix(ep, nrow=nInd, ncol=nSnps, byrow=T)
    else if(length(ep) == nSnps){
      ep = matrix(rep(ep,nSnps), nrow=nInd, ncol=nSnps, byrow=T)
    }
  } else{
    if(nrow(ep) != nInd)
      stop("Number of rows in the sequencing error matrix supplied to argument 'ep' does not equal the number of individuals in the dataset. Check input")
    if(ncol(ep) != nSnps)
      stop("Number of columns in the sequencing error matrix supplied to argument 'ep' does not equal the number of SNPs in the dataset. Check input")
    if(any(ep < 0 | ep > 1))
      stop("Argument 'ep' contains values that not between 0 and 1. check input")
    if(any(is.na(ep))){
      message("Argument 'ep' contains missing values. Setting to zero")
      ep[which(is.na(ep))] = 0
    }
  }

  if(method=="VanRaden" && (is.null(phat) || length(phat) != nSnps))
    stop("Allele frequency vector is not supplied or not equal to the number of SNPs")
  ## subset the data if required
  if(!is.null(snpsubset)){
    ep <- ep[,snpsubset]
    ## compute depth and dosage matrix
    depth <- depth[,snpsubset]
    ratio <- ratio[,snpsubset]
    if(method=="VanRaden") phat <- phat[snpsubset]
  }

  #################
  ## Compute GRM ##
  #################

  if(method == "VanRaden"){
    snpsubset <- which(phat < 1-thres & phat > thres)
    depth = depth[,snpsubset]
    nsnpsub <- length(snpsubset)
    phat <- phat[snpsubset]
    ep = ep[,snpsubset]

    ## Compute the adjusted GRM
    genon0 <- ratio[,snpsubset] - rep.int(phat, rep(nInd, nsnpsub))
    genon0[depth<1] = 0
    P0 <- matrix(phat,nrow=nInd,ncol=nsnpsub,byrow=T)
    P1 <- 1-P0
    P0[depth<1] = 0
    P1[depth<1] = 0
    ep[depth<1] = 0
    div0 <- tcrossprod(P0,P1)
    gammaMat = 1-4*ep*(1-ep)
    nuMat = 1-4*P0*P1
    GRM <- (tcrossprod(genon0/sqrt(gammaMat)) - tcrossprod(sqrt((ep^2*nuMat/gammaMat))) )/div0

    ## Adjustment for self-diagonals
    genon0[depth==1] = 0
    gammaMat[depth==1] = 1
    nuMat[depth==1] = 1
    P0[depth==1] = 0
    P1[depth==1] = 0
    ep[depth==1] = 0
    deltaMat = 1 - 1/depth
    deltaMat[depth < 2] = 1
    adj <- P0*P1*(1 - gammaMat*deltaMat) + ep*(nuMat-(1-ep)*deltaMat)
    diag(GRM) <- rowSums(  (genon0^2 - adj)/(gammaMat*deltaMat))/rowSums(P0*P1)
    GRM = sqrt(tcrossprod(ploid))*GRM
  }
  else if(method == "WG"){
    warning("The Weir-Goudet estimator for polyploids has yet to be tested. Use with caution.")
    ratio[which(depth < 1)] = NA
    snpsubset = which(!((colMeans(ratio, na.rm=T) == 1) | (colMeans(ratio, na.rm=T) == 0)))
    ratio = ratio[, snpsubset]
    depth = depth[,snpsubset]
    nSnps = length(snpsubset)
    na_indx = !is.na(ratio)
    na_mat = tcrossprod(na_indx,na_indx)
    ratio[which(is.na(ratio))] = 1/2

    ## Compute adjustment for off-diagonals
    ep = ep[,snpsubset]
    kappaMat =  ep*(1 - ep)
    gammaMat = 1 - 4*kappaMat
    kappaMat[which(depth < 1)] = 0
    gammaMat[which(depth < 1)] = 1
    mat <- 1/2 + 2*(tcrossprod((ratio - 1/2)/sqrt(gammaMat)) -
                                  tcrossprod(sqrt(((1/4)*na_indx)/gammaMat)) + tcrossprod(sqrt(1/4*na_indx)) +
                                  tcrossprod(sqrt(kappaMat/gammaMat)))/na_mat

    ## Compute adjustment to diagonals
    deltaMat = 1 - 1/depth
    deltaMat[which(depth < 2)] = 1
    kappaMat[which(depth == 1)] = 0
    gammaMat[which(depth == 1)] = 1
    na_diag = rowSums(depth > 1)
    ratio[which(depth < 2)] = 0
    #diag(mat) = 1/2 + rowSums(2*((((ratio-1/2)^2 - 1/4)/deltaMat + kappaMat)/gammaMat + 1/4))/diag(na_mat)
    diag(mat) = 1 + rowSums( (2*(ratio^2 - ratio)/deltaMat + 2*kappaMat)/gammaMat)/na_diag

    ## Construct GRM
    mat_sum = (sum(mat) - sum(diag(mat)))/(nrow(mat)*(nrow(mat)-1))
    GRM = sqrt(tcrossprod(ploid))*(mat - mat_sum)/(1-mat_sum)
  }
  return(GRM)
}
