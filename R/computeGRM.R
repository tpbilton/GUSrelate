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
#' @param ep 
#' @param snpsubset 
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
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
  else if(GUSbase::checkVector(ploid, type = "pos_integer", minv=1) || length(ploid) != nInd)
    stop("Input vector ploidy level is invalid.")
  
  depth <- ref + alt
  ratio <- ref/depth
  if(any(!(length(ep) %in% c(1,nSnps)))){
    warning("Vector for sequencing error parameter is not equal to 1 or the number of SNPs.\nSetting to zero.")
    ep = 0
  }
  if(method=="VanRaden" && (is.null(phat) || length(phat) != nSnps))
    stop("Allele frequency vector is not supplied or not equal to the number of SNPs")
  ## subset the data if required
  if(!is.null(snpsubset)){
    if(length(ep) == nSnps)
      ep <- ep[snpsubset]
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
    if(length(ep) > 1) ep = matrix(ep[snpsubset], nrow=nInd, ncol=nsnpsub, byrow=T)
    else               ep = matrix(rep(ep,nsnpsub), nrow=nInd, ncol=nsnpsub, byrow=T)

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
    return(GRM)
  }
  else if(method == "WG"){
    ## currently not implemented
    ratio[which(depth < 2)] <- NA
    snpsubset <- which(!((colMeans(ratio, na.rm=T) == ploid) | (colMeans(ratio, na.rm=T) == 0)))
    ratio <- ratio[, snpsubset]
    nSnps <- length(snpsubset)
    na_indx <- !is.na(ratio)
    na_mat <- tcrossprod(na_indx,na_indx)
    ratio[which(is.na(ratio))] <- ploid/2

    drat <- 1/depth[,snpsubset]
    drat[which(depth[,snpsubset] < 2)] <- 0

    epMat <- matrix(ep, nrow=nInd, ncol=nSnps)
    epMat[which(depth[,snpsubset] < 2)] <- 0

    mat <- 1/2 + 2/(ploid^2)*(tcrossprod((ratio-ploid/2)/sqrt(1-4*epMat*(1-epMat))) -
                                  tcrossprod(sqrt(((ploid^2/4)*na_indx)/(1-4*epMat*(1-epMat)))) + tcrossprod(sqrt(ploid^2/4*na_indx)) +
                                  tcrossprod(sqrt(ploid^2*epMat*(1-epMat)/(1-4*epMat*(1-epMat)))))/na_mat
    diag(mat) <- 1/2 + 2/(ploid^2)*rowSums( (((ratio-ploid/2)^2 - (ploid)^2/4)/(1-drat) + ploid^2*epMat*(1-epMat))/(1-4*epMat*(1-epMat)) + ploid^2/4)/diag(na_mat)

    mat_sum <- (sum(mat) - sum(diag(mat)))/(nrow(mat)*(nrow(mat)-1))
    GRM <- ploid*(mat - mat_sum)/(1-mat_sum)
    return(GRM)
  }
}
