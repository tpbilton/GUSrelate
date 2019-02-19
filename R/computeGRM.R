##########################################################################
# Genotyping Uncertainty with Sequencing data for RELATEdness (GUSrelate)
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
#' @export computeGRM

computeGRM <- function(ref, alt, ploid, snpsubset=NULL, method="VanRaden", phat, ep=0, ...){

  if(!exists("thres")) thres = 0.001

  ## do some checks
  if(length(method) != 1 || !is.character(method) || any(!(method %in% c("VanRaden","WG"))))
    stop("Method argument must be either 'VanRaden' or 'WG'")
  if(all(dim(ref) != dim(alt)))
    stop("The matrix of reference alleles has different dimensions to matrix of alternate alleles")
  else{
    nSnps = ncol(ref)
    nInd  = nrow(ref)
  }
  depth <- ref + alt
  ratio <- (ploid*ref/depth)
  if(any(!(length(ep) %in% c(1,nSnps)))){
    warning("Vector for sequecning error parameter is not equal to 1 or the number of SNPs.\nSetting to zero.")
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
    nsnpsub <- length(snpsubset)
    phat <- phat[snpsubset]
    if(length(ep) > 1) ep = matrix(ep[snpsubset], nrow=nInd, ncol=nsnpsub, byrow=T)
    else               ep = matrix(rep(ep,nsnpsub), nrow=nInd, ncol=nsnpsub, byrow=T)

    ## Compute the adjusted GRM
    genon0 <- ratio[,snpsubset] - ploid*rep.int(phat, rep(nInd, nsnpsub))
    genon0[is.na(ratio[,snpsubset])] <- 0
    genon0[depth[,snpsubset]<2] <- 0
    P0 <- matrix(phat,nrow=nInd,ncol=nsnpsub,byrow=T)
    P1 <- 1-P0
    P0[depth[,snpsubset]<2] <- 0
    P1[depth[,snpsubset]<2] <- 0
    ep[depth[,snpsubset]<2] <- 0
    div0 <- ploid*tcrossprod(P0,P1)
    GRM <- (tcrossprod(genon0/sqrt(1-4*ep*(1-ep))) - tcrossprod(sqrt(((ploid*ep)^2*(1-4*P0*P1)/(1-4*ep*(1-ep))))) )/div0
    depth.temp <- depth
    depth.temp[which(depth < 2)] <- 0
    depth.temp <- depth.temp[,snpsubset]
    depth.temp <- 1/depth.temp
    depth.temp2 <- depth.temp
    depth.temp[is.infinite(depth.temp)] <- 1
    depth.temp2[is.infinite(depth.temp2)] <- 0
    adj <- ploid^2*(P0*P1*(depth.temp + 4*ep*(1-ep)*(1-depth.temp)) +
                          ep*(ep+(1-ep)*depth.temp - 4*P0*P1))
    diag(GRM) <- rowSums(  (genon0^2 - adj)/((1-depth.temp2)*(1-4*ep*(1-ep))))/diag(div0)
    return(GRM)
  }
  else if(method == "WG"){
    snpsubset <- which(!((colMeans(ratio, na.rm=T) == ploid) | (colMeans(ratio, na.rm=T) == 0)))
    ratio <- ratio[, snpsubset]
    nSnps <- length(snpsubset)
    na_mat <- !is.na(ratio)
    na_mat <- tcrossprod(na_mat,na_mat)
    ratio[which(is.na(ratio))] <- ploid/2

    mat <- 1/2 + 2/ploid^2*tcrossprod(ratio - ploid/2)/na_mat
    mat_sum <- (sum(mat) - sum(diag(mat)))/(nrow(mat)*(nrow(mat)-1))
    GRM <- ploid*(mat - mat_sum)/(1-mat_sum)
    return(GRM)
  }
}
