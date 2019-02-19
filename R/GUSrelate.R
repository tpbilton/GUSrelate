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

#' GRM object
#'
#' Class for storing RA data and associated functions for analyzing unstructured populations (e.g.,
#' populations with no known structure).
#'
#' An US object is created from the \code{\link{makeUS}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (or methods)
#' for analyzing the data. Information in an US are specific to unrelated populations (or
#' populations with no known relationships).
#' @usage
#' USobj <- makeUS()
#' @format NULL
#' @author Timothy P. Bilton
#' @seealso \code{\link{makeUR}}
#' @name GRM
#' @export
### R6 class for creating a data format for unstructured populations

GRM <- R6Class("GRM",
                     inherit = RA,
                     public = list(
                       initialize = function(GRMobj, ploid, indsubset){
                         private$ref       <- GRMobj$.__enclos_env__$private$ref[indsubset,]
                         private$alt       <- GRMobj$.__enclos_env__$private$alt[indsubset,]
                         private$chrom     <- GRMobj$.__enclos_env__$private$chrom
                         private$pos       <- GRMobj$.__enclos_env__$private$pos
                         private$SNP_Names <- GRMobj$.__enclos_env__$private$SNP_Names
                         private$indID     <- GRMobj$.__enclos_env__$private$indID[indsubset]
                         private$nSnps     <- GRMobj$.__enclos_env__$private$nSnps
                         private$nInd      <- length(indsubset)
                         private$gform     <- GRMobj$.__enclos_env__$private$gform
                         private$AFrq      <- GRMobj$.__enclos_env__$private$AFrq
                         private$infilename<- GRMobj$.__enclos_env__$private$infilename
                         private$ploid     <- as.integer(ploid)
                       },
                       computeGRM = function(method="VanRaden", ep=0, filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL), ...){
                         ## do some checks
                         if(is.null(filter$MAF)) filter$MAF <- 0
                         else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
                           stop("Minor allele frequency filter is invalid")
                         if(is.null(filter$MISS)) filter$MISS <- 1
                         else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
                           stop("Proportion of missing data filter is invalid")
                         if(is.null(filter$PVALUE)) filter$PVALUE <- 0
                         else if( length(filter$PVALUE) != 1 || !is.numeric(filter$PVALUE) || filter$PVALUE<0 || filter$PVALUE>1 )
                           stop("P-value for Hardy-Weinberg equilibrium filter is invalid.")
                         if(length(method) != 1 || !is.character(method) || any(!(method %in% c("VanRaden","WG"))))
                           stop("Method argument must be either 'VanRaden' or 'WG'")


                         ## compute the subset of SNPs
                         snpsubset <- rep(TRUE, private$nSnps)
                         if(!is.null(private$pfreq)){
                           maf <- pmin(private$pfreq,1-private$pfreq)
                           snpsubset[which(maf < filter$MAF)] <- FALSE
                         }
                         if(!is.null(private$pvalue))
                           snpsubset[which(private$pvalue < filter$PVALUE)] <- FALSE
                         if(!is.null(private$miss))
                           snpsubset[which(private$miss > filter$MISS)] <- FALSE
                         ## compute the GRM
                         GRMmat <- GUSrelate::computeGRM(private$ref, private$alt, private$ploid, which(snpsubset), method,
                                               private$pfreq, ep=ep, ...)
                         ## work which GRM to save
                         if(method == "VanRaden") private$GRM_VR <- GRMmat
                         else if(method == "WG")  private$GRM_WG <- GRMmat
                         return(invisible())
                       },
                       p_est = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                         temp <- GUSbase::p_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                                                            indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         private$pfreq <- temp$p
                       },
                       HWE_est = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                         pest <- GUSbase::p_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         gest <- GUSbase::g_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         private$pvalue <- 1-pchisq(-2*(pest$loglik - gest$loglik), df=private$ploid-1)
                       }
                     ),
                     private = list(
                       ratio  = NULL,
                       ploid  = NULL,
                       pfreq  = NULL,
                       pvalue = NULL,
                       miss   = NULL,
                       ep     = NULL,
                       gfreq  = NULL,
                       GRM_VR = NULL,
                       GRM_WG = NULL
                     )
)
