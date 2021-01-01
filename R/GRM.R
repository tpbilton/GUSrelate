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
#' A GRM object is created from the \code{\link{makeGRM}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (or methods)
#' for analyzing the data. Information in an GRM are specific to constructing a GRM
#' @usage
#' GRMobj <- makeGRM()
#' @format NULL
#' @author Timothy P. Bilton
#' @seealso \code{\link{makeGRM}}
#' @name GRM
#' @export
### R6 class for creating a data format for genomic relationship matrices

GRM <- R6::R6Class("GRM",
                     inherit = GUSbase::RA,
                     public = list(
                       initialize = function(GRMobj, ploid, indsubset, saminfo){
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
                         private$samInfo   <- saminfo
                       },
                       computeGRM = function(name, method="VanRaden", ep=0, snpsubset=NULL, filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL),...){
                         ## do some checks
                         if(!is.vector(name) || !is.character(name) || length(name) != 1)
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         if(is.null(filter$MAF)) filter$MAF <- 0
                         else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
                           stop("Minor allele frequency filter is invalid")
                         if(is.null(filter$MISS)) filter$MISS <- 1
                         else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
                           stop("Proportion of missing data filter is invalid")
                         if(is.null(filter$PVALUE)) filter$PVALUE <- 0
                         else if( length(filter$PVALUE) != 1 || !is.numeric(filter$PVALUE) || filter$PVALUE<0 || filter$PVALUE>1 )
                           stop("P-value for Hardy-Weinberg equilibrium filter is invalid.")
                         match.arg(method, choices = c("VanRaden","WG"))
                         if(!is.null(snpsubset) & GUSbase::checkVector(snpsubset, type="pos_integer", maxv=private$nSnps))
                           stop("SNP subset index is invalid")

                         ## compute the subset of SNPs
                         subset <- rep(TRUE, private$nSnps)
                         if(!is.null(snpsubset))
                           subset[-snpsubset] <- FALSE
                         subset[which(private$maf < filter$MAF)] <- FALSE
                         subset[which(private$miss > filter$MISS)] <- FALSE
                         
                         if(!is.null(private$pvalue))
                           subset[which(private$pvalue < filter$PVALUE)] <- FALSE
                         else
                           warning("Ignoring p-value filter as p-values from the Hardy-Weinberg equilibrium test have not been computed. Use the `$HWEtest` function to compute p-values")
                    
                         ## compute the GRM
                         GRMmat <- GUSrelate::computeGRM(ref=private$ref, alt=private$alt, ploid=private$ploid, snpsubset=which(subset), method=method,
                                               phat=private$pfreq, ep=ep, ...)
                         fil_val <- list(MAF=max(private$filter$MAF,filter$MAF),
                                         MISS=min(private$filter$MISS,filter$MISS),
                                         PVALUE=filter$PVALUE)
                         GRMlist <- list(GRM=GRMmat, method=method, filter=fil_val, indID=private$indID, snpsubset=snpsubset, ep=ep, freq=private$pfreq[snpsubset])
                         ## work which GRM to save
                         private$GRM[[name]] <- GRMlist
                         return(invisible())
                       },
                       #p_est = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                       #   temp <- GUSbase::p_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                       #                                      indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                       #   private$pfreq <- temp$p
                       #},
                       HWEtest = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                         if(length(unique(private$ploid[indsubset])!=1))
                           stop("HWE test is not yet implemented for mixed ploidy populations.")
                         pest <- GUSbase::p_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         gest <- GUSbase::g_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         private$pvalue <- 1-pchisq(-2*(pest$loglik - gest$loglik), df=private$ploid-1)
                       },
                       ############# Plots for the GRM
                       PCA = function(name, npc=3, group1=NULL, group2=NULL, group.hover=NULL, interactive=FALSE){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1)
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         else if(!(name %in% names(private$GRM)))
                           stop("GRM not found. Check the name of the GRM group.")
                         if(!is.null(group1) && (!is.vector(group1) || !is.character(group1) || length(group1) != 1 || !(group1 %in% names(private$ginfo))))
                           stop("Information for 'group1' variable not found in the Sample information. Check the name of the sample information variable")
                         if(!is.null(group2) && (!is.vector(group2) || !is.character(group2) || length(group2) != 1 || !(group2 %in% names(private$ginfo))))
                           stop("Information for 'group1' variable not found in the Sample information. Check the name of the Sample information variable")
                         if(!is.null(group.hover) && (!is.vector(group.hover) || !is.character(group.hover) || any(!(group.hover %in% names(private$ginfo)))))
                           stop("Information for 'group1' variable not found in the Sample information. Check the name of the Sample information variable")
                         GRMinfo <- private$GRM[[name]]
                         GRM <- GRMinfo$GRM
                         ## compte PCs components
                         PC <- svd(GRM - matrix(colMeans(GRM), nrow = nrow(GRM), ncol = ncol(GRM), byrow = TRUE), nu = npc)
                         eval <- sign(PC$d) * PC$d^2/sum(PC$d^2)
                         PC$x <- PC$u %*% diag(PC$d[1:npc],nrow=npc) # nrow to get correct behaviour when npc=1
                         ## Sort the groups
                         if(!is.null(group1)) g1 <- private$ginfo[[group1]]
                         else g1 <- NULL
                         if(!is.null(group2)) g2 <- private$ginfo[[group2]]
                         else g2 <- NULL
                         ## Produce the plots
                         if(interactive){
                           if(!is.null(group.hover)) hover.info <- apply(sapply(group.hover,function(x) paste0(x,": ",private$ginfo[[x]]),simplify = TRUE),1,paste0,collapse="<br>")
                           else hover.info <- NULL
                           temp_p <- plotly::plot_ly(y=PC$x[, 2],x=PC$x[, 1], type="scatter", mode="markers",
                                             hoverinfo="text", text=hover.info, width=640, height=640,
                                             marker=list(size=6), color=g1, symbol=g2) %>%
                             plotly::layout(xaxis=list(title="Principal component 1",zeroline=FALSE), yaxis=list(title="Principal component 2",zeroline=FALSE))
                           temp_p
                         } else{
                           df <- data.frame(x=PC$x[, 1], y=PC$x[, 2])
                           ggplot2::ggplot(df, ggplot2::aes(x=x,y=y, col=g1, shape=g2)) + ggplot2::geom_point() + ggplot2::theme_bw() +
                             ggplot2::ylab("Principal component 2") + ggplot2::xlab("Principal component 1")
                         }
                       },
                       #### add group information
                       addSampleInfo = function(name, info){
                         info <- as.character(info)
                         if(!is.vector(name) || !is.character(name) || length(name) != 1 || (name %in% names(private$ginfo)))
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         if(!is.vector(info) || length(info) != private$nInd)
                           stop(paste0("Argument 'info' needs to be a character vector of length ",private$nInd,"."))
                         private$ginfo[[name]] <- info
                       },
                       deleteSampleInfo = function(name){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1 || !(name %in% names(private$ginfo)))
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         else if(!(name %in% names(private$ginfo)))
                           stop("Information not found. Check the name of the group.")
                         private$ginfo[[name]] <- NULL
                       },
                       writeGRM = function(name, filename, IDvar=NULL){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1 || !(name %in% names(private$GRM)))
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         
                         GRM <- private$GRM[[name]]$GRM
                         if(!is.null(IDvar)){
                           match.arg(IDvar, names(private$sam))
                           colnames(GRM) <- rowname(GRM) <- private$samInfo[[IDvar]]
                         } else colnames(GRM) <- rowname(GRM) <- private$indID
                         
                         write.table(GRM, file = filename)
                       }
                     ),
                     private = list(
                       ratio  = NULL,
                       ploid  = NULL,
                       pfreq  = NULL,
                       pvalue = NULL,
                       miss   = NULL,
                       maf    = NULL,
                       snpdepth = NULL,
                       samdepth = NULL,
                       ep     = NULL,
                       gfreq  = NULL,
                       GRM    = NULL,
                       filter = NULL,
                       ginfo  = NULL,
                       samInfo = NULL
                     )
)
