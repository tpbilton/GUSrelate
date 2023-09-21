##########################################################################
# Genotyping Uncertainty with Sequencing data and RELATEdness (GUSrelate)
# Copyright 2019-2023 Timothy P. Bilton
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
#' A genomic relationship matrix (GRM) object is created from the \code{\link{makeGRM}} function and contains RA data,
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
                       ########### Print function:
                       print = function(what = NULL, ...){
                         cat(unlist(private$summaryInfo))

                         if(!is.null(private$GRM)){
                           cat("GRMs in object:\n\n")
                           listOfGRM = names(private$GRM)
                           junk = sapply(listOfGRM, function(x) {
                             cat(private$GRM[[x]]$sumInfo, sep="")
                             return(invisible())
                           })
                         }
                       },
                       ########### Function to construct a genomic relationship matrix (GRM)
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
                         ## compute summary information:
                         summaryInfo = c("Name: ", name, "\n  Number of individuals: ",nrow(GRMmat), "\n  Number of SNPs: ",sum(subset),
                                           "\n  Mean Depth: ", round(mean(private$ref[, which(subset)] + private$alt[, which(subset)]),2),
                                           "\n  Mean Self-relatedness: ", round(mean(diag(GRMmat)),2),"\n  Filtering (retained SNPs):",
                                           "\n    MAF: > ",fil_val$MAF, "\n    % Missing: < ", fil_val$MISS*100, "%",
                                           "\n    p-value: > ", filter$PVALUE, "\n")
                         ## Same GRM and associated information to GRM object
                         GRMlist <- list(GRM=GRMmat, method=method, filter=fil_val, indID=private$indID, snpsubset=which(subset), ep=ep, freq=private$pfreq[snpsubset], sumInfo = summaryInfo)
                         ## work which GRM to save
                         private$GRM[[name]] <- GRMlist
                         ## print out summary information
                         cat("GRM computed:\n", summaryInfo, sep="")
                         return(invisible(NULL))
                       },
                       #p_est = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                       #   temp <- GUSbase::p_est_em(private$ref, private$alt, private$ploid, snpsubset=snpsubset,
                       #                                      indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                       #   private$pfreq <- temp$p
                       #},
                       HWEtest = function(snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
                         if(is.null(indsubset)) indsubset = 1:private$nInd
                         else if(GUSbase::checkVector(indsubset, type = "pos_integer", minv = 0, maxv = private$nInd))
                           stop("Index for individuals is invalid")
                         ploidy = unique(private$ploid[indsubset])
                         if(length(ploidy) != 1)
                           stop("HWE test is not yet implemented for mixed-ploidy populations.")
                         pest <- GUSbase::p_est_em(private$ref, private$alt, ploidy, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         gest <- GUSbase::g_est_em(private$ref, private$alt, ploidy, snpsubset=snpsubset,
                                                   indsubset=indsubset, nThreads=nThreads, para=para, EMpara=EMpara)
                         private$pvalue <- 1-pchisq(-2*(pest$loglik - gest$loglik), df=ploidy-1)
                       },
                       ############# Delete a computed GRM:
                       removeGRM = function(name){
                         if(!is.vector(name) || !is.character(name))
                           stop("Argument 'name' needs to be a character vector")
                         if(any(!(name %in% names(private$GRM))))
                           stop("At least one GRM cannot be found: Check the names")
                         private$GRM[name] = NULL
                         if(length(private$GRM) == 1)
                           private$GRM = NULL
                       },
                       ############# Plot relatedness estimates
                       plotGRM = function(type=c("histrogram","scatter"),values=c("both","self-relatedness","relatednes"), xaxis=NULL,
                                          interactive=FALSE){
                         stop("yet to be implemeted")
                       },
                       ############ Compare different GRMs:
                       CompareGRM = function(){
                         stop("yet to be implemeted")
                       },
                       ############# Plots for the GRM
                       PCA = function(name, npc=3, colour=NULL, shape=NULL, group.hover=NULL, interactive=FALSE){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1)
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         else if(!(name %in% names(private$GRM)))
                           stop("GRM not found. Check the name of the GRM group.")
                         if(!is.null(colour) && (!is.vector(colour) || !is.character(colour) || length(colour) != 1 || !(colour %in% colnames(private$samInfo))))
                           stop("Information for 'colour' variable not found in the Sample information. Check the name of the sample information variable")
                         if(!is.null(shape) && (!is.vector(shape) || !is.character(shape) || length(shape) != 1 || !(shape %in% names(private$samInfo))))
                           stop("Information for 'shape' variable not found in the Sample information. Check the name of the Sample information variable")
                         if(!is.null(group.hover) && (!is.vector(group.hover) || !is.character(group.hover) || any(!(group.hover %in% names(private$samInfo)))))
                           stop("Information for 'group.hover' variable not found in the Sample information. Check the name of the Sample information variable")
                         GRMinfo <- private$GRM[[name]]
                         GRM <- GRMinfo$GRM
                         ## compte PCs components
                         PC <- svd(GRM - matrix(colMeans(GRM), nrow = nrow(GRM), ncol = ncol(GRM), byrow = TRUE), nu = npc)
                         eval <- sign(PC$d) * PC$d^2/sum(PC$d^2)
                         PC$x <- PC$u %*% diag(PC$d[1:npc],nrow=npc) # nrow to get correct behaviour when npc=1
                         ## Sort the groups
                         if(!is.null(colour)) g1 <- private$samInfo[[colour]]
                         else g1 <- NULL
                         if(!is.null(shape)) g2 <- private$samInfo[[shape]]
                         else g2 <- NULL
                         ## Produce the plots
                         lab1 = paste0("Principal component 1 (",round(eval[1]*100,2),"%)")
                         lab2 = paste0("Principal component 2 (",round(eval[2]*100,2),"%)")
                         if(interactive){
                           if(!is.null(group.hover)) hover.info <- apply(sapply(group.hover,function(x) paste0(x,": ",private$samInfo[[x]]),simplify = TRUE),1,paste0,collapse="<br>")
                           else hover.info <- NULL
                           temp_p <- plotly::plot_ly(y=PC$x[, 2],x=PC$x[, 1], type="scatter", mode="markers",
                                             hoverinfo="text", text=hover.info, width=640, height=640,
                                             marker=list(size=6), color=g1, symbol=g2) %>%
                             plotly::layout(xaxis=list(title=lab1,zeroline=FALSE), yaxis=list(title=lab2,zeroline=FALSE))
                           temp_p
                         } else{
                           df <- data.frame(x=PC$x[, 1], y=PC$x[, 2])
                           p = ggplot2::ggplot(df, ggplot2::aes(x=x,y=y, col=g1, shape=g2)) + ggplot2::geom_point() + ggplot2::theme_bw() +
                             ggplot2::ylab(lab2) + ggplot2::xlab(lab1)
                           if(!is.null(colour))
                             p = p + guides(color=guide_legend(title=colour))
                           if(!is.null(shape))
                             p = p + guides(shape=guide_legend(title=shape))
                           p
                         }
                       },
                       #### add group information
                       addSampleInfo = function(samfile){
                         if(!is.vector(samfile) || !is.character(samfile) || length(samfile) != 1)
                           stop("Argument `samfile` is invalid. Must be a character of length 1.")
                         else if(!file.exists(samfile))
                           stop("File for sample information is not found")

                         ## read in the same information file:
                         saminfo = as.data.frame(data.table::fread(samfile, header=T))

                         # Checks on the sample file
                         if(all((names(saminfo) != "ID")))
                           stop("No column for sample IDs in sample file")
                         if(any(colnames(samfile) == "Ploidy"))
                           stop("Ploidy column in sample file. Please remove")
                         if(ncol(saminfo) < 2)
                           stop("Only one column present in sample file: No extra information to add.")
                         existingVar = colnames(saminfo)[-which(colnames(saminfo) == "ID")] %in% colnames(private$samInfo)[-which(colnames(private$samInfo) == "ID")]
                         if(any(existingVar))
                           stop(paste("Variables already present in the GRM object:",
                                      paste(colnames(saminfo)[-which(colnames(saminfo) == "ID")][existingVar], collapse = "\n  "), sep="\n"))
                         if(nrow(saminfo) != nrow(private$samInfo))
                           stop(paste0("different number of individuals in sample file than in GRM object. Expected ", nrow(private$samInfo), " rows."))
                         if(any(!(saminfo$ID %in% private$samInfo$ID)))
                           stop(paste("The following sample IDs not found in the GRM object:",
                                       saminfo$ID[!(saminfo$ID %in% private$samInfo$ID)],"\n", sep="\n  "))
                         if(any(duplicated(saminfo$ID)))
                           stop(paste("Sample ID(s) duplicated in the sample information file:",
                                      paste(saminfo$ID[duplicated(saminfo$ID)], collapse = "\n  "), sep="\n  "))

                         ## Add sample information
                         private$samInfo = merge(private$samInfo, saminfo, by = "ID", sort=FALSE)
                         return(invisible(NULL))
                       },
                       ############# Delete group information
                       removeSampleInfo = function(name){
                         if(!is.vector(name) || !is.character(name))
                           stop("Argument 'name' needs to be a character vector")
                         else if(!(name %in% colnames(private$samInfo)))
                           stop("Variable not found. Check the name of the group.")
                         else if(any(name %in% c("ID","Ploidy")))
                           stop("Cannot drop the 'ID' or 'Ploidy' columns")
                         private$samInfo = subset(private$samInfo, select = setdiff(colnames(private$samInfo), name))
                       },
                       ########### Extract GRM to file
                       extractGRM = function(name, IDvar=NULL){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1)
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         if(!(name %in% names(private$GRM)))
                           stop("Cannot find GRM. Chech the 'name' argument")
                         
                         ## extract GRM and put the specified IDs on rows and columns
                         GRM <- private$GRM[[name]]$GRM
                         if(!is.null(IDvar)){
                           match.arg(IDvar, names(private$samInfo))
                           colnames(GRM) <- rownames(GRM) <- private$samInfo[[IDvar]]
                         } else colnames(GRM) <- rownames(GRM) <- private$indID
                         
                         ## write GRM to file
                         return(GRM)
                       },
                       ########### Write GRM to file
                       writeGRM = function(name, filename, IDvar=NULL){
                         if(!is.vector(name) || !is.character(name) || length(name) != 1)
                           stop("Argument 'name' needs to be a character vector of length 1.")
                         if(!(name %in% names(private$GRM)))
                           stop("Cannot find GRM. Chech the 'name' argument")

                         ## extract GRM and put the specified IDs on rows and columns
                         GRM <- private$GRM[[name]]$GRM
                         if(!is.null(IDvar)){
                           match.arg(IDvar, names(private$samInfo))
                           colnames(GRM) <- rownames(GRM) <- private$samInfo[[IDvar]]
                         } else colnames(GRM) <- rownames(GRM) <- private$indID

                         ## write GRM to file
                         data.table::fwrite(GRM, file = filename)
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
                       samInfo = NULL
                     )
)
