##########################################################################
# Genotyping Uncertainty with Sequencing data and RELATEdness (GUSrelate)
# Copyright 2019-2024 Timothy P. Bilton
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

#' @title GRM object
#' 
#' @description
#' Class for storing RA data and associated functions for analyzing unstructured populations (e.g.,
#' populations with no known structure).
#' 
#' @details  
#' A genomic relationship matrix (GRM) object is created from the \code{\link{makeGRM}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (or methods)
#' for analyzing the data. Information in a GRM object are specific to constructing a genomic relationship matrix
#'
#' @author Timothy P. Bilton
#' @seealso \code{\link{makeGRM}}
GRM <- R6::R6Class("GRM",
                     inherit = GUSbase::RA,
                     public = list(
                       #' @description
                       #' Method for initializing GRM object.
                       #' 
                       #' @param RAobj Existing RA object to use in initializing the GRM object.
                       #' @param ploid Vector of integer values indicating the ploid levels of each individual.
                       #' @param indsubset Vector of indices specifying which individuals from the RA object to retain.
                       #' @param saminfo Data frame of sample information.
                       #'
                       initialize = function(RAobj, ploid, indsubset, saminfo){
                         private$ref       <- RAobj$.__enclos_env__$private$ref[indsubset,]
                         private$alt       <- RAobj$.__enclos_env__$private$alt[indsubset,]
                         private$chrom     <- RAobj$.__enclos_env__$private$chrom
                         private$pos       <- RAobj$.__enclos_env__$private$pos
                         private$SNP_Names <- RAobj$.__enclos_env__$private$SNP_Names
                         private$indID     <- RAobj$.__enclos_env__$private$indID[indsubset]
                         private$nSnps     <- RAobj$.__enclos_env__$private$nSnps
                         private$nInd      <- length(indsubset)
                         private$gform     <- RAobj$.__enclos_env__$private$gform
                         private$AFrq      <- RAobj$.__enclos_env__$private$AFrq
                         private$infilename<- RAobj$.__enclos_env__$private$infilename
                         private$ploid     <- as.integer(ploid)
                         private$samInfo   <- saminfo
                       },
                       ########### Print function:
                       #' @description
                       #' Print output from a GRM object
                       #' 
                       #' @details 
                       #' Prints the summary information of the dataset and if any GRMs have been constructed, it 
                       #' also prints out summary information of the GRM runs.
                       #' 
                       print = function(){
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
		                   #' @description
		                   #' Construct a genomic relationship matrix (GRM)
		                   #' 
		                   #' @details 
		                   #' The filtering criteria currently implemented are:
		                   #' \itemize{
		                   #' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is \code{NULL})}
		                   #' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
		                   #'  is greater than the threshold value (default is \code{NULL}).}
		                   #' \item{P-value from a HWE test(\code{PVALUE}):}{ SNPs are discarded if the p-value from the Hardy-Weinberg equilibrium test is less than
		                   #' the threshold. (default is \code{NULL})}
		                   #' }
		                   #' If a filtering criteria is set to \code{NULL}, then no filtering in regard to
		                   #' that threshold is applied.
		                   #' 
		                   #' Not that for the \code{PVALUE} filter, the HWE test must first be run using the \code{$HWEtest} method.
		                   #' 
		                   #' @param name Specific name given to the run constructing the GRM.
		                   #' @param method Character value indicating which method to used to construct the GRM.
		                   #' @param ep Sequencing error rate. Can be a single value, a vector equal to the number of SNPs
		                   #' or a matrix the same dimension as the data
		                   #' @param snpsubset Vector of indices indicated which SNPs to retain to construct the GRM.
		                   #' @param filter Name list specifying the filtering criteria to be applied.
		                   #' @param ... Not used currently.
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
		                   #' @description
		                   #' Perform a Hardy–Weinberg Equilibrium (HWE) test for autopolyploids.
		                   #' 
		                   #' @details 
		                   #' The function calls the \code{\link[GUSbase]{p_est_em}} and \code{\link[GUSbase]{g_est_em}} functions in the 
		                   #' \code{\link[GUSbase]{GUSbase}} package to perform the HWE test. 
		                   #' 
		                   #' @param snpsubset Vector of SNP indices indicating which SNPs to do the HWE test on.
		                   #' Excluded SNPs are given a p-value of 0.
		                   #' @param indsubset Vector of indices of the samples indicating which samples to use in the HWE test.
		                   #' @param nThreads Integer value specifying the number of threads to use in the parallelization.
		                   #' @param para Starting values passed to the \code{\link[GUSbase]{p_est_em}} and \code{\link[GUSbase]{g_est_em}} functions.
		                   #' @param EMpara Convergence criteria passed to the \code{\link[GUSbase]{p_est_em}} and \code{\link[GUSbase]{g_est_em}} functions.
		                   #' 
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
                         ## if only performed on a subset, then resize the vector to the number of SNPs
                         if(is.null(snpsubset)) pvalue = 1-pchisq(-2*(pest$loglik - gest$loglik), df=ploidy-1)
                         else{
                           pvalue = rep(0, private$nSnps)
                           pvalue[snpsubset] = 1-pchisq(-2*(pest$loglik - gest$loglik), df=ploidy-1)
                         }
                         private$pvalue <- pvalue
                       },
                       ############# Delete a computed GRM:
		                   #' @description
		                   #' Delete a GRM run within the GRM object.
		                   #' 
		                   #' @param name Name of the GRM run to delete.
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
                       #plotGRM = function(type=c("histrogram","scatter"),values=c("both","self-relatedness","relatednes"), xaxis=NULL,
                       #                    interactive=FALSE){
                       #   stop("yet to be implemeted")
                       #},
                       ############ Compare different GRMs:
                       #CompareGRM = function(){
                       #   stop("yet to be implemeted")
                       #},
                       ############# Plots for the GRM
		                   #' @description
		                   #' Produce a PCA plot of a constructed GRM
		                   #' 
		                   #' @details 
		                   #' The sample information used by this function is the information added by the `$addSampleInfo` method.
		                   #' 
		                   #' @param name Character value giving the name of the GRM run to extract the constructed GRM from.
		                   #' @param npc Integer number specifying the number of PCs to plot
		                   #' @param colour Character string specifying which column of the sample information to use to colour the points
		                   #' @param shape Character string specifying which column of the sample information to use to group the points based on point shape
		                   #' @param group.hover Character vector specifying which column(s) of the sample information to include in the hover information for the
		                   #' interactive plots. Only used if \code{interactive=TRUE}.
		                   #' @param interactive Logical value: If \code{TRUE}, then an interactive plot using the \code{plotly} package is produced.
		                   #' Otherwise, a standard \code{ggplot} is constructed.
		                   #' 
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
                             p = p + ggplot2::guides(color=ggplot2::guide_legend(title=colour))
                           if(!is.null(shape))
                             p = p + ggplot2::guides(shape=ggplot2::guide_legend(title=shape))
                           p
                         }
                       },
                       #### add group information
		                   #' @description
		                   #' Add sample information to the GRM object
		                   #' 
		                   #' @details 
		                   #' The sample information must contain a column 'ID' that has the ID of the samples and must match
		                   #' the IDs of the samples already in the GRM object. Any extra columns will be added to the existing 
		                   #' sample information data.
		                   #' 
		                   #' Note that a column named 'Ploidy' can not be present in the file. If there are any columns already present
		                   #' then the function throw and error. 
		                   #' 
		                   #' @param samfile Path to the file that contains the sample information to be loaded to GRM object.
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
		                   #' @description
		                   #' Remove sample information from a GRM object
		                   #' 
		                   #' @details 
		                   #' Note that columns 'ID' and 'Ploidy' cannot be deleted.
		                   #' 
		                   #' @param name Character vector of column names to remove from the sample information data.
		                   #' 
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
		                   #' @description
		                   #' Extract a GRM from a run.
		                   #' 
		                   #' @param name Character value specifying which run to extract the GRM from.
		                   #' @param IDvar Character value specifying which column of the sample information to use for 
		                   #' IDs in the GRM matrix. If \code{IDvar=NULL}, then the 'ID' column is used.
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
		                   #' @description
		                   #' Write GRM to file.
		                   #' 
		                   #' @details
		                   #' The file format is currently a symmetric matrix with the IDs on the rows and columns.
		                   #' 
		                   #' @param name Character value specifying which run to write the GRM from.
		                   #' @param filename Character value giving the name of the output file.
		                   #' @param IDvar Character value specifying which column of the sample information to use for 
		                   #' IDs in the GRM matrix. If \code{IDvar=NULL}, then the 'ID' column is used.
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
