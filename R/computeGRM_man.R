##########################################################################
# Genotyping Uncertainty with Sequencing data for RELATEdness (GUSrelate)
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
##' GRM method: Construct a genomic relationship matrix (GRM)
##'
##' Method for constructing a genomic relationship matrix (GRM) for a diploid or autopolyploid population.
##'
##' The filtering criteria currently implemented are:
##' \itemize{
##' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is \code{NULL})}
##' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
##'  is greater than the threshold value (default is \code{NULL}).}
##' \item{h (\code{PVALUE}): }{}
##' }
##' If a filtering criteria is set to \code{NULL}, then no filtering in regard to
##' that threshold is applied.
##'
##' @usage
##' GRMobj$computeGRM(name, method="VanRaden", ep=0, snpsubset=NULL,
##'                   filter=list(MAF=NULL, MISS=NULL, PVALUE=NULL), ...)
##'
##' @param name A character string giving the name of the GRM analysis.
##' @param method A character string specifying whether the VanRaden (\code{'VanRaden'}) based estimator or
##' the Weir-Goudet (\code{'WG'}) estimator is used to construct the GRM.
##' @param ep Sequencing error value. Can be a single number (sequencing error rate same for all genotypes),
##' a vector equal to the number of SNPs (SNP specific sequencing error rate) or a matrix the same diminsion as the data (SNP and
##' individual specific sequencing error rate).
##' @param snpsubset Vector of indices of SNPs that should be retained in the analysis. Useful for subsetting the SNPs before
##' any filtering is applied.
##' @param filter Named list of thresholds for various filtering criteria.
##' See below for details.
##' @param ... Not used currently.
##'
##' @name $computeGRM
##' @seealso \code{\link{GRM}}
##' @author Timothy P. Bilton
#NULL
