% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GUSrelate.R
\name{GRM}
\alias{GRM}
\title{GRM object}
\usage{
GRMobj <- makeGRM()
}
\description{
Class for storing RA data and associated functions for analyzing unstructured populations (e.g.,
populations with no known structure).
}
\details{
A GRM object is created from the \code{\link{makeGRM}} function and contains RA data,
various statistics of the dataset that have been computed, and functions (or methods)
for analyzing the data. Information in an GRM are specific to constructing a GRM
}
\seealso{
\code{\link{makeGRM}}
}
\author{
Timothy P. Bilton
}
\section{Super class}{
\code{\link[GUSbase:RA]{GUSbase::RA}} -> \code{GRM}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{GRM$new()}}
\item \href{#method-computeGRM}{\code{GRM$computeGRM()}}
\item \href{#method-p_est}{\code{GRM$p_est()}}
\item \href{#method-HWE_est}{\code{GRM$HWE_est()}}
\item \href{#method-PCA}{\code{GRM$PCA()}}
\item \href{#method-addSampleInfo}{\code{GRM$addSampleInfo()}}
\item \href{#method-deleteSampleInfo}{\code{GRM$deleteSampleInfo()}}
\item \href{#method-clone}{\code{GRM$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="cometPlot">}\href{../../GUSbase/html/RA.html#method-cometPlot}{\code{GUSbase::RA$cometPlot()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="extractVar">}\href{../../GUSbase/html/RA.html#method-extractVar}{\code{GUSbase::RA$extractVar()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="print">}\href{../../GUSbase/html/RA.html#method-print}{\code{GUSbase::RA$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="RDDPlot">}\href{../../GUSbase/html/RA.html#method-RDDPlot}{\code{GUSbase::RA$RDDPlot()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="GUSbase" data-topic="RA" data-id="writeVCF">}\href{../../GUSbase/html/RA.html#method-writeVCF}{\code{GUSbase::RA$writeVCF()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$new(GRMobj, ploid, indsubset)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-computeGRM"></a>}}
\subsection{Method \code{computeGRM()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$computeGRM(
  name,
  method = "VanRaden",
  ep = 0,
  snpsubset = NULL,
  filter = list(MAF = NULL, MISS = NULL, PVALUE = NULL),
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-p_est"></a>}}
\subsection{Method \code{p_est()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$p_est(
  snpsubset = NULL,
  indsubset = NULL,
  nThreads = 1,
  para = NULL,
  EMpara = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-HWE_est"></a>}}
\subsection{Method \code{HWE_est()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$HWE_est(
  snpsubset = NULL,
  indsubset = NULL,
  nThreads = 1,
  para = NULL,
  EMpara = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PCA"></a>}}
\subsection{Method \code{PCA()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$PCA(
  name,
  npc = 3,
  group1 = NULL,
  group2 = NULL,
  group.hover = NULL,
  interactive = FALSE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addSampleInfo"></a>}}
\subsection{Method \code{addSampleInfo()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$addSampleInfo(name, info)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deleteSampleInfo"></a>}}
\subsection{Method \code{deleteSampleInfo()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$deleteSampleInfo(name)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GRM$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
