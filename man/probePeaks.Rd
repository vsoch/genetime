% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/probePeaks.R
\name{probePeaks}
\alias{probePeaks}
\title{probePeaks}
\usage{
probePeaks(genelist, thresh = 0.95, out.pdf = FALSE,
  out.name = "genetime.pdf")
}
\arguments{
\item{genelist}{A list of gene symbol identifiers}

\item{thresh}{The quantile of expression values to take per gene.  Default is 0.95}

\item{out.pdf}{output pdf of temporal enrichment patterns, defaults to FALSE.  Default output name (set with out.name) is "genetime.pdf"}

\item{out.name}{name of output pdf of temporal enrichment patterns, default is "genetime.pdf"}
}
\description{
This function will extract temporal features for a gene set of interest to predict times of temporal enrichment.
}
\examples{
mygenes = "ADNP" "AGAP2" "ANK2" "APH1A"
thresh = 0.99
probePeaks(genelist=mygenes,threshold=thresh,out.pdf=TRUE,out.name="/home/myfile.pdf")
}
\keyword{enrichment}
\keyword{gene}
\keyword{temporal}

