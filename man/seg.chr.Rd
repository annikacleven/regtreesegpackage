% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seg.chr.R
\name{seg.chr}
\alias{seg.chr}
\title{Segment and Visualize a Chromosome}
\usage{
seg.chr(df, chromid)
}
\arguments{
\item{df}{A dataframe with columns of Start.Pos, log2r, and Chr columns. The Chr column should have format like "chr1", "chr21", "chrY".}

\item{chromid}{The chromosome id, for example "chr5"}
}
\value{
A list containing a dataframe with all of the segmentation data (segments),
and a dataframe with the predictions from the regression tree (regtreepred),
and a plot of the segmentation and chromosome data
}
\description{
Creates a plot of chromosome segmented, a dataframe including all the segments in the chromosome,
and a dataframe including the predictions from the regression tree using the optimal cp value for the chromosome.
}
\examples{
example <- seg.chr(datafr, "chr5")
example$regtreepred
example$segments
example$plot
}
\author{
Annika Cleven
}
