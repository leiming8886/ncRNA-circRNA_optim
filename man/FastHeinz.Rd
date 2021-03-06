% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastHeinz.R
\name{FastHeinz}
\alias{FastHeinz}
\title{Calculate maximum scoring heuristically subnetwork}
\usage{
FastHeinz(network, scores)
}
\arguments{
\item{network}{A graph in igraph or graphNEL format.}

\item{scores}{A named vector, containing the scores for the nodes of the network.}
}
\value{
A subnetwork in the input network format.
}
\description{
The function uses an heuristic approach to calculate the maximum scoring subnetwork. Based on the given network and scores, the positive nodes are in the first step aggregated to meta-nodes between which minimum spanning trees are calculated.
}
\references{
Daniela Beisser, Gunnar W. Klau, Thomas Dandekar et al. (2010) BioNet: an R-Package for the functional analysis of biological networks
}
