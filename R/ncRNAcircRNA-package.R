#' ncRNAcircRNA: Search for modules in a node-weighted ncRNA network
#'
#' ncRNAcircRNA constructs a node-weighted ncRNA network, performs module searching, generates simulation data from random networks,
#' normalizes module scores using simulation data, removes un-qualified modules, and orders resultant modules according to their significance.
#'
#' @details
#' This package takes three types of data as input: a list of genes with
#' association p-values and logFC, a human ncRNA network. generate_graph constructs a
#' node-weighted ncRNA network. runmodule performs  module
#' search upon the node-weighted ncRNA network.
#'
#' @references
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @references Peilin Jia, Siyuan Zheng, Jirong Rong, Wei Zheng, Zhongming Zhao. (2011) Bioformatics. dmGWAS: dense module searching for genome-wide association studies in protein-protein interaction networks.
#' @docType package
#' @name ncRNAcircRNA-package
#' @aliases ncRNAcircRNA
#' @import igraph
#' @importFrom methods is
#' @importFrom stats pnorm sd qnorm na.omit
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot
NULL
