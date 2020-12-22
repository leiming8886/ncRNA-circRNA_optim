#' @title Plot of the subnetwork
#'
#' @description The function plots a network from graphNEL or igraph format. It is used to visualize the modules. For further plotting options use the plot.igraph
#' function of the igraph package. The shapes of the nodes can be changed according to the scores argument, then negative scores appear squared
#' The color of the nodes can be changed according to the diff.expr argument. Negative(positive) values lead to green(red) nodes.
#'
#' @param network A graph in igraph or graphNEL format.
#' @param layout Layout algorithm, e.g. layout.fruchterman.reingold or layout.kamada.kawai.
#' @param labels Labels for the nodes of the network
#' @param diff.expr Named numerical vector of log2FC of the nodes in the network for coloring of the nodes.
#' @param scores Named numerical vector of scores of the nodes for the shape of the node in the network.
#' @param main Main title of the plot.
#' @param vertex.size Numerical value or verctor for the size of the vertices.

#' @examples
#' data(node_attr_Example)
#' data("network_Example")
#' gene2weight <- combinp(node_attr_Example)
#' gene2weight = gene2weight[1:100,]
#' res.list_Heinz <- runmodule(network = network_Example, gene2weight, method = "Heinz")
#' logFC_temp <- node_attr_Example$logFC
#' names(logFC_temp) <- rownames(node_attr_Example)
#' weigth_temp <- gene2weight$weight
#' names(weigth_temp) <- rownames(weigth_temp)
#' plotSub(res.list_Heinz$module, scores = weigth_temp, diff.expr = logFC_temp)

#'

#' @references Daniela Beisser, Gunnar W. Klau, Thomas Dandekar et al. (2010) BioNet: an R-Package for the functional analysis of biological networks
#'
#' @export

plotSub <- function (network, layout = layout.fruchterman.reingold, labels = NULL,
                            diff.expr = NULL, scores = NULL, main = NULL, vertex.size = NULL)
{
  if (is(network, "graphNEL")) {
    network <- igraph.from.graphNEL(network)
  }
  if (is.null(V(network)$name)) {
    V(network)$name <- as.character(V(network))
  }
  if (is.null(labels)) {
    if ("geneSymbol" %in% list.vertex.attributes(network)) {
      labels <- V(network)$geneSymbol
    }
    else {
      labels <- V(network)$name
    }
  }
  shapes <- rep("circle", length(V(network)))
  names(shapes) <- V(network)$name
  if (!is.null(scores) && !is.null(names(scores))) {
    shapes[intersect(names(which(scores < 0)), V(network)$name)] <- "csquare"
  }
  if (is.null(scores) && "score" %in% list.vertex.attributes(network)) {
    scores <- V(network)$score
    names(scores) <- V(network)$name
    shapes[names(which(scores < 0))] <- "csquare"
  }
  if (!is.null(diff.expr) && !is.null(names(diff.expr))) {
    coloring <- node.color(network, diff.expr)
  }
  else {
    coloring <- "SkyBlue2"
  }
  if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network)) {
    diff.exprs = V(network)$diff.expr
    names(diff.exprs) <- V(network)$name
    coloring <- node.color(network, diff.exprs)
  }
  max.labels <- max(nchar(labels))
  network.size = length(V(network))
  vertex.size2 <- 8
  cex = 0.6
  if (network.size < 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      vertex.size2 <- 15
      labels.dist <- 0
    }
  }
  if (network.size < 100 && network.size >= 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (network.size >= 100) {
    if (max.labels > 3) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (!is.null(vertex.size)) {
    vertex.size2 <- vertex.size
    labels.dist <- vertex.size/15
  }
  plot(network, layout = layout, vertex.size = vertex.size2,
       vertex.label = labels, vertex.label.cex = cex, vertex.label.dist = labels.dist,
       vertex.color = coloring, vertex.label.family = "sans",
       vertex.shape = shapes, main = main)
}
