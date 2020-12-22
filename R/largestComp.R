#' @title Extract largest component of network
#'
#' @description The function extracts the largest component of a network.
#' @param network A graph in igraph format
#'
#' @return A new graph object that represents the largest component of the given network.
#'
#' @examples
#' data("network_Example")
#' largestComp(network_Example)
#'
#' @export


largestComp <- function (network)
{
    #if (!require(igraph)) {
    #    stop("igraph must be pre-installed!\n")
    #}
    if (is(network, "graphNEL")) {
        cc <- RBGL::connectedComp(network)
        idx <- which.max(Biobase::listLen(cc))
        return(graph::subGraph(cc[[idx]], network))
    }
    else if (is(network, "igraph")) {
        clust <- clusters(network)
        cid <- which.max(clust$csize)
        lg <- induced.subgraph(network, V(network)[clust$membership ==
            cid])
        return(lg)
    }
}
