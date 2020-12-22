#' @title Run module search function
#'
#' @description runmodule constructs a node-weighted ncRNA network, performs module searching, generates simulation data from random networks,
#' @description normalizes module scores using simulation data, removes un-qualified modules, and orders resultant modules according to their significance.
#'
#' @param network A data frame containing a symbolic edge list of the ncRNA
#'   network
#' @param gene2weight weight A data frame containing two columns: the first is unique
#'   gene identifier (should be coordinate with the node symbol used in ncRNA network);
#'   the second is gene-based corrected p-value derived from differentially gene analysis
#' @param method	a character string indicating which the search method is to be computed
#'  . One of "Heinz" (default), "GS": can be abbreviated
#' @param d An integer used to define the order of neighbour genes to be
#'   searched if the method is "GS". This parameter is always set up as 2
#' @param r A float indicating the cut-off for increment during module expanding
#'   process if the method is "GS". Greater r will generate smaller module. Default is 0.1.
#'
#' @return \code{runmodule} returns a list containing relevant data and results,
#'   including:
#'
#'  \tabular{ll}{
#'    \code{GWPI} \tab the edge-weighted network used for searching \cr
#'    \code{graph.g.weight} \tab the edge-weighted of each gene in the edge-weighted network \cr
#'    \code{module} \tab list of genes comprising each module, named for the seed gene \cr
#'    \code{module.score.matrix} \tab contains Zm, Zn and Zcount \cr
#'    \code{Zi} \tab contains Zn \cr
#'  }
#'
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @references Peilin Jia, Siyuan Zheng, Jirong Rong, Wei Zheng, Zhongming Zhao. (2011) Bioformatics. dmGWAS: dense module searching for genome-wide association studies in protein-protein interaction networks.
#' @examples
#' \dontrun{
#'  data("node_attr_Example")
#'  data("network_Example")
#'  geneweight <- combinp(node_attr_Example)
#'  gene2weight = gene2weight[1:100,]
#'  res.list <- runmodule(network = network_Example, gene2weight, method = "Heinz")
#' }
#'
#' @export



runmodule <- function (network, gene2weight, d = 2, r = 0.1, method = c("Heinz","GS")) #c("GS","Heinz")
{
    #library(igraph)
    #library(stats)
    if (min(gene2weight[, 2]) <= 0 | max(gene2weight[, 2]) >=
        1) {
        stop("P values out of range 0 < p < 1")
    }
   if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
    method <- match.arg(method)
    cat("genes used: ", length(gene2weight[, 1]), "\n", sep = "")
    rawG <- graph.data.frame(network, directed = F)
    g.weight <- sapply(as.numeric(gene2weight[, 2]), function(x) qnorm(1 -
        x))
    intG <- integGM(rawG, as.character(gene2weight[, 1]), g.weight)
    GWPI <- simplify(intG)
    cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
	if( method == "Heinz"){
		scores <- V(GWPI)$weight
		names(scores)<- V(GWPI)$name
		module <- FastHeinz(GWPI, scores)
		genes.idx <- seq(1, length(V(GWPI)$name))
		graph.g.weight = data.frame(GWPIene = V(GWPI)$name, gain.weight = V(GWPI)$weight)
		genes <- V(module)$name
		idx <- match(genes, graph.g.weight[, 1])
		idx <- idx[!is.na(idx)]
		ZM <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
		l.zperm <- c()
		for (j in 1:1e+05) {
            idx.pseudo = sample(genes.idx, size = length(V(module)))
            l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
                2])/sqrt(length(idx.pseudo)))
        }
		k.mean <- mean(l.zperm)
		k.sd <- sd(l.zperm)
		ZN = (ZM - k.mean)/k.sd
		res.list <- list()
		res.list[["GWPI"]] = GWPI
		res.list[["module"]] = module
		res.list[["zi.ordered"]] = ZN
		cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
		return(res.list)

	}
  if( method == "GS"){
	sublist = list()
	for (node in V(GWPI)$name) {
		ng <- seedQueryJ_GS(GWPI, node,search_r = 2, r = 0.1)
		if (vcount(ng) >= 5)
			sublist[[node]] <- ng
		}
	dm.result <- sublist
    cat("extracting modules...\n", sep = "")
    genesets <- list()
    for (k in 1:length(dm.result)) {
        node = names(dm.result[k])
        g = dm.result[[k]]
        genesets[[node]] <- V(g)$name
    }
    seed.genes <- names(genesets)
    cat("removing identical modules...\n", sep = "")
    identical.idx <- list()
    for (k in 1:length(genesets)) {
        tmp.idx <- c(k)
        for (kt in 1:length(genesets)) {
            if (kt == k)
                (next)()
            genesk <- genesets[[k]]
            genest <- genesets[[kt]]
            if (length(genesk) != length(genest))
                (next)()
            overlap = intersect(genesk, genest)
            if (length(overlap) == length(genest)) {
                tmp.idx <- c(tmp.idx, kt)
            }
        }
        if (length(tmp.idx) > 1) {
            tmp.idx <- sort(tmp.idx)
            identical.idx[[seed.genes[k]]] <- tmp.idx
        }
    }
    toremove.idx <- c()
    for (k in 1:length(identical.idx)) {
        tmp.idx <- identical.idx[[k]]
        toremove.idx <- c(toremove.idx, tmp.idx[-1])
    }
    toremove.idx <- unique(toremove.idx)
    genesets.clear <- genesets[-toremove.idx]
    cat("permutation on random network...\n", sep = "")
    genesets.length <- c()
    for (k in 1:length(genesets.clear)) {
        genes <- genesets.clear[[k]]
        genesets.length <- c(genesets.length, length(genes))
    }

    genesets.length <- unique(genesets.length)
    genes.idx <- seq(1, length(V(GWPI)$name))
    graph.g.weight = data.frame(GWPIene = V(GWPI)$name, gain.weight = V(GWPI)$weight)
    genesets.length.null.dis <- list()
    length.max = max(genesets.length)
	length.min = min(genesets.length)
    for (k in length.min:length.max) {
        l.zperm <- c()
        for (j in 1:1e+05) {
            idx.pseudo = sample(genes.idx, size = k)
            l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
                2])/sqrt(length(idx.pseudo)))
        }
        genesets.length.null.dis[[as.character(k)]] = l.zperm
        cat(k, ".", sep = "")
    }
    genesets.length.null.stat <- list()
    for (k in length.min:length.max) {
        l.zperm <- genesets.length.null.dis[[as.character(k)]]
        k.mean <- mean(l.zperm)
        k.sd <- sd(l.zperm)
        genesets.length.null.stat[[as.character(k)]] = c(k.mean,
            k.sd)
    }
    zim <- data.frame(gene = names(genesets.clear), Zm = -9,
        Zn = -9, zcount = -9)
    for (k in 1:length(genesets.clear)) {
        genes <- genesets.clear[[k]]
        idx <- match(genes, graph.g.weight[, 1])
        idx <- idx[!is.na(idx)]
        zim[k, 2] <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
        tmp <- genesets.length.null.stat[[as.character(length(idx))]]
        zim[k, 3] = (zim[k, 2] - tmp[1])/tmp[2]
        zim[k, 4] = sum(genesets.length.null.dis[[as.character(length(idx))]] >=
            zim[k, 2])
    }
    zom = zim[order(zim[, 3], decreasing = T), ]
    res.list <- list()
    res.list[["GWNE"]] = GWPI
    res.list[["graph.g.weight"]] = graph.g.weight
    res.list[["module"]] = genesets.clear
    #res.list[["genesets.length.null.dis"]] = genesets.length.null.dis
    #res.list[["genesets.length.null.stat"]] = genesets.length.null.stat
    #res.list[["zi.matrix"]] = zim
    res.list[["module.score.matrix"]] = zom
    save(res.list, file = "RESULT.list.RData")
    cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
    return(res.list)
	}
}
