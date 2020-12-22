#' @title Generate corrected p-value based on p-value and logFC in the expression matrix
#'
#' @description based on the formula: corrected p-value = 2*( 1-pnorm( (-log10( p-value )) * abs( log2FC )) ), corrected p-value was generated
#'
#' @param node_attr logFC and p value in the two columns
#'
#' @return the matrix of two columns: gene, weight(corrected p-value)
#'
#' @examples
#' data("node_attr_Example")
#' result <- combinp(node_attr_Example)
#'
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @export
combinp <- function(node_attr = NULL){
	#dir <- "D:/biaoshu/cirRNA_add/算法subnetwork/BioNet/pval.txt"
	#node_attr <- read.table(dir,row.names=1, sep= "\t")
	#gene2weight <- seq(1, length(node_attr[,1]))
	node_p2 <-apply(node_attr,1,function(x) 2*( 1-pnorm( (-log10(unlist(x[2]))) * abs(unlist(x[1]))) ))
	gene2weight <- data.frame(gene=names(node_p2),weigth=node_p2)
	colnames(gene2weight) <- c("gene","weight")
	return(gene2weight)
}
