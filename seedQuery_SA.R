seedQuery_SA = function (G, seed, T = 1, T_min = 0.5, r = 0.98) 
{
  if (!require(igraph)) {
    stop("igraph must be pre-installed!\n")
  }
  net <- G
  d <- search_r
  if (!is.element("name", list.vertex.attributes(net))) {
    stop("Graph must have 'name' attribute")
  }
  if (!is.element("weight", list.vertex.attributes(net))) {
    stop("Graph must have 'weight' attribute")
  }
  subG <- induced.subgraph(net, seed)
  if (!is.connected(subG)) {
    stop("Input seeds are disjoint")
  }
  in.nodes <- V(subG)$name
  ############################################
  
  ###########################################	
  while( T > T_min){
    subx <- V(subG)$name
    subsum <- sum(V(subG)$weight)/sqrt(length(subx))
    tmp.neigh <- unlist(neighborhood(net, order = 1, 
                                     nodes = V(subG)$name))
    pot.nodes <- V(net)[tmp.neigh]$name
    pot.nodes <- setdiff(pot.nodes, in.nodes)
    if (length(pot.nodes) == 0) 
      break
    sub.weg <- V(net)[pot.nodes]$weight
    best.nodes <- pot.nodes[which(sub.weg == max(sub.weg))]
    subsum.u <- (sum(V(subG)$weight) + V(net)[best.nodes[1]]$weight)/sqrt(length(subx) + 
                                                                            1)
    dE = subsum.u - subsum
    if (dE >= 0) {
      tmp <- unlist(lapply(best.nodes, function(x) node2treePath(net, 
                                                                 V(subG)$name, x)))
      in.nodes <- c(tmp, V(subG)$name)
      subG <- induced.subgraph(net, in.nodes)
    }
    else {
      if ( exp( dE/T ) > runif(1,0,1) ){
        tmp <- unlist(lapply(best.nodes, function(x) node2treePath(net, 
                                                                   V(subG)$name, x)))
        in.nodes <- c(tmp, V(subG)$name)
        subG <- induced.subgraph(net, in.nodes) ;  #接受从Y(i)到Y(i+1)的移动
      } 
    }
    T = r * T
  }
  return(subG)
}
