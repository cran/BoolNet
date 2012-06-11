# Plot a wiring graph of the network <network> with the supplied
# graphical parameters.
# Requires igraph.
# Returns the igraph structure representing the wiring graph. 
plotNetworkWiring <- function(network,layout=layout.fruchterman.reingold,plotIt=TRUE,...)
{
  stopifnot(inherits(network,"ProbabilisticBooleanNetwork") | inherits(network,"BooleanNetworkCollection")
            | inherits(network,"BooleanNetwork"))
  
  if (!require(igraph))
    stop("Please install the igraph package before using this function!")
  
  if (installed.packages()["igraph","Version"] < package_version("0.6"))
    bias <- 1
  else
    bias <- 0 
    
  edgeList <- c()
  
  # construct list of edges from interactions

  if (inherits(network,"BooleanNetwork"))
  # deterministic network
  {
    for (i in 1:length(network$genes))
    {
      if (network$interactions[[i]]$input[1] != 0)
      # no edges for constant genes
      {
        edgeList <- rbind(edgeList,
                  cbind(network$interactions[[i]]$input,
                  rep(i,length(network$interactions[[i]]$input))))
      }
    }
  }
  else
  # probabilistic network
  {
    for (i in 1:length(network$genes))
    {
      for (j in 1:length(network$interactions[[i]]))
      {
        if (network$interactions[[i]][[j]]$input[1] != 0)
        # no edges for constant genes
        {
          edgeList <- rbind(edgeList,
                            cbind(network$interactions[[i]][[j]]$input,
                            rep(i,length(network$interactions[[i]][[j]]$input))))
        }
      }
    }
  }

  # build graph from edge list
  res <- graph.data.frame(edgeList-bias,directed=TRUE,vertices=as.data.frame((1:length(network$genes)) - bias))
  res <- set.vertex.attribute(res,"name",value=network$genes)
  
  args <- list(...)
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.color))
    args$vertex.color <- "grey"
    
  if (is.null(args$edge.arrow.size))
    args$edge.arrow.size <- 0.5
    
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.7

  if (is.null(args$vertex.size))
    args$vertex.size <- 18    

  if (plotIt)
  {
    plot(res,vertex.label=network$genes,vertex.label.cex=args$vertex.label.cex,
         vertex.size=args$vertex.size,vertex.color=args$vertex.color,
         edge.arrow.size=args$edge.arrow.size,
         layout=layout,...)
  }
  return(invisible(res))
}
