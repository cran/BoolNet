# Plot a wiring graph of the network <network> with the supplied
# graphical parameters.
# Requires igraph.
# Returns the igraph structure representing the wiring graph. 
plotNetworkWiring <- function(network,layout=layout.fruchterman.reingold,plotIt=TRUE,...)
{
	stopifnot(inherits(network,"BooleanNetwork"))
	
	if (!require(igraph))
		stop("Please install the igraph package before using this function!")
		
	edgeList <- c()
	
	# construct list of edges from interactions
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

  # build graph from edge list
	res <- graph.data.frame(edgeList-1,directed=TRUE,vertices=0:(length(network$genes)-1))
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
