plotSequence <- function(network, startState, 
                         includeAttractorStates = c("all","first","none"), 
                         sequence,
                         title = "", mode=c("table","graph"),
                         plotFixed = TRUE, grouping = list(),
                         onColor="green",offColor="red",
                         layout, drawLabels=TRUE, drawLegend=TRUE, ...)
{
  if (!missing(network))
  {
    stopifnot(inherits(network,"BooleanNetwork"))
    if (missing(startState) || !missing(sequence))
      stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
    
    sequence <- getPathToAttractor(network = network,
                                   state = startState, 
                                   includeAttractorStates = includeAttractorStates)

    numGenes <- length(network$genes)
    if (match.arg(mode,c("table","graph")) == "table")
    {                                       
      whichFixed <- which(network$fixed != -1)
      if (plotFixed | (length(whichFixed) == 0))
        plotIndices <- 1:numGenes
      else
        plotIndices <- (1:numGenes)[-whichFixed]
        
      sequence <- sequence[,plotIndices]
    }
  }
  else
  {
    if (missing(sequence) || !missing(startState))
        stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")        
  }

  switch(match.arg(mode,c("table","graph")),
  table =
  {
    # initialize with empty plot
      plot(c(),c(),xlim=c(0,nrow(sequence)),ylim=c(-2,ncol(sequence)),xlab="",ylab="",
           axes=FALSE, ...)
  
      # build accumulated matrix
      totalMatrix <- t(sequence)
      
      if(length(grouping)>0)
           # reorder genes according to the supplied groups
          totalMatrix = totalMatrix[unlist(grouping$index),]
      
      colnames(totalMatrix) <- 1:ncol(totalMatrix)
    
      axis(3,c(0,(1:ncol(totalMatrix))-0.5),c("t=",colnames(totalMatrix)), 
           lty="blank", yaxt='s', xaxt='s', xaxs="i")
      axis(2,(1:nrow(totalMatrix))-0.5,rownames(totalMatrix),
           yaxt='s', xaxt='s', xaxs="i", las=2)

        # plot active and inactive states
      for(i in 1:ncol(totalMatrix))
        for(j in 1:nrow(totalMatrix))
        {
          if(totalMatrix[j,i] == 1)
            rect(i-1,j-1,i,j,col=onColor,border="gold")
          else
            rect(i-1,j-1,i,j,col=offColor,border="gold")
        }
        
        if(length(grouping)>0)
        # draw separators between groups, and print group names
        {
          sepPos = cumsum(sapply(grouping$index,length))
          abline(h=sepPos[-length(sepPos)],col="black",lwd=3)
          text(ncol(totalMatrix)/2,sepPos-0.5,grouping$class,cex=0.9)
        }
      
      if (drawLegend)  
        legend(x="bottomright",pch=c(15,15),
               col=c(onColor,offColor),
               legend = c("active","inactive"),
               cex=0.7,
               horiz=T)   
    
    return(totalMatrix) 
 },
  
  graph = 
  {
    if (!require(igraph))
      stop("Please install the igraph package before using the \"graph\" mode of this function!")
    
   args <- list(...)
      
    if (is.null(args$vertex.size))
      args$vertex.size <- 2
    
    if (is.null(args$edge.arrow.mode))
      args$edge.arrow.mode <- 0
      
    if (is.null(args$rescale))
      args$rescale <- !missing(layout)  
    
    if (is.null(args$vertex.label.cex))
      args$vertex.label.cex <- 0.75
    
    if (is.null(args$vertex.label.dist))
      args$vertex.label.dist <- 0.25
      
      if (is.null(args$vertex.label.degree))
        args$vertex.label.degree <- -pi/2        
      
     if (is.null(args$vertex.color))
      args$vertex.color <- "grey"
    
    if (is.null(args$edge.arrow.size))
      args$edge.arrow.size <- 0.5
      
    if (missing(layout))
      layout <- matrix(c(seq(-1,1,length.out=nrow(sequence)), rep(0,nrow(sequence))), ncol=2)        
      
    states <- apply(sequence,1,paste,collapse="")
    nodes <- data.frame(unique(states),stringsAsFactors=FALSE)
  
    if (length(states) > 1)
    {
      initialStates <- states[1:(length(states) - 1)]
      nextStates <- states[2:length(states)]
      edgeMatrix <- data.frame(initialStates,nextStates)
    }
    else
    {
      edgeMatrix <- data.frame(matrix(nrow=0,ncol=2))
    }
    
    graph <- graph.data.frame(edgeMatrix,vertices=nodes,directed=TRUE)
    
    if (drawLabels)
      labels <- nodes[,1]
    else
      labels <- NA
    
     plot(graph,layout=layout,vertex.label=labels,vertex.label.cex=args$vertex.label.cex,
       vertex.size=args$vertex.size, vertex.color=args$vertex.color,
       vertex.label.dist=args$vertex.label.dist,
       vertex.label.degree=args$vertex.label.degree,
       edge.arrow.size=args$edge.arrow.size,
       rescale=args$rescale,
       main=title,...)
    return(graph)
  })
}
