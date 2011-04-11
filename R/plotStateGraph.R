# Plots a graph that visualizes the state transitions and attractor basins. <attractorInfo> is an object
# of class AttractorInfo. This requires the igraph package.
# If <highlightAttractors> is set, attractor edges are drawn bold.
# If <colorBasins> is true, each basin is drawn in a different color. 
# Colors can be provided in <colorSet>.
# <layout> specifies the graph layouting function.
# If <piecewise> is true, subgraphs are layouted separately.
# <basin.lty> and <attractor.lty> specify the line types used to draw states in the basins
# and in the attractors (if <highlightAttractor> is set).
# If <plotIt> is not set, only the igraph object is returned, but no graph is plotted.
# ... provides further graphical parameters for the plot.
# Returns an object of class igraph
plotStateGraph <- function(attractorInfo,highlightAttractors=TRUE,colorBasins=TRUE,colorSet,
             drawLegend=TRUE,drawLabels=FALSE,layout=layout.fruchterman.reingold,
             piecewise=FALSE,basin.lty=2,attractor.lty=1,plotIt=TRUE,...)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))

  if (!require(igraph))
    stop("Please install the igraph package before using this function!")

  if (is.null(attractorInfo$stateInfo$table))
    stop(paste("This AttractorInfo structure does not contain transition table information.",
           "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))
  
  graphStruct <- getStateGraphStructure(attractorInfo)
  
  res <-   graph.edgelist(el=graphStruct$edges - 1)
  
  res <- set.vertex.attribute(res,"name",value=graphStruct$vertices)
  
  # determine nodes and edges that belong to attractors
  attractorIndices <- which(attractorInfo$stateInfo$stepsToAttractor == 0)
  
  attractorEdgeIndices <- which(apply(graphStruct$edges,1,
             function(edge)((edge[1] %in% attractorIndices) & (edge[2] %in% attractorIndices)))) - 1

  # set default edge width and line type
  res <- set.edge.attribute(res,"width",value=0.9)
  res <- set.edge.attribute(res,"lty",value=basin.lty)
  
  if (highlightAttractors)
  {
    # set different edge width and line type for attractor edges
    res <- set.edge.attribute(res,"width",index=attractorEdgeIndices,value=2)
    res <- set.edge.attribute(res,"lty",index=attractorEdgeIndices,value=attractor.lty)
  }

  if (missing(colorSet))
  {
    # define default colors
    colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
        "purple","orange","seagreen","tomato","darkgray","chocolate",
        "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
        "burlywood2")
  }
  args <- list(...)
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.size))
    args$vertex.size <- 2
    
  if (is.null(args$edge.arrow.mode))
    args$edge.arrow.mode <- 0
    
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.5
    
  if (is.null(args$vertex.label.dist))
    args$vertex.label.dist <- 1
  
  if (colorBasins)
  {  
    for (attractor in 1:length(attractorInfo$attractors))
    {
      # determine nodes and edges belonging to the basin of <attractor>
      basinIndices <- which(attractorInfo$stateInfo$attractorAssignment == attractor)
      
      # change vertex color
      res <- set.vertex.attribute(res,"color",basinIndices - 1,
                value=colorSet[(attractor-1) %% length(colorSet) + 1])
      if (drawLabels)
        res <- set.vertex.attribute(res,"label.color",basinIndices - 1,
            value=colorSet[(attractor-1) %% length(colorSet) + 1])
      basinEdgeIndices <- which(apply(graphStruct$edges,1,
                 function(edge)((edge[1] %in% basinIndices) 
                         & (edge[2] %in% basinIndices)))) - 1
      
      # change edge color
      res <- set.edge.attribute(res,"color",index=basinEdgeIndices,
              value=colorSet[(attractor-1) %% length(colorSet) + 1])
    }
  }

  if(plotIt)
  {
    if (drawLabels)
      labels <- graphStruct$vertices
    else
      labels <- NA
    if (piecewise)
      layout <- piecewise.layout(res, layout)
    plot(res,vertex.size=args$vertex.size,layout=layout,
         edge.arrow.mode=args$edge.arrow.mode,
         vertex.label=labels,vertex.label.cex=args$vertex.label.cex,
         vertex.label.dist=args$vertex.label.dist,
         ...)
    if (colorBasins & drawLegend)
      legend(x="bottomleft",pch=15,ncol=1,
             col=colorSet[0:(length(attractorInfo$attractors) - 1) %% length(colorSet) + 1],
             legend = paste("Attractor",1:length(attractorInfo$attractors)),
             cex=0.5)
  }
  return(invisible(res))
}
