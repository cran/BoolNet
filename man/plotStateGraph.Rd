\name{plotStateGraph}
\Rdversion{1.1}
\alias{plotStateGraph}
\title{
Visualize state transitions and attractor basins
}
\description{
Plots a graph containing all states visited in \code{stateGraph}, and optionally highlights attractors and basins of attraction. This requires the \pkg{igraph} package.
}
\usage{
plotStateGraph(stateGraph, highlightAttractors = TRUE, 
               colorBasins = TRUE, colorSet, 
               drawLegend = TRUE, drawLabels = FALSE, 
               layout = layout.fruchterman.reingold,
               piecewise = FALSE,  
               basin.lty = 2, attractor.lty = 1, 
               plotIt = TRUE, ...)
}

\arguments{
  \item{stateGraph}{An object of class \code{AttractorInfo} or \code{SymbolicSimulation}, as returned by \code{\link{getAttractors}} and \code{\link{simulateSymbolicModel}} respectively. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE. Similarly, \code{simulateSymbolicModel} must be called with \code{returnGraph=TRUE}. Alternatively, \code{stateGraph} can be an object of class \code{TransitionTable}, which can be extracted using the functions \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, or \code{\link{getStateSummary}}
} 
  \item{highlightAttractors}{
	If this parameter is true, edges in attractors are drawn bold and with a different line type (which can be specified in \code{attractor.lty}). Defaults to TRUE.
}
  \item{colorBasins}{
	If set to true, each basin of attraction is drawn in a different color. Colors can be specified in \code{colorSet}. Defaults to TRUE.
}
  \item{colorSet}{
	An optional vector specifying the colors to be used for the different attractor basins. If not supplied, a default color set is used.
}
  \item{drawLegend}{
	If set to true and \code{colorBasins} is true, a legend for the colors of the basins of attraction is drawn. Defaults to TRUE.
}
  \item{drawLabels}{
	If set to true, the binary encodings of the states are drawn beside the vertices of the graph. As this can be confusing for large graphs, the default value is FALSE.
}
  \item{layout}{
	A layouting function that determines the placement of the nodes in the graph. Please refer to the \code{\link[igraph]{layout}} manual entry in the \pkg{igraph} package for further details. By default, the Fruchterman-Reingold algorithm is used.
}

  \item{piecewise}{
  If set to true, a piecewise layout is used, i.e. the subgraphs corresponding to different basins of attraction are separated and layouted separately.  
}

  \item{basin.lty}{
	The line type used for edges in a basin of attraction. Defaults to 2 (dashed).
}
  \item{attractor.lty}{
	If \code{highlightAttractors} is true, this specifies the line type for edges in an attractor. Defaults to 1 (straight).
}
  \item{plotIt}{
	If this is true, a plot is generated. Otherwise, only an object of class \code{igraph} is returned, but no plot is drawn.
}
  \item{\dots}{
	Further graphical parameters to be passed to \code{\link[igraph:plot.graph]{plot.igraph}}.
}
}
\details{
	This function uses the \code{\link[igraph:plot.graph]{plot.igraph}} function from the \pkg{igraph} package. The plots are customizeable using the \code{\dots} argument. For details on possible parameters, please refer to \code{\link[igraph:plot.common]{igraph.plotting}}.
}
\value{
	Returns an invisible object of class \code{igraph} containing the state graph, including color and line attributes.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, \code{\link{getStateSummary}}, \code{\link{plotNetworkWiring}}, \code{\link[igraph:plot.common]{igraph.plotting}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# plot state graph
\dontrun{
plotStateGraph(attractors, main="Cell cycle network")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph
	plot}
