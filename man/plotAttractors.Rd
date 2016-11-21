\name{plotAttractors}
\Rdversion{1.1}
\alias{plotAttractors}

\title{
Plot state tables or transition graphs of attractors
}
\description{
Visualizes attractors, either by drawing a table of the involved states in two colors, or by drawing a graph of transitions between the states of the attractor.
}
\usage{
plotAttractors(attractorInfo, 
               subset, 
               title = "", 
               mode = c("table","graph"),
               grouping = list(), 
               plotFixed = TRUE, 
               onColor = "#4daf4a",
               offColor = "#e41a1c",
               layout = layout.circle, 
               drawLabels = TRUE,
               drawLegend = TRUE,
               ask = TRUE,
               reverse = FALSE,
               borderColor = "black",
               eps = 0.1, 
               allInOnePlot = FALSE,
               ...) 
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}, or an object of class \code{SymbolicSimulation}, as returned by \code{\link{simulateSymbolicModel}}.
}
  \item{subset}{An subset of attractors to be plotted. This is a vector of attractor indices in \code{attractorInfo}.}
  
   
  \item{title}{
	An optional title for the plot
}
  
  \item{mode}{Switches between two kinds of attractor plots. See Details for more information. Default is "table".}
  
  \item{grouping}{This optional parameter is only used if \code{mode="table"} and specifies a structure to form groups of genes in the plot. This is a list with the following elements: 
 \describe{
  	\item{class}{A vector of names for the groups. These names will be printed in the region belonging to the group in the plot.}  	
  	\item{index}{A list with the same length as \code{class}. Each element is a vector of gene names or gene indices belonging to the group.}}
}

  \item{plotFixed}{This optional parameter is only used if \code{mode="table"}.
	If this is true, genes with fixed values are included in the plot. Otherwise, these genes are not drawn.
}
  \item{onColor}{This optional parameter is only used if \code{mode="table"} and specifies the color value for the 1/ON values in the table. Defaults to greenish color.
}
  \item{offColor}{This optional parameter is only used if \code{mode="table"} and specifies the color value for the 0/OFF values in the table. Defaults to reddish color.
}

  \item{layout}{If \code{mode="graph"}, this parameter specifies a layouting function that determines the placement of the nodes in the graph. Please refer to the \code{\link[igraph]{layout}} manual entry in the \pkg{igraph} package for further details. By default, the circle layout is used.
}

  \item{drawLabels}{This parameter is only relevant if \code{mode="graph"}. It determines whether the nodes of the graph are annotated with the corresponding values of the genes in the attractor states.
}

  \item{drawLegend}{Specifies whether a color key for the ON/OFF states is drawn if \code{mode="table"}. Defaults to \code{TRUE}.
} 
 
   \item{ask}{
If set to true, the plot function will prompt for a user input for each new plot that is shown on an interactive device (see \code{link{par("ask")}}). 
}
  \item{reverse}{
Specifies the order of the genes in the plot. By default, the first gene is placed in the first row of the plot. If \code{reverse=TRUE} (which was the default until \pkg{BoolNet} version 2.0.2), the first gene in the network is placed in the bottom row of the plot. 
}
  \item{borderColor}{
Specifies the border or seprating color of states in an attractor. Defaults to \code{"black"}.
}
  \item{eps}{
Specifies plotting margin for the sequence of states. Defaults to \code{0.1}.
}
  \item{allInOnePlot}{
If this is \code{TRUE} then all attractors, with \code{mode = "table"}, are plotted in one plot as specified internally by \code{par$mfrow} parameter. Previous value of the \code{par$mfrow} parameter is preserved. Defaults to \code{FALSE}, meaning the plots for more than one attractor will be switched interactively or all plotted in an non-interactive graphical device.
}

    \item{\dots}{
Further graphical parameters to be passed to \code{\link[igraph:plot.graph]{plot.igraph}} if \code{mode="graph"}.
}

}
\details{
This function comprises two different types of plots: 

The "table" mode visualizes the gene values of the states in the attractor and is only suited for synchronous or steady-state attractors. Complex asynchronous attractors are omitted in this mode. Attractors in \code{attractorInfo} are first grouped by length. Then, a figure is plotted to the currently selected device for each attractor length (i.e. one plot with all attractors consisting of 1 state, one plot with all attractors consisting of 2 states, etc.). If \code{ask=TRUE} and the standard X11 output device is used, the user must confirm that the next plot for the next attractor size should be shown.
The figure is a table with the genes in the rows and the states of the attractors in the columns. Cells of the table are (by default) red for 0/OFF values and green for 1/ON values. If \code{grouping} is set, the genes are rearranged according to the indices in the group, horizontal separation lines are plotted between the groups, and the group names are printed.

The "graph" mode visualizes the transitions between different states. It creates a graph in which the vertices are the states in the attractor and the edges are state transitions among these states. This mode can visualize all kinds of attractors, including complex/loose attractors. One plot is drawn for each attractor. As before, this means that on the standard output device, only the last plot is displayed unless you set \code{par(mfrow=c(...))} accordingly.

}
\value{
If \code{mode="table"}, a list of matrices corresponding to the tables is returned. Each of these matrices has the genes in the rows and the states of the attractors in the columns.

If \code{mode="graph"}, a list of objects of class \code{igraph} is returned. Each of these objects describes the graph for one attractor.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{attractorsToLaTeX}}, \code{\link{plotSequence}}, \code{\link{sequenceToLaTeX}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# calculate number of different attractor lengths, 
# and plot attractors side by side in "table" mode
par(mfrow=c(1, length(table(sapply(attractors$attractors,
                          function(attractor)
                          {
                            length(attractor$involvedStates)
                          })))))
plotAttractors(attractors)

# plot attractors in "graph" mode
par(mfrow=c(1, length(attractors$attractors)))
plotAttractors(attractors, mode="graph")

# identify asynchronous attractors
attractors <- getAttractors(cellcycle, type="asynchronous")

# plot attractors in "graph" mode
par(mfrow=c(1, length(attractors$attractors)))
plotAttractors(attractors, mode="graph")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{
	Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph
	plot}
