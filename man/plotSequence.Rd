\name{plotSequence}
\alias{plotSequence}
\title{
Plot a sequence of states
}
\description{
Visualizes sequences of states in synchronous Boolean networks, either by drawing a table of the involved states in two colors, or by drawing a graph of transitions between the successive states.
}
\usage{
plotSequence(network, 
             startState, 
             includeAttractorStates = c("all","first","none"), 
             sequence,
             title = "", 
             mode=c("table","graph"),
             plotFixed = TRUE, grouping = list(),
             onColor="#4daf4a", 
             offColor = "#e41a1c",
             layout, 
             drawLabels=TRUE, 
             drawLegend=TRUE, 
             highlightAttractor=TRUE,
             reverse = FALSE,
             borderColor = "black", 
             eps=0.1,
             attractor.sep.lwd = 2, 
             attractor.sep.col = "blue",
             ...)
}

\arguments{

  \item{network}{
	An object of class \code{BooleanNetwork} or \code{SymbolicBooleanNetwork} for which a sequence of state transitions is calculated
}

  \item{startState}{
  The start state of the sequence
  }
  
  \item{includeAttractorStates}{
  Specifies whether the actual attractor states are included in the plot or not (see also \code{\link{getPathToAttractor}}). If \code{includeAttractorStates = "all"} (which is the default behaviour), the sequence ends when the attractor was traversed once. If \code{includeAttractorStates = "first"}, only the first state of attractor is added to the sequence. If {includeAttractorStates = "none"}, the sequence ends with the last non-attractor state. 
  }
  
  \item{sequence}{
  The alternative call to \code{plotSequence} requires the specification of the sequence itself instead of the network and the start state. The sequence must be provided as a data frame with the genes in the columns and the successive states in the rows. For example, sequences can be obtained using \code{\link{getPathToAttractor}} or \code{\link{getAttractorSequence}} (however, the specialized plot \code{\link{plotAttractors}} exists for attractors).}
  
  \item{title}{
	An optional title for the plot
}
  
  \item{mode}{Switches between two kinds of attractor plots. See Details for more information. Default is "table".}
  
  \item{plotFixed}{This optional parameter is only used if \code{mode="table"}.
	If this is true, genes with fixed values are included in the plot. Otherwise, these genes are not drawn.
}
  
  \item{grouping}{This optional parameter is only used if \code{mode="table"} and specifies a structure to form groups of genes in the plot. This is a list with the following elements: 
 \describe{
  	\item{class}{A vector of names for the groups. These names will be printed in the region belonging to the group in the plot.}  	
  	\item{index}{A list with the same length as \code{class}. Each element is a vector of gene names or gene indices belonging to the group.}}
}

  \item{onColor}{This optional parameter is only used if \code{mode="table"} and specifies the color value for the 1/ON values in the table. Defaults to green.
}
  \item{offColor}{This optional parameter is only used if \code{mode="table"} and specifies the color value for the 0/OFF values in the table. Defaults to red.
}

  \item{layout}{If \code{mode="graph"}, this parameter specifies a layouting function that determines the placement of the nodes in the graph. Please refer to the \code{\link[igraph]{layout}} manual entry in the \pkg{igraph} package for further details. By default, the nodes are placed in a horizontal line.
}

  \item{drawLabels}{This parameter is only relevant if \code{mode="graph"}. It determines whether the nodes of the graph are annotated with the corresponding values of the genes in the attractor states.
}

  \item{drawLegend}{Specifies whether a color key for the ON/OFF states is drawn if \code{mode="table"}. Defaults to \code{TRUE}.
}

  \item{highlightAttractor}{
If set to true, the attractor states are highlighted in the plot. If \code{mode="table"}, a line is drawn at the begin of the attractor, and the states are labeled correspondingly. If \code{mode="graph"}, the attractor transitions are drawn as bold lines. Information on the attractor must be supplied in the attribute \code{attractor} of the sequence, which is a vector of indices of the states that belong to the attractor. This attribute is usually present if the sequence was obtained using \code{\link{getPathToAttractor}}.
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
  \item{attractor.sep.lwd}{
Specifies the line width of the attractor separator. Defaults to \code{2}.
}
  \item{attractor.sep.col}{
Specifies the line color of the attractor separator. Defaults to \code{"blue"}.
}

    \item{\dots}{
Further graphical parameters to be passed to \code{\link[igraph:plot.graph]{plot.igraph}} if \code{mode="graph"}.
}
    
}
\details{
This function comprises two different types of plots: 

The "table" mode visualizes the gene values of the states in the sequence. The figure is a table with the genes in the rows and the successive states of the sequence in the columns. Cells of the table are (by default) red for 0/OFF values and green for 1/ON values. If \code{grouping} is set, the genes are rearranged according to the indices in the group, horizontal separation lines are plotted between the groups, and the group names are printed.

The "graph" mode visualizes the transitions between different states. It creates a graph in which the vertices are the states in the sequence and the edges are state transitions among these states.

The function can be called with different types of inputs: The user can specify the parameters \code{network}, \code{startState} and \code{includeAttractorStates}), in which case \code{\link{getPathToAttractor}} is called to obtain the sequence. Alternatively, the sequence can be supplied directly as a data frame in the \code{sequence} parameter.
}
\value{
If \code{mode="table"}, a matrix corresponding to the table is returned. The matrix has the genes in the rows and the states of the attractors in the columns. If \code{sequence} was supplied, this corresponds to the transposed input whose rows may be rearranged if \code{grouping} was set.

If \code{mode="graph"}, an object of class \code{igraph} describing the graph for the sequence is returned.  
}

\seealso{
\code{\link{sequenceToLaTeX}}, \code{\link{plotAttractors}}, \code{\link{attractorsToLaTeX}}, \code{\link{getPathToAttractor}}, \code{\link{getAttractorSequence}}, \code{\link{simulateSymbolicModel}}
}
\examples{
# load example data
data(cellcycle)

# alternative 1: supply network and start state
# and plot sequence as a table
plotSequence(network=cellcycle, 
             startState=rep(1,10), 
             includeAttractorStates="all")

# alternative 2: calculate sequence in advance
sequence <- getPathToAttractor(cellcycle, 
                               state=rep(1,10), 
                               includeAttractorStates="all")

# plot sequence as a graph
plotSequence(sequence=sequence, 
             mode="graph")
}
\keyword{Boolean network
         state
         path
         sequence
         plot}
