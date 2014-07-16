\name{toPajek}
\Rdversion{1.1}
\alias{toPajek}
\title{
Export a network to the Pajek file format
}
\description{
Exports a network to the Pajek file format to visualize transition trajectories. For more information on Pajek, please refer to \url{http://pajek.imfm.si/doku.php}
}
\usage{
toPajek(stateGraph, file = "boolean.net", includeLabels=FALSE, ...)
}
\arguments{
  \item{stateGraph}{
	An object of class \code{AttractorInfo} or \code{SymbolicSimulation}, as returned by \code{\link{getAttractors}} and \code{\link{simulateSymbolicModel}} respectively. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE. Similarly, \code{simulateSymbolicModel} must be called with \code{returnGraph=TRUE}. Alternatively, \code{stateGraph} can be an object of class \code{TransitionTable}, which can be extracted using the functions \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, or \code{\link{getStateSummary}} .
}
  \item{file}{
	The name of the output file for Pajek. Defaults to "boolean.net".
}
 \item{includeLabels}{
	If set to true, the vertices of the graph in the output file are labeled with the binary encodings of the states. Defaults to FALSE.
}

 \item{\dots}{
This is only for compatibility with previous versions and should not be used. 
}
}
\value{
This function has no return value.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, \code{\link{getStateSummary}}, \code{\link{toSBML}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# export to Pajek
toPajek(attractors, file="pajek_export.net")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pajek
	Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph}
