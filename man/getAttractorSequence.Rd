\name{getAttractorSequence}
\alias{getAttractorSequence}
\title{
Decode the state sequence of a synchronous attractor
}
\description{
Obtains the sequence of states belonging to a single synchronous attractor from the encoded data in an \code{AttractorInfo} structure or in a \code{SymbolicSimulation} structure.
}
\usage{
getAttractorSequence(attractorInfo, attractorNo)
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}, or of class \code{SymbolicSimulation}, as returned by \code{\link{simulateSymbolicModel}}. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE. Similarly, \code{simulateSymbolicModel} must be called with \code{returnGraph=TRUE}.
}

  \item{attractorNo}{
The index of the attractor in \code{attractorInfo} whose state sequence should be obtained
}
}

\value{
Returns a data frame with the genes in the columns. The rows are the successive states of the attractor. The successor state of the last state (i.e. the last row) is the first state (i.e. the first row). 
}

\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{getPathToAttractor}}, \code{\link{plotSequence}}, \code{\link{sequenceToLaTeX}}
}
\examples{
\dontrun{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# print basin of 7-state attractor
print(getAttractorSequence(attractors, 2))
}
}
