\name{getAttractorSequence}
\alias{getAttractorSequence}
\title{
Decode the state sequence of a synchronous attractor
}
\description{
Obtains the sequence of states belonging to a single synchronous attractor from the encoded data in an \code{AttractorInfo} structure.
}
\usage{
getAttractorSequence(attractorInfo, attractorNo)
}
\arguments{
  \item{attractorInfo}{
An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}. \code{getAttractors} must be called in synchronous mode.
}
  \item{attractorNo}{
The index of the attractor in \code{attractorInfo} whose state sequence should be obtained
}
}

\value{
Returns a data frame with the genes in the columns. The rows are the successive states of the attractor. The successor state of the last state (i.e. the last row) is the first state (i.e. the first row). 
}

\seealso{
\code{\link{getAttractors}}, \code{\link{getPathToAttractor}}, \code{\link{plotSequence}}, \code{\link{sequenceToLaTeX}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# print basin of 7-state attractor
print(getAttractorSequence(attractors, 2))

}

\keyword{Boolean network
	attractor
	cycle
	sequence
	state}
