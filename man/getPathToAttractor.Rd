\name{getPathToAttractor}
\alias{getPathToAttractor}
\title{
Get state transitions between a state and its attractor
}
\description{
Lists the states in the path from a specified state to the corresponding synchronous attractor.
}
\usage{
getPathToAttractor(network, state)
}
\arguments{
  \item{network}{Either a structure of class \code{BooleanNetwork}, or a structure of class \code{AttractorInfo}. In the former case, a synchronous attractor search starting from \code{state} is conducted. In the latter case, \code{network} must be the result of a call to \code{\link{getAttractors}} with \code{returnTable=TRUE}, and its transition table must include \code{state}.
    
}
  \item{state}{
  A binary vector with exactly one entry per gene in the network.
}
}
\value{
Returns a data frame with the genes in the columns. The rows are the successive states from \code{state} to the first state of the corresponding attractor.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}
}
\examples{
library(BoolNet)

# load example network
data(cellcycle)

# get path from a state to its attractor
path <- getPathToAttractor(cellcycle, rep(0,10))
print(path)
}

\keyword{Boolean network
        path
        basin
        attractor}
