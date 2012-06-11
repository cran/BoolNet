\name{getPathToAttractor}
\alias{getPathToAttractor}
\title{
Get state transitions between a state and its attractor
}
\description{
Lists the states in the path from a specified state to the corresponding synchronous attractor.
}
\usage{
getPathToAttractor(network, 
                   state,
                   includeAttractorStates = c("all","first","none"))
}
\arguments{
  \item{network}{Either a structure of class \code{BooleanNetwork}, or a structure of class \code{AttractorInfo}. In the former case, a synchronous attractor search starting from \code{state} is conducted. In the latter case, \code{network} must be the result of a call to \code{\link{getAttractors}} with \code{returnTable=TRUE}, and its transition table must include \code{state}.
    
}
  \item{state}{
  A binary vector with exactly one entry per gene in the network.
}
  \item{includeAttractorStates}{
  Specifies whether the actual attractor states are included in the resulting table or not. If \code{includeAttractorStates = "all"} (which is the default behaviour), the sequence ends when the attractor was traversed once. If \code{includeAttractorStates = "first"}, only the first state of attractor is added to the sequence. This corresponds to the behaviour prior to \pkg{BoolNet} version 1.5. If {includeAttractorStates = "none"}, the sequence ends with the last non-attractor state. In this case, the sequence can be empty if the start state is an attractor state.
}
}
\value{
Returns a data frame with the genes in the columns. The rows are the successive states from \code{state} to the  the corresponding attractor. Depending on \code{includeAttractorStates}, attractor states are included or not.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}
}
\examples{
# load example network
data(cellcycle)

# get path from a state to its attractor
# include all attractor states
path <- getPathToAttractor(cellcycle, rep(1,10), 
                           includeAttractorStates="all")
print(path)

# include only the first attractor state
path <- getPathToAttractor(cellcycle, rep(1,10), 
                           includeAttractorStates="first")
print(path)

# exclude attractor states
path <- getPathToAttractor(cellcycle, rep(1,10), 
                           includeAttractorStates="none")
print(path)
}

\keyword{Boolean network
        sequence
        path
        basin
        attractor}
