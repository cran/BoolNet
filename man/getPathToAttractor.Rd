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
  \item{network}{Either a network structure of class \code{BooleanNetwork} or \code{SymbolicBooleanNetwork} , or an attractor search result of class \code{AttractorInfo}. In the former case, a synchronous attractor search starting from \code{state} is conducted. In the latter case, \code{network} must be the result of a call to \code{\link{getAttractors}} with \code{returnTable=TRUE}, and its transition table must include \code{state}.    
}

  \item{state}{
  A binary vector with exactly one entry per gene in the network. If \code{network} is of class \code{SymbolicBooleanNetwork} and makes use of more than one predecessor state, this can also be a matrix with the genes in the columns and multiple predecessor states in the rows.
}

  \item{includeAttractorStates}{
  Specifies whether the actual attractor states are included in the resulting table or not. If \code{includeAttractorStates = "all"} (which is the default behaviour), the sequence ends when the attractor was traversed once. If \code{includeAttractorStates = "first"}, only the first state of attractor is added to the sequence. This corresponds to the behaviour prior to \pkg{BoolNet} version 1.5. If {includeAttractorStates = "none"}, the sequence ends with the last non-attractor state. In this case, the sequence can be empty if the start state is an attractor state.
}
}
\value{
Returns a data frame with the genes in the columns. The rows are the successive states from \code{state} to the  the corresponding attractor. Depending on \code{includeAttractorStates}, attractor states are included or not. The data frame has an attribute \code{attractor} specifying the indices of the states that belong to the attractor. If \code{includeAttractorStates} is \code{"first"} or \code{"none"}, these indices may correspond to states that are not included in the sequence itself. This attribute is used by \code{\link{plotSequence}} to highlight the attractor states.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, \code{\link{plotSequence}}, \code{\link{attributes}}
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
