\name{generateState}
\alias{generateState}

\title{
Generate a state vector from single gene values
}
\description{
This function provides a simple interface to generate full state vectors by specifying only the genes of interest. For example, only those genes that are active can be specified, while the others are set to a default value.
}
\usage{
generateState(network, 
              specs, 
              default = 0)
}
\arguments{
  \item{network}{
An network of class \code{BooleanNetwork}, \code{SymbolicBooleanNetwork} or \code{ProbabilisticBooleanNetwork} for which a state is generated.
}
  \item{specs}{
A named vector or list specifying the genes to be set. Here, the names of the elements correspond to the gene names, and the elements correspond to the gene values. The function can also generate a matrix of states if the elements of \code{specs} are vectors of values (of the same length).
}
  \item{default}{
The default value used for the unspecified genes (usually 0).
}
}
\value{
Returns a full state vector with one entry for each gene of the network, or a matrix with one state in each row if \code{specs} contains vectors of state values.
}

\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{stateTransition}}
}
\examples{
# load cell cycle network
data(cellcycle)

# generate a state in which only CycD and CycA are active
state <- generateState(cellcycle, c("CycD"=1, "CycA"=1))
print(state)

# use the state as a start state for attractor search
print(getAttractors(cellcycle, startStates=list(state)))
}
\keyword{state 
         Boolean network
         probabilistic Boolean network
         PBN}
