\name{getStateSummary}
\Rdversion{1.1}
\alias{getStateSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve summary information on a state
}
\description{
Returns information on the supplied state, i.e. the successive state after a transition, the attractor to which the state leads, and the distance to this attractor.
}
\usage{
getStateSummary(attractorInfo, state)
}

\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}
}
  \item{state}{
	A 0-1 vector with n elements (where n is the number of genes in the underlying networks) describing the state.
}
}
\value{
Returns a generic dataframe of the class \code{TransitionTable}. For n genes, the first n columns code for the original state (in this case, the \code{state} parameter), i.e. each column represents the value of one gene. The next n columns code for the successive state after a transition. The column \code{attractorAssignment} indicates the attractor to the state is assigned, and the column \code{stepsToAttractor} indicates how many transitions are needed from the original state to the attractor. In this case, the table has only one row describing the supplied state.
The \code{TransitionTable} class supports pretty printing using the \code{\link{print}} method.
}
\seealso{
\code{\link{getBasinOfAttraction}}, \code{\link{getTransitionTable}}, \code{\link{getAttractors}}
}
\examples{
library(BoolNet)

# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(net)

# print information for an arbitrary state
print(getStateSummary(attractors,c(1,1,1,1,1,1,1,1,1,1)))
}
\keyword{Boolean network
	transition
	state}
