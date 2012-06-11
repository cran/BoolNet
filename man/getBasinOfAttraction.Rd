\name{getBasinOfAttraction}
\Rdversion{1.1}
\alias{getBasinOfAttraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get states in basin of attraction
}
\description{
Extracts information on all states in the basin of a supplied attractor
}
\usage{
getBasinOfAttraction(attractorInfo, attractorNo)
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE.
}
  \item{attractorNo}{
	The index of the attractor in \code{attractorInfo} whose basin should be identified
}
}
\details{
	The function outputs a transition table containing only the states that are contained in the basin of attraction, and displays additional information on these states. If \code{attractorInfo} is the result of an exhaustive synchronous attractor search, the complete basin of attraction is returned. If \code{attractorInfo} is the result of a heuristic synchronous search, there is no guarantee that the complete basin of attraction is returned, as only the calculated states are included. Asynchronous search results are not supported, as no transition table is calculated.
}
\value{
Returns a generic dataframe of the class \code{TransitionTable}. For n genes, the first n columns code for the original state, i.e. each column represents the value of one gene. The next n columns code for the successive state after a transition. The column \code{attractorAssignment} indicates the attractor to the state is assigned (in this case, \code{attractorNo}), and the column \code{stepsToAttractor} indicates how many transitions are needed from the original state to the attractor.
The \code{TransitionTable} class supports pretty printing using the \code{\link{print}} method.
}
\seealso{
\code{\link{getStateSummary}}, \code{\link{getTransitionTable}}, \code{\link{getAttractors}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# print basin of first attractor
print(getBasinOfAttraction(attractors, 1))
}
\keyword{Boolean network
	attractor
	cycle
	basin}
