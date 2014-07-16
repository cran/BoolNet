\name{getTransitionTable}
\Rdversion{1.1}
\alias{getTransitionTable}
\title{
Retrieve the transition table of a network
}
\description{
Retrieves the transition table and additional attractor information of a network. 
}
\usage{
getTransitionTable(attractorInfo)
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}, or of class \code{SymbolicSimulation}, as returned by \code{\link{simulateSymbolicModel}}. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE. Similarly, \code{simulateSymbolicModel} must be called with \code{returnGraph=TRUE}.
}
}
\details{
Depending on the configuration of the call to \code{getAttractors} or \code{simulateSymbolicModel} that returned \code{attractorInfo}, this function either returns the complete transition table (for exhaustive synchronous search) or the part of the transition table calculated in a heuristic synchronous search. Asynchronous search is not supported, as no transition table is calculated.
}

\value{
Returns a generic dataframe of the class \code{TransitionTable}. For n genes, the first n columns code for the original state (in this case, the \code{state} parameter), i.e. each column represents the value of one gene. The next n columns code for the successive state after a transition. The column \code{attractorAssignment} indicates the attractor to the state is assigned. If this information is available,  the column \code{stepsToAttractor} indicates how many transitions are needed from the original state to the attractor. The table has a row for each possible input state.
The \code{TransitionTable} class supports pretty printing using the \code{\link{print}} method.
}
\seealso{
\code{\link{getStateSummary}}, \code{\link{getBasinOfAttraction}}, \code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# print the transition table
print(getTransitionTable(attractors))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Boolean network
	transition
	transition table
	state}
