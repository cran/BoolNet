\name{getTransitionTable}
\Rdversion{1.1}
\alias{getTransitionTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve the transition table of a network
}
\description{
Retrieves the complete transition table and additional attractor information of a network
}
\usage{
getTransitionTable(attractorInfo)
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}
}
}

\value{
Returns a generic dataframe of the class \code{TransitionTable}. For n genes, the first n columns code for the original state (in this case, the \code{state} parameter), i.e. each column represents the value of one gene. The next n columns code for the successive state after a transition. The column \code{attractorAssignment} indicates the attractor to the state is assigned, and the column \code{stepsToAttractor} indicates how many transitions are needed from the original state to the attractor. The table has a row for each possible input state.
The \code{TransitionTable} class supports pretty printing using the \code{\link{print}} method.
}
\seealso{
\code{\link{getStateSummary}}, \code{\link{getBasinOfAttraction}}, \code{\link{getAttractors}}
}
\examples{
library(BoolNet)

# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(net)

# print the transition table
print(getTransitionTable(attractors))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Boolean network
	transition
	state}
