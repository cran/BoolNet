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
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}
}
  \item{attractorNo}{
	The index of the attractor in \code{attractorInfo} whose basin should be identified
}
}
\details{
	The function outputs a transition table containing only the states that are contained in the basin of attraction, and displays additional information on these states.
}
\value{
Returns a generic dataframe of the class \code{TransitionTable}. For n genes, the first n columns code for the original state, i.e. each column represents the value of one gene. The next n columns code for the successive state after a transition. The column \code{attractorAssignment} indicates the attractor to the state is assigned (in this case, \code{attractorNo}), and the column \code{stepsToAttractor} indicates how many transitions are needed from the original state to the attractor.
The \code{TransitionTable} class supports pretty printing using the \code{\link{print}} method.
}
\seealso{
\code{\link{getStateSummary}}, \code{\link{getTransitionTable}}, \code{\link{getAttractors}}
}
\examples{
library(BoolNet)

# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(net)

# print basin of first attractor
print(getBasinOfAttraction(attractors,1))
}
\keyword{Boolean network
	attractor
	cycle
	basin}
