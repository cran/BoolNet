\name{print.TransitionTable}
\Rdversion{1.1}
\alias{print.TransitionTable}
\alias{print.BooleanStateInfo}
\title{
Print a transition table
}
\description{
Specialized print method to print a transition table with the initial state in the first column, the successor state in the second column, the basin of attraction to which the state leads in the third column, and the number of transitions to the attractor in the fourth column.
}
\usage{
print.TransitionTable(x, quote = FALSE, max.levels = NULL, width = getOption("width"), ...)

print.BooleanStateInfo(x, quote = FALSE, max.levels = NULL, width = getOption("width"), ...)
}
\arguments{
  \item{x}{
An object of class \code{TransitionTable} or \code{BooleanStateInfo} to be printed
}
  \item{quote}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{max.levels}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{width}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{\dots}{
Further parameters for the \code{\link{print}} method.
}
}
\value{
Invisibly returns the printed object
}

\seealso{
\code{\link{print}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, \code{\link{getStateSummary}}
}
\keyword{transition table
	print}