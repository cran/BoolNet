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
\method{print}{TransitionTable}(x,
      activeOnly = FALSE, 
      ...)

\method{print}{BooleanStateInfo}(x,
      activeOnly=FALSE, 
      ...)
}
\arguments{

\item{x}{
An object of class \code{TransitionTable} or \code{BooleanStateInfo} to be printed
}

\item{activeOnly}{
  If set to true, a state is represented by a list of active genes (i.e., genes which are set to 1). If set to false, a state is represented by a binary vector with one entry for each gene, specifying whether the gene is active or not. Defaults to \code{FALSE}. 
}

\item{\dots}{
Further parameters for the \code{\link{print}} method. Currently not used.
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
