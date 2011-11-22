\name{print.MarkovSimulation}
\alias{print.MarkovSimulation}
\title{
Print the results of a Markov chain simulation
}
\description{
A specialized method to print an object of class \code{MarkovSimulation}. This prints all states that have a non-zero probability to be reached after the number of iterations in the Markov simulation. If the simulation was run with \code{returnTable=TRUE}, it also prints a table of state transitions and their probabilities to be chosen in a single step.
}
\usage{
\method{print}{MarkovSimulation}(x,
      activeOnly = FALSE,
      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{
An object of class \code{MarkovSimulation} to be printed
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
\code{\link{print}}, \code{\link{markovSimulation}}
}

\keyword{print
        Markov chain simulation
        probabilistic Boolean network
        PBN}% __ONLY ONE__ keyword per line
