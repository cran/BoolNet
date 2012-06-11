\name{getTransitionProbabilities}
\alias{getTransitionProbabilities}
\title{
Get a matrix of transitions and their probabilities in probabilistic Boolean networks
}
\description{
Retrieves the state transitions and their probabilities in a probabilistic Boolean network. This takes the transition table information calculated by the \code{\link{markovSimulation}} method. 
}
\usage{
getTransitionProbabilities(markovSimulation)
}

\arguments{
  \item{markovSimulation}{
  An object of class \code{MarkovSimulation}, as returned by \code{\link{markovSimulation}}. As the transition table information in this structure is required, \code{\link{markovSimulation}} must be called with \code{returnTable} set to TRUE.
}
}

\value{
Returns a data frame with the first \code{n} columns describing the values of the genes before the transition, the next \code{n} columns describing the values of the genes after the transition, and the last column containing the probability of the transition. Here, \code{n} is the number of genes in the underlying network. Only transitions with non-zero probability are included.
}

\seealso{
\code{\link{markovSimulation}}
}
\examples{
# load example network
data(examplePBN)

# perform a Markov chain simulation
sim <- markovSimulation(examplePBN)

# print out the probability table
print(getTransitionProbabilities(sim))
}

\keyword{Boolean network
        probabilistic Boolean network
        PBN
        transition
        probability}
