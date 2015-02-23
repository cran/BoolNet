\name{markovSimulation}
\alias{markovSimulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify important states in probabilistic Boolean networks
}
\description{
Identifies important states in probabilistic Boolean networks (PBN) using a Markov chain simulation
}
\usage{
markovSimulation(network, 
                 numIterations = 1000, 
                 startStates = list(), 
                 cutoff = 0.001, 
                 returnTable = TRUE)
}

\arguments{
  \item{network}{
  An object of class \code{ProbabilisticBooleanNetwork} or \code{BooleanNetwork} whose transitions are simulated
}
  \item{numIterations}{
  The number of iterations for the matrix multiplication, which corresponds to the number of state transitions to simulate
}
  \item{startStates}{
  An optional list of start states. Each entry of the list must be a vector with a 0/1 value for each gene. If specified, the simulation is restricted to the states reachable from the supplied start states. Otherwise, all states are considered.
}
  \item{cutoff}{
  The cutoff value used to determine if a probability is 0. All output probabilities less than or equal to this value are set to 0.
}
  \item{returnTable}{
  If set to true, a transition table annotated with the probabilities for the transitions is included in the results. This is required by \code{\link{plotPBNTransitions}} and \code{\link{getTransitionProbabilities}}.
}
}
\details{
The algorithm identifies important states by performing the following steps: First, a Markov matrix is calculated from the set of transition functions, where each entry of the matrix specifies the probability of a state transition from the state belonging to the corresponding row to the state belonging to the corresponding column. A vector is initialized with uniform probability for all states (or -- if specified -- uniform probability for all start states) and repeatedly multiplied with the Markov matrix. The method returns all states with non-zero probability in this vector.
See the references for more details.
}
\value{
An object of class \code{MarkovSimulation} with the following components:
\item{reachedStates}{A data frame with one state in each row. The first columns specify the gene values of the state, and the last column holds the probability that the corresponding state is reached after \code{numIterations} transitions. Only states with a probability greater than \code{cutoff} are included in this table.}
\item{genes}{A vector of gene names of the input network}
\item{table}{If \code{returnTable=TRUE}, this structure holds a table of transitions with the corresponding probabilities that transitions are chosen. This is a list with the following components:
\describe{
    \item{initialStates}{A matrix of encoded start states of the transitions}
    \item{nextStates}{The encoded successor states of the transitions}
    \item{probabilities}{The probabilities that the transitions are chosen in a single step}
}}
}

\references{
I. Shmulevich, E. R. Dougherty, S. Kim, W. Zhang (2002), Probabilistic Boolean networks: a rule-based uncertainty model for gene regulatory networks. Bioinformatics 18(2):261--274.
}

\seealso{
\code{\link{reconstructNetwork}}, \code{\link{plotPBNTransitions}}, \code{\link{getTransitionProbabilities}}
}
\examples{
# load example network
data(examplePBN)

# perform a Markov chain simulation
sim <- markovSimulation(examplePBN)

# print the relevant states and transition probabilities
print(sim)

# plot the transitions and their probabilities
plotPBNTransitions(sim)
}
\keyword{probabilistic Boolean network
         PBN
         Markov chain
}
