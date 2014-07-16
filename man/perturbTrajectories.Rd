\name{perturbTrajectories}
\alias{perturbTrajectories}

\title{
Perturb the state trajectories and calculate robustness measures
}
\description{
Perturbs the state trajectories of a network and assesses the robustness by comparing the successor states or the attractors of a set of initial states and a set of perturbed copies of these initial states.
}
\usage{
perturbTrajectories(network, 
                    measure = c("hamming", "sensitivity", "attractor"), 
                    numSamples = 1000, 
                    flipBits = 1, 
                    updateType = c("synchronous", "asynchronous", "probabilistic"), 
                    gene, 
                    ...)
}
\arguments{
  \item{network}{
A network structure of class \code{BooleanNetwork}, \code{SymbolicBooleanNetwork} or \code{ProbabilisticBooleanNetwork} whose robustness is measured.
}
  \item{measure}{
Defines the way the robustness is measured (see Details).
}
  \item{numSamples}{
The number of randomly generated pairs of initial states and perturbed copies. Defaults to 1000.
}
  \item{flipBits}{
The number of bits that are flipped to generate a perturbed copy of an initial state. Defaults to 1.
}
  \item{updateType}{
If \code{measure="hamming"}, the type of update that is performed to calculate successor states.
}
  \item{gene}{
If \code{measure="sensitivity"}, the name or index of the gene for whose transition function the average sensitivity is calculated.
}
  \item{\dots}{
Further parameters to \code{\link{stateTransition}} and \code{\link{getAttractors}}.
}
}

\details{
The function generates a set of \code{numSamples} initial states and then applies \code{flipBits} random bit flips to each initial state to generate a perturbed copy of each initial state. For each pair of initial state and perturbed state, a robustness statistic is calculated depending \code{measure}:

If \code{measure="hamming"}, the normalized Hamming distances between the successor states of each initial state and the corresponding perturbed state are calculated. 

If \code{measure="sensitivity"}, the average sensitivity of a specific transition function (specified in the \code{gene} parameter) is approximated: The statistic is a logical vector that is \code{TRUE} if \code{gene} differs in the successor states of each initial state and the corresponding perturbed state.

If \code{measure="attractor"}, the attractors of all initial states and all perturbed states are identified. The statistic is a logical vector specifying whether the attractors are identical in each pair of initial state and perturbed initial state.
}
\value{
A list with the following items:
 \item{stat}{A vector of size \code{numSamples} containing the robustness statistic for each pair of initial state and perturbed copy.}
 \item{value}{The summarized statistic (i.e. the mean value) over all state pairs.}
}
\references{
I. Shmulevich and S. A. Kauffman (2004), Activities and Sensitivities in Boolean Network Models. Physical Review Letters 93(4):048701.
}

\seealso{
\code{\link{testNetworkProperties}}, \code{\link{perturbNetwork}}
}
\examples{
data(cellcycle)

# calculate average normalized Hamming distance of successor states
hamming <- perturbTrajectories(cellcycle, measure="hamming", numSamples=100)
print(hamming$value)

# calculate average sensitivity of transition function for gene "Cdh1"
sensitivity <- perturbTrajectories(cellcycle, measure="sensitivity", numSamples=100, gene="Cdh1")
print(sensitivity$value)

# calculate percentage of equal attractors for state pairs
attrEqual <- perturbTrajectories(cellcycle, measure="attractor", numSamples=100)
print(attrEqual$value)
}
\keyword{Boolean network
         probabilistic Boolean network,
         PBN,
         perturbation
         perturb
         trajectory
         noise
         robustness}
