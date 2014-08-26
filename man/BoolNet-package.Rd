\name{BoolNet-package}
\Rdversion{1.1}
\alias{BoolNet-package}
\alias{BoolNet}
\docType{package}
\title{Analyze Boolean networks}
\description{Tools for reconstruction, analysis and visualization of synchronous, asynchronous and probabilistic Boolean networks, in particular for the identification of attractors in gene-regulatory networks}
\details{
\tabular{ll}{
Package: \tab BoolNet\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2014-08-26\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
This package provides useful methods for the construction and generation of Boolean networks and for their analysis. In particular, it is designed for the analysis of gene-regulatory networks. The software supports four types of networks:
\describe{

\item{Synchronous Boolean networks}{These networks consist of a set of Boolean variables (genes) \var{X} and a set of transition functions, one for each variable. These transition functions map an input from the set \var{X} to a Boolean value. A state is a vector of values for each of the variables in \var{X}. Then, the next state of the network is calculated by applying \emph{all} transition functions to the state.}

\item{Asynchronous Boolean networks}{Asynchronous networks have the same structure as synchronous Boolean networks. Yet, the next state of the network is calculating by choosing only \emph{one} of the transition functions at random and updating the corresponding Boolean variable (gene). This corresponds to the assumption that in a genetic network, gene expression levels are likely to change at different points of time.}

\item{Synchronous Boolean networks with time delays}{These networks additionally include dependencies on genes at time steps other than the previous time step. That is, not only the immediate predecessor state is considered to determine the next state of the network, but earlier states can be considered as well. Furthermore, it is possible to use predicates that depend on the absolute time point, i.e. the number of transitions from an initial state. }

\item{Probabilistic Boolean networks}{Probabilistic networks allow for specifying more than one transition function per variable/gene. Each of these functions has a probability to be chosen, where the probabilities of all functions for one variable sum up to 1. Transitions are performed synchronously by choosing one transition function for each gene according to their probabilities and applying them to the current state.}
}

Networks can be assembled in several ways using \pkg{BoolNet}: The \code{\link{reconstructNetwork}} function infers Boolean networks from time series of measurements using several popular reconstruction algorithms. \code{\link{binarizeTimeSeries}} provides a means of binarizing real-valued time series for these reconstruction algorithms.
Boolean networks (synchronous, asynchronous, and probabilistic networks) can also be expressed in a description language and loaded from files using \code{\link{loadNetwork}} or stored to files using \code{\link{saveNetwork}}. Furthermore, networks can be imported from BioTapestry using \code{\link{loadBioTapestry}} and from SBML with the \code{sbml-qual} package using \code{\link{loadSBML}}. The package also includes an export to SBML (see \code{\link{toSBML}}).

Via \code{\link{generateRandomNKNetwork}} and \code{\link{perturbNetwork}}, the package supports various methods of generating random networks and perturbing existing networks for robustness analysis.

The \code{\link{getAttractors}} function identifies attractor cycles in a synchronous or asynchronous Boolean network. Attractors can be identified by exhaustive search or heuristic methods. For networks with time delays, the function \code{simulateSymbolicModel} simulates the model and identifies attractors.

The \code{\link{markovSimulation}} function identifies relevant states in probabilistic Boolean networks by performing a Markov chain simulation.

The package also provides methods to visualize state transitions and basins of attraction \if{latex}{\cr}
(\code{\link{plotPBNTransitions}}, \code{\link{plotStateGraph}}), to plot the wiring of a network \if{latex}{\cr}(\code{\link{plotNetworkWiring}}), to plot attractor state tables (\code{\link{plotAttractors}}) and sequences of states (\code{\link{plotSequence}}), and to export them to LaTeX (\code{\link{attractorsToLaTeX}} and \code{\link{sequenceToLaTeX}}) and Pajek (\code{\link{toPajek}}). 

Transition tables of the network can be analyzed using \code{\link{getTransitionTable}}. Paths from start states to their corresponding attractors are identified using \code{\link{getPathToAttractor}}.
}
\author{
Christoph Müssel, Martin Hopfensitz, Dao Zhou, Hans Kestler

Maintainer: Hans Kestler <hans.kestler@uni-ulm.de>

}
\references{
S. A. Kauffman (1969), Metabolic stability and epigenesis in randomly constructed nets. J. Theor. Biol. 22:437--467.

S. A. Kauffman (1993), The Origins of Order. Oxford University Press.

Further references are listed in the corresponding help sections.
}
\keyword{package 
	Boolean network
	probabilistic Boolean network
	PBN
	attractor
	cycle
	basin
	transition
	Markov chain simulation
	state
	graph
	reconstruction}

\examples{
##################################
# Example 1: identify attractors #
################################## 

# load example data
data(cellcycle)

# get all synchronous attractors by exhaustive search
attractors <- getAttractors(cellcycle)

# plot attractors side by side
par(mfrow=c(2,length(attractors$attractors)))
plotAttractors(attractors)

# identifies asynchronous attractors
attractors <- getAttractors(cellcycle, 
                            type="asynchronous", startStates=100)
\donttest{plotAttractors(attractors, mode="graph")}

####################################
# Example 2: reconstruct a network #
####################################

# load example data
data(yeastTimeSeries)

# perform binarization with k-means
bin <- binarizeTimeSeries(yeastTimeSeries)

# reconstruct networks from transition table
net <- reconstructNetwork(bin$binarizedMeasurements, 
                          method="bestfit", maxK=3, returnPBN=TRUE)

# analyze the network using a Markov chain simulation
print(markovSimulation(net, returnTable=FALSE))
}
