\name{simulateSymbolicModel}
\alias{simulateSymbolicModel}

\title{
Simulate a symbolic Boolean network
}
\description{
This function simulates Boolean networks in a symbolic representation, possibly with additional temporal qualifiers. The function can identify attractors, determine the state transition graph, and generate sequences of successive states.
}
\usage{
simulateSymbolicModel(network, 
              method = c("exhaustive", 
                         "random", 
                         "chosen", 
                         "sat.exhaustive", 
                         "sat.restricted"), 
              startStates = NULL,               
              returnSequences = 
                (!(match.arg(method) \%in\% 
              	  c("sat.exhaustive", "sat.restricted"))), 
              returnGraph = 
                  (!(match.arg(method) \%in\% 
              	    c("sat.exhaustive", "sat.restricted"))), 
              returnAttractors = TRUE, 
              maxTransitions = Inf,
              maxAttractorLength = Inf,
              canonical = TRUE)
}

\arguments{
  \item{network}{
A network structure of class \code{SymbolicBooleanNetwork}. These networks can be read from files by \code{\link{loadNetwork}}, \code{\link{loadBioTapestry}} or \code{\link{loadSBML}} with the \code{symbolic=TRUE} flag.
}
  \item{startStates}{
An optional parameter specifying the start states. If this is an integer value, it denotes the number of random start states to generate. Otherwise, it has to be a list of states. The list elements must either be vectors with one value for each gene in the network, or matrices with the genes in the columns and multiple predecessor states in the rows. These predecessor states may be evaluated if temporal predicates in the network have a time delay of more than one. If the number of supplied predecessor states is smaller than the maximum time delay in the network, genes are assumed to have had the same value as in the first supplied state prior to this state. In particular, if only a single state is supplied, it is assumed that the network always resided in this state prior to starting the simulation.
}
  \item{method}{
The simulation method to be used (see details). If \code{method} is not specified, the desired method is inferred from the type of \code{startStates}.
}
  \item{returnSequences}{
If set to true (and no SAT-based method is chosen), the return value has an element \code{sequences} specifying the sequences of states to the attractor.
}
  \item{returnGraph}{
If set to true (and no SAT-based method is chosen), the return value has an element \code{graph} specifying the state transition graph of the network.
}
  \item{returnAttractors}{
If set to true, the return value has an element \code{attractors} containing a list of identified attractors.
}
  \item{maxTransitions}{
The maximum number of state transitions to be performed for each start state (defaults to \code{Inf}).
}

\item{maxAttractorLength}{If \code{method="sat.restricted"}, this required parameter specifies the maximum size of attractors (i.e. the number of states in the loop) to be searched. For \code{method="sat.exhaustive"}, this parameter is optional and specifies the maximum attractor length for the initial length-restricted search phase that is performed to speed up the subsequent exhaustive search. In this case, changing this value might bring performance benefits, but does not change the results.}

 \item{canonical}{
	If set to true and \code{returnAttractors=TRUE}, the states in the attractors are rearranged such that the state whose binary encoding
makes up the smallest number is the first element of the vector. This ensures that attractors determined in runs with different start states are comparable, as the cycles may have been entered at different states.
}
}

\details{
Similarly to \code{\link{getAttractors}}, the symbolic simulator supports different simulation modes which can be specified in the \code{method} parameter:
\itemize{
\item{Exhaustive search} {If \code{method="exhaustive"}, all possible states in the network are used as start states. If the network has time delays greater than one (temporal network), this means that exhaustive search does not only cover all 2^n possible states for a network with n genes, but also all possible state histories of those genes for which longer delays are required.}
  	
\item{Heuristic search} {For \code{method="random"} or \code{method="chosen"}, a subset of states is used as start states for the simulation.

If \code{method="random"}, \code{startStates} is interpreted as an integer value specifying the number of states to be generated randomly. The algorithm is then initialized with these random start states.
  	
If \code{method="chosen"}, \code{startStates} is interpreted as a list of binary vectors, each specifying one start state (see also parameter description above for details). }

\item{SAT-based attractor search} { 	
If \code{method} is "sat.exhaustive" or "sat.restricted", the simulator transforms the network into a satisfiability problem and solves it using Armin Biere's PicoSAT solver (see also \code{\link{getAttractors}} for more details). If \code{method="sat.restricted"}, only attractors comprising up to \code{maxAttractorLength} states are identified. Otherwise, the algorithm by Dubrova and Teslenko is applied to identify all attractors. As the SAT-based approaches identify attractors directly, no state sequences and no transition graph are returned.}
}
}

\value{
Returns a list of class \code{SymbolicSimulation} containing the simulation results:

If \code{returnSequences} is true and no SAT-based method was chosen, the list contains an element \code{sequences} consisting of a list of data frames, each representing the state transitions performed from one start state (denoted as time step 0) to the attractor. Here, the columns correspond to the genes in the network, and the rows correspond to the states. Apart from the immediate start state, the sequences may also contain the supplied or assumed predecessor states of the start state (marked by a negative time step t) if the network contains time delays greater than one. 

If \code{returnGraph} is true and no SAT-based method was chosen, the list contains an element \code{graph} of class \code{TransitionTable}. Each row of the table corresponds to one state transition from an initial state to a successor state, i.e. an edge in the state transition graph.

If \code{returnAttractors} is true, the list contains an element \code{attractors}, which itself is a list of data frames. Each data frame represents one unique attractor, where each column corresponds to a gene, and each row corresponds to one state in the attractor.


If both \code{returnSequences} and \code{returnAttractors} are true, there is an additional element \code{attractorAssignment}. This integer vector specifies the indices of the attractors to which the sequences lead.

The structure supports pretty printing using the \code{\link{print}} method.}

\references{
E. Dubrova, M. Teslenko (2011), A SAT-based algorithm for finding attractors in synchronous Boolean networks. IEEE/ACM Transactions on Computational Biology and Bioinformatics 8(5):1393--1399.

A. Biere (2008), PicoSAT Essentials. Journal on Satisfiability, Boolean Modeling and Computation 4:75-97.
}

\seealso{
\code{\link{loadNetwork}}, \code{\link{loadBioTapestry}}, \code{\link{loadSBML}}, \code{\link{getAttractors}}, \code{\link{plotAttractors}}, \code{\link{attractorsToLaTeX}}, \code{\link{getTransitionTable}}, \code{\link{getBasinOfAttraction}}, \code{\link{getAttractorSequence}}, \code{\link{getStateSummary}}, \code{\link{getPathToAttractor}}, \code{\link{fixGenes}}
}

\examples{
\dontrun{
data(igf)

# exhaustive state space simulation
sim <- simulateSymbolicModel(igf)
plotAttractors(sim)

# exhaustive attractor search using SAT solver
sim <- simulateSymbolicModel(igf, method="sat.exhaustive")
plotAttractors(sim)
}
}