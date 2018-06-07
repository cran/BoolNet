\name{NEWS}
\title{NEWS file for the \pkg{BoolNet} package}

\section{Changes in version 2.1.4}{
\itemize{
\item{\code{getAttractors()} now allows for plotting all attractors inside one plot. Fixed minor bugs.}
}
}

\section{Changes in version 2.1.1}{
\itemize{
\item{Fixed compilation errors on Windows.}
}
}

\section{Changes in version 2.1.0}{

\subsection{Modified functions}{
\itemize{
\item{\code{getAttractors()} and \code{simulateSymbolicModel()} now support the identification of attractors in large networks based on a formulation as a satisfiability (SAT) problem. }

\item{\code{plotAttractors()} and \code{plotSequence()} now plot the genes in the order of the network instead of the reverse order by default. This behaviour can be controlled using a new parameter \code{reverse} which is also available in \code{attractorsToLaTeX()} and \code{sequenceToLaTeX()}.}

\item{Bugfix regarding negated temporal predicates in \code{simulateSymbolicModel()}.}
}
}
}

\section{Changes in version 2.0.2}{

\subsection{Modified functions}{
\itemize{
\item{Fixed undefined behaviour in \code{markovSimulation()}.}
}
}
}

\section{Changes in version 2.0.1}{

\subsection{Modified functions}{
\itemize{
\item{Fixed memory misalignment in \code{simulateSymbolicModel()}.}

\item{\code{loadSBML()} now accepts nodes that are constant, but have no initial value.}
}
}
}

\section{Changes in version 2.0}{

\subsection{New functions}{
\itemize{
\item{Support of temporal networks, and inclusion of a new simulator \code{simulateSymbolicModel()} to simulate these networks. Related functions include \code{truthTableToSymbolic()} and \code{symbolicToTruthTable()} to convert networks between the symbolic representation used by the new simulator and the truth table representation employed by the standard simulator.}

\item{New function \code{perturbTrajectories()} to assess the robustness of networks to noise in the states}.}

}

\subsection{Modified functions}{
\itemize{

\item{\code{loadNetwork()} can now load networks in a symbolic representation and with temporal extensions. \code{loadSBML()} and \code{loadBioTapestry()} can load symbolic networks without temporal extensions. }
 
\item{Most functions of the package have been adapted to work either with the symbolic representation or with the truth table representation of networks.}

\item{\code{plotSequence()} and \code{sequenceToLaTeX()} now visualize the attractor.}

\item{\code{reconstructNetwork()} now supports the specification of prior knowledge in form of required or excluded dependencies. Furthermore, it can now reconstruct networks from perturbation time series. By default, the function now returns incomplete functions with "don't care" values" instead of enumerating all possible functions.
}

\item{\code{generateTimeSeries()} can now generate artificial perturbation data with simulated knock-out or overexpression of genes.}

\item{\code{generateRandomNKNetwork()} can now be supplied with a user-defined generation function for the transition functions. Generation functions \code{generateCanalyzing()} and \code{generateNestedCanalyzing()} for canalyzing functions and nested canalyzing functions are included in the package.}

\item{\code{testNetworkProperties()} supports several new tests that perturb the network states instead of the networks themselves. These are available in the test functions \code{testAttractorRobustness()} and \code{testTransitionRobustness()}.
}
}
}

}

\section{Changes in version 1.63}{

\subsection{Modified functions}{

\itemize{
\item{Fixed issues preventing the use of \pkg{BoolNet} on Big Endian systems.}

\item{Eliminated some bad style code.}

\item{Fixed some valgrind notes.}

}
}
}

\section{Changes in version 1.62}{
\subsection{Modified functions}{

\itemize{

\item{Minor bugfixes in \code{loadNetwork()}.}

\item{Fixed undefined behaviour warnings for GCC 4.9.}
}
}
}

\section{Changes in version 1.61}{
\subsection{Modified functions}{
\itemize{

\item{Bugfixes in \code{plotAttractors()} and \code{plotSequence()}.}

\item{Fixed compatibility issues with R 3.0 alpha.}

}
}
}

\section{Changes in version 1.60}{

\subsection{New functions}{

\itemize{

\item{Support for SBML: \code{loadSBML()} and \code{toSBML()} import from and export to SBML documents with the \code{sbml-qual} extension package.}

\item{\code{saveNetwork()} stores networks in the \pkg{BoolNet} file format.}

}
}

\subsection{Modified functions}{
\itemize{

\item{The DNF generator employed by \code{generateRandomNKNetwork()} and \code{simplifyNetwork()} (as well as by the new functions \code{saveNetwork()} and \code{toSBML()}) now supports minimizing the canonical DNFs.}

}
}
}

\section{Changes in version 1.51}{

\subsection{Modified functions}{
\itemize{

\item{\pkg{BoolNet} now supports the modified interface of \pkg{igraph} 0.6 in all plotting functions, but is still compatible with older versions of \pkg{igraph}.}

\item{\code{loadNetwork()} supports comment lines in the network files.}

}
}
}

\section{Changes in version 1.50}{

\subsection{New functions}{

\itemize{

\item{\code{generateTimeSeries()} generates random state sequences from an existing network.}

\item{\code{plotSequence()} and \code{sequenceToLaTeX()} plot and export sequences of states similar to \code{plotAttractors()} and \code{attractorsToLaTeX()}.}

\item{\code{getAttractorSequence()} extracts the states of a single synchronous attractor from an attractor information structure as a data frame.}

\item{\code{generateState()} provides a simple way to specify network states using partial assignments.}
}
}

\subsection{Modified functions}{

\itemize{

\item{\code{getPathToAttractor()} has an additional parameter \code{includeAttractorStates} specifying which attractor states should be included in the path. The default behaviour has been changed to include all attractor states.}

\item{\code{generateRandomNKNetwork()} now supports the generation of networks using specific classes of functions. For this purpose, it has two new parameters \code{validationFunction} and \code{failureIterations}.}

\item{By default, \code{loadNetwork()} no longer changes gene names to lower case. If this behaviour is desired, it can be reactivated using the new \code{lowercaseGenes} parameter.}

\item{\code{stateTransition()} now names the state vector using the gene names.}

\item{\code{plotAttractors()} has a new parameter \code{drawLegend} to disable the legend.}

\item{The \code{randomChainLength} parameter of \code{getAttractors()} now defaults to 10000.}

\item{\code{getAttractors()}, \code{reconstructNetwork()} and \code{markovSimulation()} can now be interrupted by the user.}
}
}
}
