\name{NEWS}
\title{NEWS file for the BoolNet package}

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
