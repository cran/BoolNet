\name{NEWS}
\title{NEWS file for the \pkg{BoolNet} package}

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
