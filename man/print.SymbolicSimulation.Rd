\name{print.SymbolicSimulation}
\alias{print.SymbolicSimulation}
\title{
Print simulation results
}
\description{
Specialized print method to print the information stored in an \code{AttractorInfo} object. By default, the states of the identified attractors are printed in a binary encoding. Furthermore, the state transition graph and the sequences from the start states to the attractors can be printed. The method can print either the full states, or only the active genes of the states.
}
\usage{
\method{print}{SymbolicSimulation}(x, 
                                   activeOnly = FALSE, 
                                   sequences = FALSE, 
                                   graph = FALSE, 
                                   attractors = TRUE, 
                                   ...)
}

\arguments{

 \item{x}{
An object of class \code{SymbolicSimulation} to be printed.
}
  
 \item{activeOnly}{
If set to true, a state is represented by a list of active genes (i.e., genes which are set to 1). If set to false, a state is represented by a binary vector with one entry for each gene, specifying whether the gene is active or not. Defaults to \code{FALSE}. 
}
  \item{sequences}{
If set to true and if \code{\link{simulateSymbolicModel}} has been started with \code{returnSequences=TRUE}, the sequences from the start states to the attractors are printed. Defaults to \code{FALSE}.
}
  \item{graph}{
If set to true if \code{\link{simulateSymbolicModel}} has been started with \code{returnGraph=TRUE}, the state transition table is printed. Defaults to \code{FALSE}.
}
  \item{attractors}{
If set to true if \code{\link{simulateSymbolicModel}} has been started with \code{returnAttractor=TRUE}, the state transition table is printed. Defaults to \code{TRUE}.
}
  \item{\dots}{
Further parameters for the \code{\link{print}} method. Currently not used.
}
}

\value{
Invisibly returns the printed object
}
\seealso{
\code{\link{simulateSymbolicModel}}
}
\keyword{print
         symbolic Boolean network}
