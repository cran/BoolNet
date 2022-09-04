\name{print.AttractorInfo}
\Rdversion{1.1}
\alias{print.AttractorInfo}
\title{
Print attractor cycles
}
\description{
Specialized print method to print the attractor cycles stored in an \code{AttractorInfo} object. For simple or steady-state attractors, the states of the attractors are printed in binary encoding in the order they are reached. For asynchronous complex/loose attractors, the possible transitions of the states in the attractor are printed. The method can print either the full states, or only the active genes of the states.
}
\usage{
\method{print}{AttractorInfo}(x,
      activeOnly = FALSE, 
      ...)
}
\arguments{

\item{x}{
An object of class \code{AttractorInfo} to be printed
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
\code{\link{print}}, \code{\link{getAttractors}}
}