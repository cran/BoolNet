\name{print.ProbabilisticBooleanNetwork}
\Rdversion{1.1}
\alias{print.ProbabilisticBooleanNetwork}
\alias{print.BooleanNetworkCollection}
\title{
Print a probabilistic Boolean network
}
\description{
A specialized method to print an object of class \code{ProbabilisticBooleanNetwork}. For backward compatibility, this method also prints objects of class \code{BooleanNetworkCollection}, which have been replaced by \code{ProbabilisticBooleanNetwork}.
This prints all alternative transition functions and their probabilities. If the network is the result of a reconstruction from time series measurements, it also outputs the error the functions make on the time series. If genes are knocked-out or over-expressed, these genes are listed below the functions.
}
\usage{
print.ProbabilisticBooleanNetwork(x, ...)
                               
print.BooleanNetworkCollection(x, ...)
}
\arguments{

\item{x}{
An object of class \code{ProbabilisticBooleanNetwork} or \code{BooleanNetworkCollection} to be printed
}

\item{\dots}{
Further parameters for the \code{\link{print}} method. Currently not used.
}
}
\value{
Invisibly returns the printed object
}

\seealso{
\code{\link{print}}, \code{\link{reconstructNetwork}}, \code{\link{loadNetwork}}
}
\keyword{network
	Boolean network
	probabilistic Boolean network
	PBN
	print}
