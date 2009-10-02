\name{print.BooleanNetworkCollection}
\Rdversion{1.1}
\alias{print.BooleanNetworkCollection}
\title{
Print a collection of Boolean functions from a reconstruction algorithm
}
\description{
A specialized method to print an object of class \code{BooleanNetworkCollection}. This prints all alternative functions returned by a reconstruction algorithm for each of the genes. It also outputs the error the functions make on the input time series. If genes are knocked-out or over-expressed, these genes are listed below the functions.
}
\usage{
print.BooleanNetworkCollection(x, quote = FALSE, max.levels = NULL, width = getOption("width"), ...)
}
\arguments{
  \item{x}{
An object of class \code{BooleanNetworkCollection} to be printed
}
  \item{quote}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{max.levels}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{width}{
Standard parameter from the generic \code{\link{print}} method, not used here.
}
  \item{\dots}{
Further parameters for the \code{\link{print}} method.
}
}
\value{
Invisibly returns the printed object
}

\seealso{
\code{\link{print}}, \code{\link{reconstructNetwork}}
}
\keyword{network
	Boolean network
	print}