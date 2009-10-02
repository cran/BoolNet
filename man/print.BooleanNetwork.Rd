\name{print.BooleanNetwork}
\Rdversion{1.1}
\alias{print.BooleanNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a Boolean network
}
\description{
A specialized method to print an object of class \code{BooleanNetwork}. This prints the transition functions of all genes. If genes are knocked-out or over-expressed, these genes are listed below the functions.
}
\usage{
print.BooleanNetwork(x, quote = FALSE, max.levels = NULL, width = getOption("width"), ...)
}
\arguments{
  \item{x}{
An object of class \code{BooleanNetwork} to be printed
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
\code{\link{print}}, \code{\link{loadNetwork}}
}
\keyword{network
	Boolean network
	print}