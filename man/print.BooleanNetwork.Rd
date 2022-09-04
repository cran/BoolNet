\name{print.BooleanNetwork}
\Rdversion{1.1}
\alias{print.BooleanNetwork}
\title{
Print a Boolean network
}
\description{
A specialized method to print an object of class \code{BooleanNetwork}. This prints the transition functions of all genes. If genes are knocked-out or over-expressed, these genes are listed below the functions.
}
\usage{
\method{print}{BooleanNetwork}(x, ...)
}
\arguments{

\item{x}{
An object of class \code{BooleanNetwork} to be printed
}

\item{\dots}{
Further parameters for the \code{\link{print}} method. Currently not used.
}
}
\value{
Invisibly returns the printed object
}

\seealso{
\code{\link{print}}, \code{\link{loadNetwork}}
}