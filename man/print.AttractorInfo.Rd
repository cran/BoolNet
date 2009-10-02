\name{print.AttractorInfo}
\Rdversion{1.1}
\alias{print.AttractorInfo}
\title{
Print attractor cycles
}
\description{
Specialized print method to print the attractor cycles stored in an \code{AttractorInfo} object. The states of the attractors are printed in binary encoding in the order they are reached.
}
\usage{
print.AttractorInfo(x, quote = FALSE, max.levels = NULL, width = getOption("width"), ...)
}
\arguments{
  \item{x}{
An object of class \code{AttractorInfo} to be printed
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
\code{\link{print}}, \code{\link{getAttractors}}
}
\keyword{attractor
	print}