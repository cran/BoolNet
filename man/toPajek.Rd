\name{toPajek}
\Rdversion{1.1}
\alias{toPajek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export a network to the Pajek file format
}
\description{
Exports a network to the Pajek file format to visualize transition trajectories. For more information on Pajek, please refer to \url{http://pajek.imfm.si/doku.php}
}
\usage{
toPajek(attractorInfo, file = "boolean.net", includeLabels=FALSE)
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}. As the transition table information in this structure is required, \code{getAttractors} must be called in synchronous mode and with \code{returnTable} set to TRUE.
}
  \item{file}{
	The name of the output file for Pajek. Defaults to "boolean.net".
}
 \item{includeLabels}{
	If set to true, the vertices of the graph in the output file are labeled with the binary encodings of the states. Defaults to FALSE.
}
}
\value{
This function has no return value.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{toSBML}}
}
\examples{
# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(cellcycle)

# export to Pajek
toPajek(attractors, file="pajek_export.net")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pajek
	Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph}
