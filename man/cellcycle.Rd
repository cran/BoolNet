\name{cellcycle}
\alias{cellcycle}
\alias{net}
\docType{data}
\title{
Mammalian cell cycle network
}
\description{
The mammalian cell cycle network as described by Faure et al.
}
\usage{data(cellcycle)}
\details{
This is an object of class \code{BooleanNetwork} with 10 genes describing the four phases of the mammalian cell cycle. The network has one attractor with one gene and one attractor with 7 genes. The class \code{BooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\source{
A. Faure, A. Naldi, C. Chaouiya and D. Thieffry (2006), Dynamical analysis of a generic Boolean model for the control of the mammalian cell cycle. Bioinformatics 22(14):e124--e131.
}
\examples{
data(cellcycle)

# the network is stored in a variable called 'net'
print(net)

}
\keyword{datasets
	cell cycle
	mammalian
	Boolean network}
