\name{cellcycle}
\alias{cellcycle}
\docType{data}
\title{
Mammalian cell cycle network
}
\description{
The mammalian cell cycle network as described by Faure et al.
}
\usage{data(cellcycle)}
\details{
The data consists of a variable \code{cellcycle} of class \code{BooleanNetwork} with 10 genes describing the four phases of the mammalian cell cycle. The network has one steady-state attractor. Furthermore, it has one synchronous attractor with 7 states and one asynchronous complex/loose attractor with 112 states. The class \code{BooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\source{
A. Faure, A. Naldi, C. Chaouiya and D. Thieffry (2006), Dynamical analysis of a generic Boolean model for the control of the mammalian cell cycle. Bioinformatics 22(14):e124--e131.
}
\examples{
data(cellcycle)

# the network is stored in a variable called 'cellcycle'
print(cellcycle)

}
\keyword{datasets
	cell cycle
	mammalian
	Boolean network}
