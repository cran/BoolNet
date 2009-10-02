\name{chooseNetwork}
\alias{chooseNetwork}
\title{
Build a network from a collection of transition functions
}
\description{
Extracts the specified set of transition functions from a \code{BooleanNetworkCollection} object. This is necessary due to the fact that the reconstruction algorithms may return several transition functions per gene if no unique network could be reconstructed.
}
\usage{
chooseNetwork(networkCollection, functionIndices)
}
\arguments{
  \item{networkCollection}{
A \code{BooleanNetworkCollection} object as returned by \code{\link{reconstructNetwork}}
}
  \item{functionIndices}{
A vector of function indices with one entry for each gene
}
}
\value{
Returns an object of class \code{BooleanNetwork} consisting of the transition functions whose indices were specified in \code{functionIndices}. The class \code{BooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.

Constant genes are automatically fixed (e.g. knocked-out or over-expressed). This means that they are always set to the constant value, and states with the complementary value are not considered in transition tables etc. If you would like to change this behaviour, use \code{\link{fixGenes}} to reset the fixing.
}
\seealso{
\code{\link{reconstructNetwork}}, \code{\link{loadNetwork}}
}
\examples{
library(BoolNet)

# load example data
data(cellcycle)

# extract transition table
tt <- getTransitionTable(getAttractors(net))

# reconstruct networks from transition table
netCollection <- reconstructNetwork(tt,method="bestfit",maxK=6)

# extract a unique network
net2 <- chooseNetwork(netCollection,rep(1,length(netCollection$genes)))

# get attractors from this network
print(getAttractors(net2))

}
\keyword{Boolean network
	reconstruction}
