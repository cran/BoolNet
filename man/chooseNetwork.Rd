\name{chooseNetwork}
\alias{chooseNetwork}
\title{
Extract a single Boolean network from a probabilistic Boolean network
}
\description{
Creates a \code{BooleanNetwork} object with exactly one function per gene by extracting a specified set of transition functions from a \code{ProbabilisticBooleanNetwork} object.
}
\usage{
chooseNetwork(probabilisticNetwork, functionIndices)
}
\arguments{
  \item{probabilisticNetwork}{
A \code{ProbabilisticBooleanNetwork} object as returned by \if{latex}{\cr}\code{\link{reconstructNetwork}} or \code{\link{loadNetwork}}
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
# load example data
data(examplePBN)

# extract a unique network
# - always use the first function
net <- chooseNetwork(examplePBN, rep(1, length(examplePBN$genes)))

# get attractors from this network
print(getAttractors(net))

}
\keyword{Boolean network
  probabilistic Boolean network
  PBN
  conversion
	reconstruction}
