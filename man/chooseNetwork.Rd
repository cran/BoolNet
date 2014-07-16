\name{chooseNetwork}
\alias{chooseNetwork}
\title{
Extract a single Boolean network from a probabilistic Boolean network
}
\description{
Creates a \code{BooleanNetwork} object with exactly one function per gene by extracting a specified set of transition functions from a \code{ProbabilisticBooleanNetwork} or \code{BooleanNetworkCollection} object.
}
\usage{
chooseNetwork(probabilisticNetwork, 
              functionIndices, 
              dontCareValues=NULL, 
              readableFunctions=FALSE)
}
\arguments{
  \item{probabilisticNetwork}{
A \code{ProbabilisticBooleanNetwork} or \code{BooleanNetworkCollection} object as returned by \code{\link{reconstructNetwork}} or \code{\link{loadNetwork}}
}
  \item{functionIndices}{
A vector of function indices with one entry for each gene
}
  
  \item{dontCareValues}{
If \code{probabilisticNetwork} is of class \code{BooleanNetworkCollection}, this specifies the values to fill in for "don't care" (*) values in the truth tables of the transition functions. This is a list containing one vector of Boolean values for each gene. The lengths of the vectors must coincide with the numbers of "don't care" values in the functions.
}
  \item{readableFunctions}{
If \code{probabilisticNetwork} is of class \code{BooleanNetworkCollection}, the string representations of the transition functions must be refreshed after filling in values for the "don't care" entries. 
This parameter specifies if readable DNF representations of the transition function truth tables are generated and displayed when the network is printed. If set to FALSE, the truth table result column is displayed. If set to "canonical", a canonical Disjunctive Normal Form is generated from each truth table. If set to "short", the canonical DNF is minimized by joining terms (which can be time-consuming for functions with many inputs). If set to TRUE, a short DNF is generated for functions with up to 12 inputs, and a canonical DNF is generated for functions with more than 12 inputs.
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
