\name{fixGenes}
\alias{fixGenes}

\title{
Simulate knocked-out or over-expressed genes
}
\description{
Simulates knocked-out or over-expressed genes by fixing the values of genes to 0 or 1, or turn off knock-out or over-expression of genes.
}
\usage{
fixGenes(network, fixIndices, values)
}
\arguments{
  \item{network}{
	The original network of class \code{BooleanNetwork} or \if{latex}{\cr}\code{ProbabilisticBooleanNetwork} containing the genes to be fixed
}
  \item{fixIndices}{
	A vector of names or indices of the genes to be fixed
}
  \item{values}{
	Either one single value, or a vector with the same length as \code{fixIndices}. For each gene, a value of 1 means that the gene is always turned on (over-expressed), a value of 0 means that the gene is always turned off (knocked-out), and a value of -1 means that the gene is not fixed.
}
}
\value{
Depending on the input, an object of class \code{BooleanNetwork} or \code{ProbabilisticBooleanNetwork} containing the fixed genes is returned. The classes \code{BooleanNetwork} and \if{latex}{\cr}\code{ProbabilisticBooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\seealso{
\code{\link{loadNetwork}}}

\examples{
# load example data
data(cellcycle)

# knock out gene CycD (index 1)
net <- fixGenes(cellcycle, 1, 0)
# or
net <- fixGenes(cellcycle, "CycD", 0)

# get attractors by exhaustive search
attractors <- getAttractors(net)

print(attractors)
}
\keyword{Boolean network
    probabilistic Boolean network
    PBN
		knock-out
		over-expression
		fixed gene
		fix}
