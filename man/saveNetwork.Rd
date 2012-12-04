\name{saveNetwork}
\alias{saveNetwork}
\title{
Save a network
}
\description{
Saves synchronous, asynchronous or probabilistic networks in the \pkg{BoolNet} network file format .
}
\usage{
saveNetwork(network, 
            file, 
            generateDNFs = FALSE,
            saveFixed = TRUE)
}
\arguments{
  \item{network}{
An object of class \code{BooleanNetwork} to be exported
}
  \item{file}{
The name of the network file to be created
}

  \item{generateDNFs}{
This parameter specifies whether formulae in Disjunctive Normal Form are exported instead of the expressions that describe the transition functions. If set to FALSE, the original expressions are exported. If set to "canonical", a canonical Disjunctive Normal Form is generated from each truth table. If set to "short", the canonical DNF is minimized by joining terms (which can be time-consuming for functions with many inputs). If set to TRUE, a short DNF is generated for functions with up to 12 inputs, and a canonical DNF is generated for functions with more than 12 inputs.
}

  \item{saveFixed}{
If set to TRUE, knock-outs and overexpression of genes override their transition functions. That is, if a gene in the network is fixed to 0 or 1, this value is saved, regardless of the transition function. If set to FALSE, the transition function is saved. Defaults to TRUE.
}

}
\details{
The network is saved in the \pkg{BoolNet} file format (see \code{\link{loadNetwork}} for details).

If the expressions in the transition functions cannot be parsed or \code{generateDNFs} is true, a DNF representation of the transition functions is generated.
}

\seealso{
\code{\link{loadNetwork}}
}
\examples{
# load the cell cycle network
data(cellcycle)

# save it to a file
saveNetwork(cellcycle, file="cellcycle.txt")

# reload the model
print(loadNetwork("cellcycle.txt"))
}
\keyword{Boolean network
  probabilistic Boolean network
  PBN
	file
	export}

