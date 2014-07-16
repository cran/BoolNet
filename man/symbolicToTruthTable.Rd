\name{symbolicToTruthTable}
\alias{symbolicToTruthTable}

\title{
Convert a symbolic network into a truth table representation
}
\description{
Converts an object of class \code{SymbolicBooleanNetwork} into an object of class \code{BooleanNetwork} by generating truth tables from the symbolic expression trees.
}
\usage{
symbolicToTruthTable(network)
}

\arguments{
  \item{network}{
An object of class \code{SymbolicBooleanNetwork} to be converted.
}
}

\details{
The symbolic network \code{network} must not contain temporal operators, as these are not compatible with the truth table representation in \code{BooleanNetwork} objects.
}

\value{
Returns an object of class \code{BooleanNetwork}, as described in \code{\link{loadNetwork}}.
}

\seealso{
\code{\link{truthTableToSymbolic}}, \code{\link{loadNetwork}}
}
\examples{
# Convert a truth table representation into a 
# symbolic representation and back
data(cellcycle)

symbolicNet <- truthTableToSymbolic(cellcycle)
print(symbolicNet)

ttNet <- symbolicToTruthTable(symbolicNet)
print(cellcycle)
}
\keyword{Boolean network
         symbolic Boolean network
         conversion}
