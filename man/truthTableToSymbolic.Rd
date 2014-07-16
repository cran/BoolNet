\name{truthTableToSymbolic}
\alias{truthTableToSymbolic}
\title{
Convert a network in truth table representation into a symbolic representation
}
\description{
Converts an object of class \code{BooleanNetwork} into an object of class \code{SymbolicBooleanNetwork} by generating symbolic expression trees.
}
\usage{
truthTableToSymbolic(network, generateDNFs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{
An object of class \code{BooleanNetwork} to be converted.
}

\item{generateDNFs}{
This parameter specifies whether formulae in Disjunctive Normal Form are generated instead of the parsing the string expressions that describe the transition functions. If set to FALSE, the original expressions are parsed. If set to "canonical", a canonical Disjunctive Normal Form is generated from each truth table. If set to "short", the canonical DNF is minimized by joining terms (which can be time-consuming for functions with many inputs). If set to TRUE, a short DNF is generated for functions with up to 12 inputs, and a canonical DNF is generated for functions with more than 12 inputs.
}
}

\value{
Returns an object of class \code{SymbolicBooleanNetwork}, as described in \code{\link{loadNetwork}}.
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
