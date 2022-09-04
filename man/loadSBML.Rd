\name{loadSBML}
\alias{loadSBML}
\title{
Load an SBML document
}
\description{
Loads an SBML document that specifies a qualitative model using the \code{sbml-qual} extension package.
}
\usage{
loadSBML(file, symbolic=FALSE)
}

\arguments{
  \item{file}{
The SBML document to be imported
}

  \item{symbolic}{
If set to \code{TRUE}, the function returns an object of class \code{SymbolicBooleanNetwork} with an expression tree representation. Otherwise, it returns an object of class \code{BooleanNetwork} with a truth table representation.
}
}
\details{
The import assumes an SBML level 3 version 1 document with the \code{sbml-qual} extension package version 1.0.
\pkg{BoolNet} only supports a subset of the \code{sbml-qual} standard. The function tries to import those documents that describe a logical model with two possible values per species. It does not support general logical models with more than two values per species or Petri nets. 

Further details on the import:
\itemize{
\item{The import supports multiple function terms with the same output for a transition and interprets them as a disjunction, as proposed in the specification.}
\item{Comparison operators are converted to the corresponding Boolean expressions.}
\item{Compartments are ignored.}
}

For the import, the \pkg{XML} package is required.
}
\value{
Returns a structure of class \code{BooleanNetwork} or \code{SymbolicBooleanNetwork}, as described in \code{\link{loadNetwork}}.
}
\references{
\url{http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Qualitative_Models_(qual)}
}

\seealso{
\code{\link{toSBML}}, \code{\link{loadNetwork}}
}
\examples{
\dontrun{
# load the cell cycle network
data(cellcycle)

fil <- tempfile()
# export the network to SBML
toSBML(cellcycle, fil)

# reimport the model
print(loadSBML(fil))
}
}