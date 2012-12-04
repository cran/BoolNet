\name{loadSBML}
\alias{loadSBML}
\title{
Load an SBML document
}
\description{
Loads an SBML document that specifies a qualitative model using the \code{sbml-qual} extension package.
}
\usage{
loadSBML(file)
}

\arguments{
  \item{file}{
The SBML document to be imported
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
Returns a structure of class \code{BooleanNetwork}, as described in \code{\link{loadNetwork}}.
}
\references{
\url{http://sbml.org/Community/Wiki/SBML_Level_3_Proposals/Qualitative_Models}
}

\seealso{
\code{\link{toSBML}}, \code{\link{loadNetwork}}
}
\examples{
# load the cell cycle network
data(cellcycle)

# export the network to SBML
toSBML(cellcycle, "cellcycle.sbml")

# reimport the model
print(loadSBML("cellcycle.sbml"))
}
\keyword{SBML
        Boolean network
        file
        import}
