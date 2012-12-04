\name{toSBML}
\alias{toSBML}
\title{
Export a network to SBML
}
\description{
Exports a synchronous or asynchronous Boolean network to SBML with the \code{sbml-qual} extension package.
}
\usage{
toSBML(network, 
       file,
       generateDNFs = FALSE,
       saveFixed = TRUE)
}
\arguments{
  \item{network}{
An object of class \code{BooleanNetwork} to be exported
}
  \item{file}{
The name of the SBML file to be created
}

  \item{generateDNFs}{
This parameter specifies whether formulae in Disjunctive Normal Form are exported instead of the expressions that describe the transition functions. If set to FALSE, the original expressions are exported. If set to "canonical", a canonical Disjunctive Normal Form is generated from each truth table. If set to "short", the canonical DNF is minimized by joining terms (which can be time-consuming for functions with many inputs). If set to TRUE, a short DNF is generated for functions with up to 12 inputs, and a canonical DNF is generated for functions with more than 12 inputs.
}

  \item{saveFixed}{
If set to TRUE, knock-outs and overexpression of genes override their transition functions. That is, if a gene in the network is fixed to 0 or 1, this value is exported, regardless of the transition function. If set to FALSE, the transition function is exported. Defaults to TRUE.
}
}
\details{
The export creates an SBML file describing a general logical model that corresponds to the Boolean network. Importing tools must support the \code{sbml-qual} extension package version 1.0.

The export translates the expressions that describe the network transition functions to a MathML description. If these expressions cannot be parsed or \code{generateDNFs} is true, a DNF representation of the transition functions is generated and exported.
}

\references{
\url{http://sbml.org/Community/Wiki/SBML_Level_3_Proposals/Qualitative_Models}
}

\seealso{
\code{\link{loadSBML}}, \code{\link{loadNetwork}}, \code{\link{saveNetwork}}, \code{\link{toPajek}}
}
\examples{
# load the cell cycle network
data(cellcycle)

# export the network to SBML
toSBML(cellcycle, file="cellcycle.sbml")

# reimport the model
print(loadSBML("cellcycle.sbml"))
}
\keyword{SBML
        Boolean network
        file
        export}
