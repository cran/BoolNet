\name{igf}
\alias{igf}
\docType{data}
\title{
Boolean model of the IGF pathway
}
\description{
A small Boolean model of major components of the IGF (Insuline-like growth receptor) pathway. Through IRS, IGF activates the well-known PI3K-Akt-mTOR signalling cascade. This cascade is finally inactivated by a feedback inhibion of IRS.

The model simplifies several complex formations and cascades by representing them as single nodes and specifying time delays instead. It therefore demonstrates the usage of temporal Boolean networks in \pkg{BoolNet}.
}
\usage{data(igf)}
\format{
  This data set consists of a variable \code{igf} of class \code{SymbolicBooleanNetwork} with 5 genes. The class \code{SymbolicBooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\examples{
data(igf)

sim <- simulateSymbolicModel(igf)
plotAttractors(sim)
}