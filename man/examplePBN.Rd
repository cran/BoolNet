\name{examplePBN}
\alias{examplePBN}
\docType{data}
\title{
An artificial probabilistic Boolean network
}
\description{
An artificial probabilistic Boolean network example introduced by Shmulevich et al.
}
\usage{data(examplePBN)}

\details{
This artificial network is introduced by Shmulevich et al. for a step-by-step description of their Markov chain algorithm. It is included as a general example for a probabilistic Boolean network. The network consists of 3 genes, where gene 1 and gene 3 have two alternative transition functions, and gene 1 has a unique transition function.
}
\source{
I. Shmulevich, E. R. Dougherty, S. Kim, W. Zhang (2002), Probabilistic Boolean networks: a rule-based uncertainty model for gene regulatory networks. Bioinformatics 18(2):261--274.
}
\examples{
library(BoolNet)
data(examplePBN)

# the network is stored in a variable called 'examplePBN'
print(examplePBN)
}
\keyword{datasets
        probabilistic Boolean network
        PBN}
