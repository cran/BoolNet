\name{generationFunctions}
\alias{generateCanalyzing}
\alias{generateNestedCanalyzing}
\title{
Generation functions for biologically relevant function classes
}
\description{
These generation functions randomly generate canalyzing or nested canalyzing Boolean functions. These functions are usually not called directly, but are supplied to the \code{functionGeneration} parameter of \code{\link{generateRandomNKNetwork}}.
}
\usage{
generateCanalyzing(input)
generateNestedCanalyzing(input)
}

\arguments{
  \item{input}{
  A vector of input gene indices for the Boolean function
}
}
\value{
A binary vector corresponding to the result column of the truth table that represents the canalyzing/nested canalyzing function.
}
\references{
S. Kauffman and C. Peterson and B. Samuelsson and C. Troein (2004), Genetic networks with canalyzing Boolean rules are always stable. PNAS 101(49):7102--17107.
}


\seealso{
\code{\link{generateRandomNKNetwork}}
}
\examples{
  # generate a random network with canalyzing functions
  net1 <- generateRandomNKNetwork(n=10, k=5, 
                                  functionGeneration="generateCanalyzing")
  print(net1)
                                  
  # generate a random network with nested canalyzing functions
  net2 <- generateRandomNKNetwork(n=10, k=5, 
                                  functionGeneration="generateNestedCanalyzing")
  print(net2)                                  
}
\keyword{canalyzing function
         nested canalyzing function,
         random network}

