\name{testNetworkProperties}
\alias{testNetworkProperties}
\alias{testIndegree}
\alias{testAttractorRobustness}
\title{
Test properties of networks by comparing them to random networks
}
\description{
This is a general function designed to determine unique properties of biological networks by comparing them to a set of randomly generated networks with similar structure.
}
\usage{
testNetworkProperties(network, 
                      numRandomNets = 100, 
                      testFunction = "testIndegree", 
                      testFunctionParams = list(), 
                      accumulation = c("characteristic", "kullback_leibler"), 
                      sign.level = 0.05, 
                      drawSignificanceLevel = TRUE, 
                      klBins, 
                      klMinVal = 1e-05, 
                      linkage = c("uniform", "lattice"), 
                      functionGeneration = c("uniform", "biased"), 
                      simplify = FALSE, 
                      noIrrelevantGenes = TRUE, 
                      d_lattice = 1, 
                      zeroBias = 0.5, 
                      title = "", 
                      xlab, 
                      xlim, 
                      breaks = 30, 
                      ...)
}
\arguments{
  \item{network}{
  A network structure of class \code{BooleanNetwork}
}
  \item{numRandomNets}{
   The number of random networks to generate for comparison
}
  \item{testFunction}{
  The name of a function that calculates characteristic values that describe properties of the network. There are two built-in functions: "testIndegree" calculates the in-degrees of states in the network, and "testAttractorRobustness" counts the occurrences of attractors in perturbed copies. It is possible to supply user-defined functions here. See Details.
}
  \item{testFunctionParams}{
A list of parameters to \code{testFunction}. The elements of the list depend on the chosen function.
}
  \item{accumulation}{
If "characteristic" is chosen, the test function is required to return a single value that describes the network. In this case, a histogram of these values in random networks is plotted, and the value of the original network is inserted as a vertical line.

If "kullback_leibler" is chosen, the test function can return a vector of values which is regarded as a sample from a probability distribution. In this case, the Kullback-Leibler distances of the distributions from the original network and each of the random networks are calculated and plotted in a histogram. The Kullback-Leibler distance measures the difference between two probability distributions. In this case, the resulting histogram shows the distribution of differences between the original network and randomly generated networks.
}
  \item{sign.level}{
if \code{accumulation="characteristic"}, this specifies a significance level for a computer-intensive test. The test is significant if the characteristic value is greater than at least \code{(1-sign.level)*100}\% of the characteristic values of the random networks
}
  \item{drawSignificanceLevel}{
If \code{accumulation="characteristic"} and this is true, a vertical line is plotted for the significance level in the histogram.
}
  \item{linkage, functionGeneration, simplify, noIrrelevantGenes, d_lattice, zeroBias}{
  The corresponding parameters of \code{\link{generateRandomNKNetwork}} used to generate the random networks. This allows for customization of the network generation process. The three remaining parameters of \code{\link{generateRandomNKNetwork}} are set to values that ensure structural similarity to the original network: The parameters \code{n} and \code{k} are set to the corresponding values of the original network, and \code{topology="fixed"}. 
}

  \item{klBins}{
If \code{accumulation="kullback_leibler"}, the number of bins used to discretize the samples for the Kullback-Leibler distance calculations. By default, each unique value in the samples has its own bin, i.e. no further discretization is performed. The influence of discretization on the resulting histogram may be high.
}
  \item{klMinVal}{
If \code{accumulation="kullback_leibler"}, this defines the minimum probability for the calculation of the Kullback-Leibler distance to ensure stability of the results.  
}

  \item{title}{
  The title of the plots. This is empty by default.
}
  \item{xlab}{
Customizes label of the x axis of the histogram. For the built-in test functions, the x axis label is set automatically.
}
  \item{xlim}{
Customizes the limits of the x axis of the histogram. For the built-in test functions, suitable values are chosen automatically. 
}
  \item{breaks}{
Customizes the number of breaks in the 
}
  \item{\dots}{
Further graphical parameters for \code{\link{hist}}
}
}
\details{
This function generically compares properties of biological networks to a set of random networks. It can be extended by supplying custom functions to the parameter \code{testFunction}. Such a function must have the signature

\code{function(network,
               accumulate=TRUE,
               params)}

\describe{
\item{network}{This is the network to test. In the process of the comparison, both the original network and the random networks are passed to the function}

\item{accumulate}{If \code{accumulate=TRUE}, the function must return a single value quantifying the examined property of the network. If \code{accumulate=FALSE}, the function can return a vector of values (e.g., one value for each gene/state etc.)}

\item{params}{A list of further parameters for the function supplied by the user in \code{testFunctionParams} (see above). This can contain any type of information your test function needs.}
}

Two built-in functions for synchronous Boolean networks already exist:
\describe{
\item{testIndegree}{This function is based on the observation that, often, in biological networks, many state transitions lead to the same states. In other words, there is a small number of "hub" states. In the state graph, this means that the in-degree of some states (i.e., the number of transitions leading to it) is high, while the in-degree of many other states is 0. We observed that random networks do not show this behaviour, thus it may be a distinct property of biological networks.

The function does not require any parameter entries in \code{params}. If \code{accumulate=FALSE}, it returns the in-degrees of all synchronous states in the network. If \code{accumulate=TRUE}, the Gini index of the in-degrees is returned as a characteristic value of the network. The Gini index is a measure of inequality. If all states have an in-degree of 1, the Gini index is 0. If all state transitions lead to one single state, the Gini index is 1.

This function requires the \pkg{igraph} package for the analysis of the in-degrees.
}

\item{testAttractorRobustness}{This function tests the robustness of attractors in a network to noise. We expect attractors in a real network to be less susceptible to noise than attractors in randomly generated networks, as biological processes can be assumed to be comparatively stable. The function generates a number of perturbed copies of the network and checks whether the original attractors can still be found in the network.

\code{params} can hold a number of parameters:
\describe{
\item{copies}{The number of perturbed copies to generate}
\item{perturb, method, simplify, readableFunctions, excludeFixed, maxNumBits, numStates}{The corresponding parameters of \code{\link{perturbNetwork}} that influence the way the network is perturbed.}
}
If \code{accumulate=FALSE}, the function returns a vector of percentages of original attractors found in each of the perturbed copies of the original network. If \code{accumulate=TRUE}, the function returns the overall percentage of original attractors found in all perturbed copies.
}
} 
               
}
\value{
The function returns a list with the following elements
 \item{hist}{The histogram that was plotted. The type of histogram depends on the parameter \code{accumulation}.}
 \item{pval}{If \code{accumulation="characteristic"}, a p-value for the alternative hypothesis that the test statistic value of the supplied network is greater than the value of a randomly generated network is supplied.}
 \item{significant}{If \code{accumulation="characteristic"}, this is true for \code{pval < sign.level}.} 
}

\seealso{
\code{\link{generateRandomNKNetwork}}, \code{\link{perturbNetwork}}, \code{\link{plotStateGraph}}, \code{\link{getAttractors}}
}
\examples{
library(BoolNet)

# load mammalian cell cycle network
data(cellcycle)

# compare the in-degrees of the states in the
# cell cycle network to random networks
testNetworkProperties(cellcycle, testFunction="testIndegree")

# compare the in-degrees of the states in the
# cell cycle network to random networks,
# and plot the Kullback-Leibler distances of the 100 experiments
testNetworkProperties(cellcycle, testFunction="testIndegree",
                      accumulation = "kullback_leibler")

# compare the robustness of attractors in the cell cycle network
# to random networks
testNetworkProperties(cellcycle, testFunction="testAttractorRobustness",
                      testFunctionParams=list(copies=10))
}

\keyword{robustness analysis
         perturbation
         in-degree}
