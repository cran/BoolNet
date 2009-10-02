\name{BoolNet-package}
\Rdversion{1.1}
\alias{BoolNet-package}
\alias{BoolNet}
\docType{package}
\title{Analyze Boolean networks}
\description{Tools for reconstruction, analysis and visualization of Boolean networks, in particular for the identification of attractors in gene-regulatory networks}
\details{
\tabular{ll}{
Package: \tab BoolNet\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2009-10-02\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
This package provides useful methods for the reconstruction and generation of Boolean networks and for their analysis. In particular, it is designed for the analysis of gene-regulatory networks.

The \code{\link{reconstructNetwork}} function infers Boolean networks from time series measurements using several popular reconstruction algorithms.
Boolean networks can also be expressed in a description language and loaded from files using \code{\link{loadNetwork}}.
Via \code{\link{generateRandomNKNetwork}} and \code{\link{perturbNetwork}}, the package supports various methods of generating random networks and perturbing existing networks for robustness analysis.

The \code{\link{getAttractors}} function identifies attractor cycles in a Boolean network. Attractors can be identified by exhaustive search or heuristic methods.

The package also provides methods to visualize state transitions and basins of attraction (\code{\link{plotStateGraph}}),
to plot the wiring of a network (\code{\link{plotNetworkWiring}}), to plot the attractor state table (\code{\link{plotAttractors}}), and to export them to LaTeX (\code{\link{attractorsToLaTeX}}) and Pajek (\code{\link{toPajek}}).

Transition tables of the network can be analyzed using \code{\link{getTransitionTable}}.
}
\author{
Christoph Muessel, Dao Zhou, Hans Kestler

Maintainer: Hans Kestler <hans.kestler@uni-ulm.de>

}
\references{
S. A. Kauffman (1969), Metabolic stability and epigenesis in randomly constructed nets. J. Theor. Biol. 22:437--467.

S. A. Kauffman (1993), The Origins of Order. Oxford University Press.

Further references are listed in the corresponding help sections.
}
\keyword{package 
	Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph
	reconstruction}

\examples{
# Perform a robustness test on a network
# by counting the numbers of perturbed networks
# containing the attractors of the original net
# and measuring their basin sizes

library(BoolNet)

# load mammalian cell cycle network
data(cellcycle)

# get attractors in original network
attrs <- getAttractors(net, canonical=TRUE)

# create 100 perturbed copies of the network,
# search for attractors and measure basin sizes
perturbationResults <- sapply(1:100,function(i)
{
  # perturb network and identify attractors
	perturbedNet <- perturbNetwork(net, perturb="functions", method="bitflip", maxNumBits=1)
	perturbedAttrs <- getAttractors(perturbedNet, canonical=TRUE)
	
	# check whether the attractors in the original network exist in the perturbed network
	attractorIndices <- sapply(attrs$attractors,function(attractor1)
				{
					index <- which(sapply(perturbedAttrs$attractors,function(attractor2)
						{
							identical(attractor1,attractor2)
						}))
					if (length(index) == 0)
						NA
					else
						index
				})

  # if the attractors exist, measure the sizes of their basins
	basinSizes <- sapply(attractorIndices,function(attractorIndex)
				{
					if (is.na(attractorIndex))
						NA
					else
						perturbedAttrs$attractors[[attractorIndex]]$basinSize
				})
	return(c(attractorIndices,basinSizes))
})

# perturbationResults now contains a matrix
# with the first 2 columns specifying the indices or the 
# original attractors in the perturbed network 
# (or NA if the attractor was not found) and the next 2 
# columns counting the numbers of states
# in the basin of attraction (or NA if the attractor was not found)

# measure the total numbers of occurrences in the perturbed copies
# for each of the original attractors
numOccurrences <- apply(perturbationResults[1:length(attrs$attractors),,drop=FALSE], 1,
                        function(row)sum(!is.na(row)))


# print information
cat("Number of occurrences of the original attractors",
	"in 100 perturbed copies of the network:\n")
for (i in 1:length(attrs$attractors))
{
  cat("Attractor ",i,": ",numOccurrences[i],"\n",sep="")
}

# plot a histogram of the basin sizes
# - in this case mostly unspectacular, as the attractor
# usually either disappears or does not change.

par(mfrow=c(length(attrs$attractors),1))
for (i in 1:length(attrs$attractors))
{
  basinSizes <- perturbationResults[i + length(attrs$attractors),]
  basinSizes <- basinSizes[!is.na(basinSizes)]
  hist(basinSizes,
       breaks=seq(0,1024,length.out=20),
       xlab="Size of basin of attraction",
       main=paste("Basin sizes of attractor",i,"in 100 perturbed networks"))
  
  # for comparison, draw basin size of original attractor
  abline(v=attrs$attractors[[i]]$basinSize,col="red")
}
}
