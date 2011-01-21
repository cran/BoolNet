\name{reconstructNetwork}
\alias{reconstructNetwork}
\title{
Reconstruct a Boolean network from time series of measurements
}
\description{
Reconstructs a Boolean network from a set of time series or from a transition table using the best-fit extension algorithm or the REVEAL algorithm.
}
\usage{
reconstructNetwork(measurements, 
                   method = c("bestfit", "reveal"), 
                   maxK = 5, 
                   readableFunctions = FALSE,
                   allSolutions = FALSE)
}
\arguments{
  \item{measurements}{
	This can either be an object of class \code{TransitionTable} as returned by \if{latex}{\cr}\code{\link{getTransitionTable}}, or a set of time series of measurements. In this case, \code{measurements} must be a list of matrices, each corresponding to one time series. Each row of these matrices contains measurements for one gene on a time line, i. e. column \code{i+1} contains the successor states of column \code{i}. The genes must be the same for all matrices in the list. Real-valued time series can be binarized using \code{\link{binarizeTimeSeries}}.
}
  \item{method}{
	This specifies the reconstruction algorithm to be used. If set to "bestfit", Laehdesmaeki's Best-Fit Extension algorithm is employed. This algorithm is an improvement of the algorithm by Akutsu et al. with a lower runtime complexity. It determines the functions with a minimal error for each gene. If set to "reveal", Liang's REVEAL algorithm is used. This algorithm searches for relevant input genes using the mutual information between the input genes and the output gene. 
}
  \item{maxK}{
	The maximum number of input genes for one gene to be tested. Defaults to 5.
}
  \item{readableFunctions}{
	If this is true, readable DNF representations of the truth tables of the functions are generated. These DNF are displayed when the network is printed. The DNF representations are not minimized and can thus be very long. If set to FALSE, the truth table result column is displayed.
}

  \item{allSolutions}{
  If this is true, all solutions with the minimum error and up to \code{maxK} inputs are returned. By default, \code{allSolutions=FALSE}, which means that only the solutions with both minimum error and minimum \code{k} are returned.
  }
}
\details{
	Both algorithms iterate over all possible input combinations. While Best-Fit Extension is capable of returning functions that do not perfectly explain the measurements (for example, if there are inconsistent measurements or if \code{maxK} was specified too small), REVEAL only finds functions that explain all measurements. For more information, please refer to the cited publications.
}
\value{
Returns an object of class \code{ProbabilisticBooleanNetwork}, with each alternative function of a gene having the same probability. The structure is described in detail in \code{\link{loadNetwork}}. In addition to the standard components, each alternative transition function has a component \code{error} which stores the error of the function on the input time series data.
If a non-probabilistic network of class \code{BooleanNetwork} is required, a specific Boolean network can be extracted using \code{\link{chooseNetwork}}.

}
\references{
H. Laehdesmaeki, I. Shmulevich and O. Yli-Harja (2003), On Learning Gene-Regulatory Networks Under the Boolean Network Model. Machine Learning 52:147--167.

T. Akutsu, S. Miyano and S. Kuhara (2000). Inferring qualitative relations in genetic networks and metabolic pathways. Bioinformatics 16(8):727--734.

S. Liang, S. Fuhrman and R. Somogyi (1998), REVEAL, a general reverse engineering algorithm for inference of genetic network architectures. Pacific Symposium on Biocomputing 3:18--29.

}


\seealso{
\code{\link{binarizeTimeSeries}}, \code{\link{chooseNetwork}}
}
\examples{
library(BoolNet)

# load example data
data(yeastTimeSeries)

# perform binarization with k-means
bin <- binarizeTimeSeries(yeastTimeSeries)

# reconstruct networks from transition table
net <- reconstructNetwork(bin$binarizedMeasurements, method="bestfit", maxK=3)

# print reconstructed net
print(net)

# plot reconstructed net
\donttest{plotNetworkWiring(net)}
}

\keyword{probabilistic Boolean network
  PBN
	reconstruction
	REVEAL
	best-fit extension}
