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
                   readableFunctions = FALSE)
}
\arguments{
  \item{measurements}{
	This can either be an object of class \code{TransitionTable} as returned by \code{\link{getTransitionTable}}, or a set of time series of measurements. In this case, \code{measurements} must be a list of matrices, each corresponding to one time series. Each row of these matrices contains measurements for one gene on a time line, i. e. column \code{i+1} contains the successor states of column \code{i}. The genes must be the same for all matrices in the list. Real-valued time series can be binarized using \code{\link{binarizeTimeSeries}}.
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
}
\details{
	Both algorithms iterate over all possible input combinations. While Best-Fit Extension is capable of returning functions that do not perfectly explain the measurements (for example, if there are inconsistent measurements or if \code{maxK} was specified too small), REVEAL only finds functions that explain all measurements. For more information, please refer to the cited publications.
}
\value{
Returns an object of class \code{BooleanNetworkCollection}, which is similar to \code{BooleanNetwork} except that it may contain more than one function per gene. A specific Boolean network can be extracted using \code{\link{chooseNetwork}}.
The return value consists of the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the alternative transition functions for the i-th gene. Each element is a list of transition functions with the following sub-components:
  \describe{
  	\item{input}{A vector of indices of the genes that serve as the input of the Boolean transition function. If the function has no input (i.e. the gene is constant), the vector consists of a zero element.}
  	\item{func}{The transition function in truth table representation. This vector has \code{2^length(input)} entries, one for each combination of input variables. If the gene is constant, the function is -1.}
  	\item{expression}{A string representation of the underlying Boolean expression}
  }}

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
data(cellcycle)

# extract transition table
tt <- getTransitionTable(getAttractors(net))

# reconstruct networks from transition table
net2 <- reconstructNetwork(tt,method="bestfit",maxK=6)
net3 <- reconstructNetwork(tt,method="reveal",maxK=6)

# print original net and reconstructed nets
print(net$interactions)
print(net2$interactions)
print(net3$interactions)

}

\keyword{Boolean network
	reconstruction
	REVEAL
	best-fit extension}
