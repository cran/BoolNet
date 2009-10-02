\name{binarizeTimeSeries}
\alias{binarizeTimeSeries}
\title{
Binarize a set of real-valued time series
}
\description{
Binarizes a set of real-valued time series using k-means clustering.
}
\usage{
binarizeTimeSeries(measurements, nstart = 100, iter.max = 1000, 
                   dropInsignificant = TRUE, alpha = 0.05)
}
\arguments{
  \item{measurements}{
A list of matrices, each corresponding to one time series. Each row of these matrices contains real-valued measurements for one gene on a time line, i. e. column \code{i+1} contains the successor states of column \code{i+1}. The genes must be the same for all matrices in the list. 
}
  \item{nstart}{
The number of restarts for k-means. See \code{\link{kmeans}} for details.
}
  \item{iter.max}{
The maximum number of iterations for k-means. See \code{\link{kmeans}} for details.
}
  \item{dropInsignificant}{
If this is set to true, genes whose binarizations are insignificant in the Wilcoxon tests (see Details) are removed from the binarized time series. Otherwise, a warning is printed if such genes exist.
}
  \item{alpha}{
The significance level for the tests (see Details).
}
}
\details{
For each gene, k-means clusterings are performed to determine a good separation of groups. The values belonging to the cluster with the smaller centroid are set to 0, and the values belonging to the greater centroid are set to 1. The difference between the two clusters is tested for each gene using a two-sample Wilcoxon test. If \code{dropInsignificant} is true, all genes with a p-value greater than or equal to \code{alpha} are removed.
}
\value{
Returns a list with the following elements:
\item{binarizedMeasurements}{A list of matrices with the same structure as \code{measurements} containing the binarized time series measurements}
\item{pvals}{The p-values corresponding to the remaining genes}
\item{thresholds}{The thresholds used for binarization}
}
\seealso{
\code{\link{reconstructNetwork}}
}
\examples{
library(BoolNet)

# generate random measurements
measurements <- lapply(1:5,function(i)
			{
				m <- matrix(runif(n=150),ncol=15,nrow=10)
				rownames(m) <- paste("Gene",1:10)
				return(m)
			})
			
# perform binarization
bin <- binarizeTimeSeries(measurements)
print(bin)

# reconstruct a network from the data
reconstructed <- reconstructNetwork(bin$binarizedMeasurements,method="bestfit")
print(reconstructed)
}
\keyword{Boolean network
	binarization
	binarize
	k-means}
