\name{binarizeTimeSeries}
\alias{binarizeTimeSeries}
\title{
Binarize a set of real-valued time series
}
\description{
Binarizes a set of real-valued time series using k-means clustering, edge detection, or scan statistics.
}
\usage{
binarizeTimeSeries(measurements, 
                   method = c("kmeans","edgeDetector","scanStatistic"), 
                   nstart = 100, 
                   iter.max = 1000, 
                   edge = c("firstEdge","maxEdge"),                    
                   scaling = 1, 
                   windowSize = 0.25, 
                   sign.level = 0.1,
                   dropInsignificant = FALSE)
}
\arguments{
  \item{measurements}{
A list of matrices, each corresponding to one time series. Each row of these matrices contains real-valued measurements for one gene on a time line, i. e. column \code{i+1} contains the successor states of column \code{i+1}. The genes must be the same for all matrices in the list. 
}
  \item{method}{The employed binarization technique. "kmeans" uses k-means clustering for binarization. "edgeDetector" searches for a large gradient in the sorted measurements. "scanStatistic" searches for accumulations in the measurements. See Details for descriptions of the techniques.}

  \item{nstart}{
If \code{method="kmeans"}, this is the number of restarts for k-means. See \code{\link{kmeans}} for details.
}
  \item{iter.max}{
If \code{method="kmeans"}, the maximum number of iterations for k-means. See \code{\link{kmeans}} for details.
}

  \item{edge}{If \code{method="edgeDetector"}, this decides which of the edges is used as a threshold for binarization. If set to "firstEdge",the binarization threshold is the first combination of two successive sorted values whose difference exceeds a predefined value (average gradient * \code{scaling}). The parameter \code{scaling} can be used to adjust this value. 
  
  If set to "maxEdge", the binarization threshold is the position of the edge with the overall highest gradient.}
  
  \item{scaling}{If \code{method="edgeDetector"} and \code{edge="firstEdge"}, this holds the scaling factor used for adjustment of the average gradient.}
  
  \item{windowSize}{If \code{method="scanStatistic"}, this specifies the size of the scanning window (see Details). The size is given as a fraction of the whole range of input values for a gene. Default is 0.25.}
  
  \item{sign.level}{If \code{method="scanStatistic"}, the significance level used for the scan statistic (see Details).
}

  \item{dropInsignificant}{
If this is set to true, genes whose binarizations are insignificant in the scan statistic (see Details) are removed from the binarized time series. Otherwise, a warning is printed if such genes exist.
}

}
\details{

This method supports three binarization techniques:
\describe{
\item{k-means clustering}{For each gene, k-means clusterings are performed to determine a good separation of groups. The values belonging to the cluster with the smaller centroid are set to 0, and the values belonging to the greater centroid are set to 1.} 

\item{Edge detector}{This approach first sorts the measurements for each gene. In the sorted measurements, the algorithm searches for differences of two successive values that satisfy a predefined condition: If the "firstEdge" method was chosen, the pair of values whose difference exceeds the scaled average gradient of all values is chosen and used as maximum and minimum value of the two groups. If the "maxEdge" method was chosen, the largest difference between two successive values is taken. For details, see Shmulevich et al.}

\item{Scan statistic}{The scan statistic assumes that the measurements for each gene are uniformly and independently distributed independently over a certain range. The scan statistic shifts a scanning window across the data and decides for each window position whether there is an unusual accumulation of data points based on an approximated test statistic (see Glaz et al.). The window with the smallest p-value is remembered. The boundaries of this window form two thresholds, from which the value that results in more balanced groups is taken for binarization. Depending on the supplied significance level, gene binarizations are rated according to the p-value of the chosen window.}

}
}
\value{
Returns a list with the following elements:
\item{binarizedMeasurements}{A list of matrices with the same structure as \code{measurements} containing the binarized time series measurements}
\item{reject}{If \code{method="scanStatistic"}, a Boolean vector indicating for each gene whether the scan statistic algorithm was able to find a significant binarization window (FALSE) or not (TRUE). Rejected genes should probably be excluded from the data.}
\item{thresholds}{The thresholds used for binarization}
}
\references{
I. Shmulevich and W. Zhang (2002), Binary analysis and optimization-based normalization of gene expression data. Bioinformatics 18(4):555--565.

J. Glaz, J. Naus, S. Wallenstein (2001), Scan Statistics. New York: Springer.
}
\seealso{
\code{\link{reconstructNetwork}}
}
\examples{
# load test data
data(yeastTimeSeries)
			
# perform binarization with k-means
bin <- binarizeTimeSeries(yeastTimeSeries)
print(bin)

# perform binarization with scan statistic
# - will find and remove 2 insignificant genes!
bin <- binarizeTimeSeries(yeastTimeSeries, method="scanStatistic",
                          dropInsignificant=TRUE, sign.level=0.2)
print(bin)

# perform binarization with edge detector
bin <- binarizeTimeSeries(yeastTimeSeries, method="edgeDetector")
print(bin)

# reconstruct a network from the data
reconstructed <- reconstructNetwork(bin$binarizedMeasurements,
                                    method="bestfit", maxK=4)
print(reconstructed)
}
\keyword{Boolean network
	binarization
	binarize
	k-means
	edge detector
	scan statistic}
