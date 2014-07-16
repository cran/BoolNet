\name{generateTimeSeries}
\alias{generateTimeSeries}

\title{
Generate time series from a network
}
\description{
Generates time series by simulating successive state transitions from random start states. In addition, the resulting matrices can be perturbed by Gaussian noise.
}
\usage{
generateTimeSeries(network, 
                   numSeries, 
                   numMeasurements, 
                   type = c("synchronous","asynchronous","probabilistic"),
                   geneProbabilities, 
                   perturbations = 0,
                   noiseLevel = 0)
}

\arguments{
  \item{network}{
An object of class \code{BooleanNetwork} or \code{SymbolicBooleanNetwork} that contains the network for which time series are generated
}

  \item{numSeries}{
The number of random start states used to generate successive series of states, that is, the number of time series matrices to generate
}

  \item{numMeasurements}{
The number of states in each of the time series matrices. The first state of each time series is the randomly generated start state. The remaining \code{numMeasurements - 1} states are obtained by successive state transitions.
}

  \item{type}{
The type of state transitions to be performed (see \code{\link{stateTransition}})
}

\item{geneProbabilities}{
An optional vector of probabilities for the genes if \code{type="asynchronous"}. By default, each gene has the same probability to be chosen for the next state transition. These probabilities can be modified by supplying a vector of probabilities for the genes which sums up to one.
}

  \item{perturbations}{
If this argument has a value greater than 0, artificial perturbation experiments are generated. That is, \code{perturbations} genes in each time series are knocked out or overexpressed artificially using the \code{\link{fixGenes}} function.   
}

  \item{noiseLevel}{
If this is non-zero, it specifies the standard deviation of the Gaussian noise which is added to all entries of the time series matrices. By default, no noise is added to the time series.
}
}
\value{
A list of matrices, each corresponding to one time series. Each row of these matrices contains measurements for one gene on a time line, i. e. column \code{i+1} contains the successor states of column \code{i+1}. If \code{noiseLevel} is non-zero, the matrices contain real values, otherwise they contain only 0 and 1. 

If \code{perturbations>0}, the result list contains an additional matrix \code{perturbations} specifying the artificial perturbations applied to the different time series. This matrix has \code{numSeries} columns and one row for each gene in the network. A matrix entry is 0 for a knock-out of the corresponding gene in the corresponding time series, 1 for overexpression, and NA for no perturbation.

The result format is compatible with the input parameters of \code{\link{binarizeTimeSeries}} \if{latex}{\cr} and \code{\link{reconstructNetwork}}.
}

\seealso{
\code{\link{stateTransition}}, \code{\link{binarizeTimeSeries}}, \code{\link{reconstructNetwork}}
}
\examples{
# generate noisy time series from the cell cycle network
data(cellcycle)
ts <- generateTimeSeries(cellcycle, numSeries=50, numMeasurements=10, noiseLevel=0.1)

# binarize the noisy time series
bin <- binarizeTimeSeries(ts, method="kmeans")$binarizedMeasurements

# reconstruct the network
print(reconstructNetwork(bin, method="bestfit"))
}

\keyword{time series
        Boolean network
      	synchronous update
      	asynchronous update }
