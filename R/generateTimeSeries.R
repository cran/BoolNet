
# Generate time series from a Boolean network <network>.
# <numSeries> is the number of start states for which time series are generated.
# <numMeasurements> is the number of time points for each time series.
# <type> is the type of transitions used (synchronous or asynchronous)
# If <type> is "asynchronous", <geneProbabilities> describes the 
# transition probabilities for the genes.
# If <noiseLevel> is not 0, Gaussian noise is added to the result.
generateTimeSeries <- function(network, numSeries, numMeasurements, 
                            type = c("synchronous","asynchronous","probabilistic"),
                            geneProbabilities, 
                            noiseLevel = 0.0)
{
  if (missing(geneProbabilities))
    geneProbabilities <- NULL
    
  ts <- lapply(seq_len(numSeries), function(i)
  {
    startState <- round(runif(length(network$genes)))
    res <- startState
    for (j in 2:numMeasurements)
    {
      startState <- stateTransition(network, startState, 
                                    type=type, geneProbabilities=geneProbabilities)
      res <- cbind(res,startState)
    }
    
    colnames(res) <- NULL

    if (noiseLevel != 0)
    {
      res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))
    }
    
    rownames(res) <- network$genes
    colnames(res) <- seq_len(ncol(res))
    return(res)
  })
  return(ts)
}
