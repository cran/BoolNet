# Sets the genes in <fixIndices> to the values in <values>
# and returns the customized copy of <network>
fixGenes <- function(network,fixIndices,values)
{
  stopifnot(inherits(network,"BooleanNetwork") | inherits(network,"ProbabilisticBooleanNetwork") | inherits(network,"BooleanNetworkCollection"))

  if (length(fixIndices) != length(values) && length(values) != 1)
    stop("fixIndices and values must have the same number of elements, or values must have 1 element!")

  if (any(is.na(network$fixed[fixIndices])))
    stop("fixIndices contains invalid indices!")
    
  if (any(values != 0 & values != 1 & values != -1))
    stop("Please supply only 0, 1, or -1 in values!")

  network$fixed[fixIndices] <- values
  return(network)
}
