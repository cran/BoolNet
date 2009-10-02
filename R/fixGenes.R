# Sets the genes in <fixIndices> to the values in <values>
# and returns the customized copy of <network>
fixGenes <- function(network,fixIndices,values)
{
	stopifnot(inherits(network,"BooleanNetwork"))
	if (length(fixIndices) != length(values) && length(values) != 1)
		stop("fixIndices and values must have the same number of elements, or values must have 1 element!")
	network$fixed[fixIndices] <- values
	return(network)
}
