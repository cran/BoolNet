
chooseNetwork <- function(networkCollection,functionIndices)
{
	stopifnot(inherits(networkCollection,"BooleanNetworkCollection"))
	interactions <- mapply(function(interaction,index)
				{
					interaction[[index]]
				},
				networkCollection$interactions,functionIndices,SIMPLIFY=FALSE)
	res <- list(genes=networkCollection$genes,
		    interactions=interactions,
		    fixed=networkCollection$fixed)
	class(res) <- "BooleanNetwork"
	return(res)
}
