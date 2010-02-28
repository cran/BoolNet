# Creates a deterministic Boolean network from a Probabilistic Boolean network
# by choosing the functions specified in <functionIndices> from the
# network <probabilisticNetwork>.
chooseNetwork <- function(probabilisticNetwork,functionIndices)
{
  stopifnot(inherits(probabilisticNetwork,"ProbabilisticBooleanNetwork") 
        | inherits(probabilisticNetwork,"BooleanNetworkCollection"))

  if (length(functionIndices) != length(probabilisticNetwork$genes))
    stop("Please provide a vector of function indices for each gene!")

  interactions <- mapply(function(interaction,index)
        {
          list(input=interaction[[index]]$input,
            func=interaction[[index]]$func,
            expression=interaction[[index]]$expression)
        },
        probabilisticNetwork$interactions,functionIndices,SIMPLIFY=FALSE)

  res <- list(genes=probabilisticNetwork$genes,
        interactions=interactions,
        fixed=probabilisticNetwork$fixed)

  class(res) <- "BooleanNetwork"
  return(res)
}
