generateState <- function(network, specs, default=0)
{
  stopifnot(inherits(network,"ProbabilisticBooleanNetwork") | inherits(network,"BooleanNetworkCollection")
            | inherits(network,"BooleanNetwork"))
  
  if (!all(names(specs) %in% network$genes))
    stop(paste("Undefined gene names:",
               paste(setdiff(names(specs), network$genes), collapse=", ")))
  
  if (!all(specs %in% c(0,1)))
    stop("Please provide only Boolean values!")
    
 state <- rep(default, length(network$genes))
 names(state) <- network$genes
 state[names(specs)] <- specs
 return(state)                
}
