getPathToAttractor <- function(network, state)
{
  stopifnot(inherits(network,"BooleanNetwork") | inherits(network, "AttractorInfo"))
  
  if (inherits(network,"BooleanNetwork"))
  {
    table <- getTransitionTable(getAttractors(network, startStates=list(state)))
  }
  else
  {
    if (is.null(network$stateInfo$table))
      stop(paste("This AttractorInfo structure does not contain transition table information.",
           "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))
    table <- getTransitionTable(network)
  }
  
  numGenes <- (ncol(table) - 2) / 2
  initialStates <- apply(table[,1:numGenes,drop=FALSE],1,function(x)paste(x,collapse=""))
  
  currentState <- state
  res <- state
  repeat
  {
    currentStateIdx <- which(initialStates == paste(currentState,collapse=""))

    if (length(currentStateIdx) == 0)
      stop(paste("Could not find state",paste(currentState,collapse=""),"in the transition table!"))

    if (table[currentStateIdx,"transitionsToAttractor"] == 0)
      break;

    currentState <- table[currentStateIdx,(numGenes+1):(2*numGenes)]
    res <- rbind(res,currentState)
  }
  class(res) <- "data.frame"  
  colnames(res) <- sapply(colnames(table)[1:numGenes],function(n)strsplit(n,".",fixed=TRUE)[[1]][2])
  rownames(res) <- NULL
  return(res);
}
