getPathToAttractor <- function(network, state, includeAttractorStates=c("all","first","none"))
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
  
  includeAttractorStates <- match.arg(includeAttractorStates, c("all","first","none"))
  
  numGenes <- (ncol(table) - 2) / 2
  initialStates <- apply(table[,1:numGenes,drop=FALSE],1,function(x)paste(x,collapse=""))
  
  currentState <- state
  
  res <- data.frame(matrix(state,nrow=1))
  
  stateCount <- 1
  repeat
  {
    currentStateIdx <- which(initialStates == paste(currentState,collapse=""))

    if (length(currentStateIdx) == 0)
      stop(paste("Could not find state",paste(currentState,collapse=""),"in the transition table!"))

    # stop depending on "includeAttractorStates" option
    if ((includeAttractorStates == "all" && stateCount == nrow(table))
        || (includeAttractorStates == "first" && table[currentStateIdx,"transitionsToAttractor"] == 0)
        || (includeAttractorStates == "none" && table[currentStateIdx,"transitionsToAttractor"] <= 1))
      break
      
    currentState <- as.integer(table[currentStateIdx,(numGenes+1):(2*numGenes)])
    res <- rbind(res,currentState)
    
    stateCount <- stateCount + 1
  }
  
  # special case: start state is attractor state and we do not want to include attractor states
  # => return empty data frame
  if (includeAttractorStates == "none" && table[currentStateIdx,"transitionsToAttractor"] == 0)
  {
    res <- data.frame(matrix(nrow=0,ncol=numGenes))
  }

  colnames(res) <- sapply(colnames(table)[1:numGenes],function(n)strsplit(n,".",fixed=TRUE)[[1]][2])
  
  rownames(res) <- NULL
  return(res)
}
