# Determine the basin of attraction of attractor <attractorNo> in <attractorInfo>.
getBasinOfAttraction <- function(attractorInfo,attractorNo)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))

  if (is.null(attractorInfo$stateInfo$table))
    stop(paste("This AttractorInfo structure does not contain transition table information.",
           "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))
  
  if (missing(attractorNo) || attractorNo < 0 || attractorNo > length(attractorInfo$attractors))
    stop("Please provide a valid attractor number!")
  
  table <- getTransitionTable(attractorInfo)
  return(table[which(table$attractorAssignment == attractorNo),,drop=FALSE])
}
