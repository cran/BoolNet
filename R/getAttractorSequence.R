
getAttractorSequence <- function(attractorInfo, attractorNo)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))
  
  if (missing(attractorNo) || attractorNo <= 0 || attractorNo > length(attractorInfo$attractors))
    stop("Please provide a valid attractor number!")
  
  if (!is.null(attractorInfo$attractors[[attractorNo]]$initialStates))
    stop("A sequence can be obtained for synchronous attractors only!")
  
  numGenes <- length(attractorInfo$stateInfo$genes)
  
  # decode binary representation of states involved in the attractors
  binMatrix <- t(matrix(apply(attractorInfo$attractors[[attractorNo]]$involvedStates,2,function(state)
                        dec2bin(state,numGenes)),nrow=numGenes))
                        
  colnames(binMatrix) <- attractorInfo$stateInfo$genes
  return(as.data.frame(binMatrix))
}
