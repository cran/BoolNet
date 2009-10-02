# Determine the basin of attraction of attractor <attractorNo> in <attractorInfo>.
getBasinOfAttraction <- function(attractorInfo,attractorNo)
{
	stopifnot(inherits(attractorInfo,"AttractorInfo"))
	table <- getTransitionTable(attractorInfo)
	return(table[which(table$attractorAssignment == attractorNo),,drop=FALSE])
}
