# Custom print function for class AttractorInfo
print.AttractorInfo <- function(x, quote = FALSE, max.levels = NULL,width = getOption("width"), ...)
{
	numGenes <- length(x$stateInfo$genes)
	
	# decode binary representation of states involved in the attractors
	binMatrices <- lapply(x$attractors,function(attractor)
			{
				res <- matrix(apply(attractor$involvedStates,2,function(state)
					dec2bin(state,numGenes)),nrow=numGenes)
			})
	
	# order attractors according to their lengths
	attractorLengths <- sapply(x$attractors,function(attractor)nrow(attractor$involvedStates))	
	reordering <- order(attractorLengths)
	binMatrices <- binMatrices[reordering]
	attractors <- x$attractors[reordering]
	
	lapply(1:length(attractors),function(i)
	{
		# print general information on the attractor
		cat("Attractor ",i," consists of ",ncol(attractors[[i]]$involvedStates),
		    " state(s) and has a basin of ",attractors[[i]]$basinSize, " state(s):\n\n",sep="")
		cat("Genes are encoded in the following order: ",paste(x$stateInfo$genes,collapse=" "),"\n\n",sep="")
		
		# print a graphical representation of the attractor cycle
		cat(" |--<",paste(rep("-",numGenes-1),collapse=""),"|\n",sep="")
		cat(" V ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
		apply(binMatrices[[i]],2,function(col)
		{
			cat(" ",col,"   |\n",sep="")
			cat(" | ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
			
		})
		cat(" V ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
		cat(" |-->",paste(rep("-",numGenes-1),collapse=""),"|\n\n",sep="")
	})	
	return(invisible(x))
}
