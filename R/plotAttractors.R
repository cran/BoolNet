# Plot state tables of all attractors in <attractorInfo>.
# Genes are grouped according to <grouping>.
# An additional title can be supplied in <title>.
# If <plotFixed> is set, fixed variables are included in the plot.
# <onColor> and <offColor> specify the colors of ON/1 and OFF/0 states.
plotAttractors <- function (attractorInfo, grouping = list(), title = "", plotFixed = TRUE,
			onColor="green",offColor="red") 
{
	stopifnot(inherits(attractorInfo,"AttractorInfo"))
	numGenes <- length(attractorInfo$stateInfo$genes)
	
	# determine list of genes to be plotted
	whichFixed <- which(attractorInfo$stateInfo$fixed != -1)
	if (plotFixed | (length(whichFixed) > 0))
		plotIndices <- 1:numGenes
	else
		plotIndices <- (1:numGenes)[-whichFixed]
	
	# convert decimal state numbers to binary state matrices (one for each attractor)
	binMatrices <- lapply(attractorInfo$attractors,function(attractor)
					{
						res <- matrix(apply(attractor$involvedStates,2,function(state)
							dec2bin(state,numGenes)[plotIndices]),nrow=length(plotIndices))
					})

	# count the numbers of attractors with equal lengths
	attractorLengths <- sapply(attractorInfo$attractors,function(attractor)ncol(attractor$involvedStates))	
	lengthTable <- table(attractorLengths)
	
	res <- lapply(1:length(lengthTable),function(i)
	# accumulate all attractors with equal length in one matrix and plot them
	{
		len <- as.integer(names(lengthTable)[i])
		cnt <- lengthTable[i]

		# initialize with empty plot
		plot(c(),c(),xlim=c(0,len*cnt),ylim=c(-2,length(plotIndices)+1),xlab="",ylab="",
		     axes=FALSE,main=paste(title, "Attractors with ",len," state(s)",sep=""))
		
		# build accumulated matrix     
		totalMatrix <- c()
		for (mat in binMatrices[which(attractorLengths == len)])
		{
			totalMatrix <- cbind(totalMatrix,mat)
		}
		rownames(totalMatrix) <- attractorInfo$stateInfo$genes[plotIndices]
		colnames(totalMatrix) <- sapply(1:cnt,function(i)paste("Attr",i,".",1:len,sep=""))
		
		if(length(grouping)>0)
	 		# reorder genes according to the supplied groups
			totalMatrix = totalMatrix[unlist(grouping$index),]
		
		par(yaxt='s',las=2)
    		axis(2,(1:length(plotIndices))-0.5,rownames(totalMatrix))

		# plot active and inactive states
		for(i in 1:ncol(totalMatrix))
			for(j in 1:nrow(totalMatrix))
			{
				if(totalMatrix[j,i] == 1)
					rect(i-1,j-1,i,j,col=onColor,border="gold")
				else
					rect(i-1,j-1,i,j,col=offColor,border="gold")
			}
		
		# draw vertical separators between attractors	
	    	sep = seq(0,ncol(totalMatrix),by=len)    
	    	abline(v = sep[-1],col="white",lwd=3)
	    	
	    	# output frequency of attractor (basin size / number of states)
		freq <- sapply(attractorInfo$attractors[which(attractorLengths == len)],
				function(attractor)attractor$basinSize/length(attractorInfo$stateInfo$table)) * 100

		text(sep[1:(length(sep)-1)] + len/2, rep(0.4+nrow(totalMatrix),ncol(totalMatrix)),
			paste(freq,"%",sep=""),cex=.75,font=3)
			
		if(length(grouping)>0)
		# draw separators between groups, and print group names
		{
			sepPos = cumsum(sapply(grouping$index,length))
			abline(h=sepPos[-length(sepPos)],col="black",lwd=3)
			text(ncol(totalMatrix)/2,sepPos-0.5,grouping$class,cex=0.9)
      		}

	  	legend(x="bottomright",pch=c(15,15),col=c(onColor,offColor),legend = c("active","inactive"),cex=0.7,horiz=T)
		totalMatrix
	})
	
	# return a list of accumulated matrices
	names(res) <- names(lengthTable)
	return(res)

}
