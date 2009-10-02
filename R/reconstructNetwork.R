# Reconstruct a Boolean network from a transition table or a list of time series in <measurements>.
# If <method> is "bestfit", Lähdesmäki's best-fit extension algorithm is called.
# If <method> is "reveal", Liang's REVEAL algorithm is called.
# <maxK> specifies the maximum number of input genes of a function.
# If <readableFunctions> is true, DNF representations of the functions are generated.
reconstructNetwork <- function(measurements,method=c("bestfit","reveal"),maxK=5,readableFunctions=FALSE)
{
	# determine method to use
	meth <- switch(match.arg(method,c("bestfit","reveal")),
			bestfit=0,
			reveal=1,
			stop("'method' must be one of \"bestfit\",\"reveal\""))
	if (inherits(measurements,"TransitionTable"))
	# preprocess transition table and call algorithm
	{
		numGenes <- (ncol(measurements) - 2) / 2
		
		if (numGenes < maxK)
		{
		  maxK <- numGenes
		  warning(paste("maxK was chosen greater than the total number of input genes and reset to ",numGenes,"!",sep=""))
		}
		
		# call C code
		res <- .Call("reconstructNetwork_R",
			    as.integer(t(as.matrix(measurements[,1:numGenes]))),
			    as.integer(t(as.matrix(measurements[,(numGenes+1):(2*numGenes)]))),
			    as.integer(nrow(measurements)),
			    as.integer(maxK),
			    as.integer(meth))
		genenames <- sapply(colnames(measurements)[1:numGenes],function(x)strsplit(x,".",fixed=TRUE)[[1]][2])
	}
	else
	# the measurements are time series
	{
		if (!is.null(dim(measurements)))
		# only one time series => create list
			measurements <- list(measurements)
			
		numGenes <- nrow(measurements[[1]])
		
		if (numGenes < maxK)
		{
		  maxK <- numGenes
		  warning(paste("maxK was chosen greater than the total number of input genes and reset to ",numGenes,"!",sep=""))
		}

		
		genenames <- rownames(measurements[[1]])
		if (is.null(genenames))
		  genenames <- paste("Gene",1:numGenes)
		inputStates <- c()
		outputStates <- c()
		
		for (measurement in measurements)
		# iterate over all time series and build state vectors
		{
			if (numGenes != nrow(measurement))
				stop("All measurement matrices must contain the same genes!")
			inputStates <- c(inputStates,as.integer(as.matrix(measurement[,1:(ncol(measurement)-1)])))
			outputStates <- c(outputStates,as.integer(as.matrix(measurement[,2:ncol(measurement)])))
			
		}
		
		# call C code
		res <- .Call("reconstructNetwork_R",
			    inputStates,
			    outputStates,
			    as.integer(length(inputStates) / numGenes),
			    as.integer(maxK),
			    as.integer(meth))
		
	}
	
	if (any(sapply(res,function(interaction)length(interaction)==0)))
	# some function lists are empty
		warning("Some functions could not be inferred. Possibly the input data is noisy or maxK was chosen too small!")
	
	# prepare result object
	res <- list(genes=genenames,
		    interactions=lapply(res,function(gene)
		    			lapply(gene,function(interaction)
		    			{
		    				interaction$expression <- 
		    					getInteractionString(readableFunctions,
							        	     interaction$func,
									    genenames[interaction$input])
						interaction
		    			})),
		    fixed=sapply(res,function(gene)
		    		{
		    			if (length(gene) == 0)
		    				-1
		    			else
		    			if (gene[[1]]$input[1] == 0)
		    				gene[[1]]$func[1]
		    			else
		    				-1
		    		}))
	
	names(res$interactions) <- res$genes
	
	class(res) <- "BooleanNetworkCollection"
	return(res)
}
