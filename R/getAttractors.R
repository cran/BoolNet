# Identify attractors in a Boolean network.
# <network> is a BooleanNetwork structure specifying the network.
# <method> specifies what kind of search is conducted: "exhaustive" performs a full search over all states,
# "random" generates <startStates> random initial states, and "chosen" uses the states supplied in <startStates>.
# <genesON> and <genesOFF> are lists of genes to be set to ON/1 or OFF/0 respectively.
# If <canonical> is true, states in the attractors are reordered such that the "smallest" state is the first
getAttractors <- function (network, method=c("exhaustive","random","chosen"), startStates=list(),
			   genesON = c(), genesOFF = c(), canonical=TRUE) 
{
	stopifnot(inherits(network,"BooleanNetwork"))
	
	if (length(method) == 3)
	# if no method is supplied, infer method from the type of <startStates>
	{
		if (is.numeric(startStates))
			method <- "random"
		else
		if (is.list(startStates) & (length(startStates) > 0))
			method = "chosen"
		else
			method = "exhaustive"
	}
			
	if ((length(network$genes) > 29) & (match.arg(method) == "exhaustive"))
		stop(paste("An exhaustive attractor search is not supported for networks with more than 29 genes!",
		 	   "Please use a heuristic method or provide start states!"))
		
	fixedPositions <- which(network$fixed != -1)
	nonFixedPositions <- which(network$fixed == -1)
	
	startStates <- switch(match.arg(method,c("exhaustive","random","chosen")),
		exhaustive = list(),
		random = {
				if (!is.numeric(startStates))
					stop("Please supply the number of random states in startStates!")

				if (startStates > (2 ^ length(nonFixedPositions)))
				# more start states than in the full network
				{
					warning(paste("The number of random states is set larger than the total",
						      "number of states. Performing an exhaustive search!"))
					list()					
				}
				else
				# generate random matrix
				{
					mat <- matrix(nrow=startStates,ncol=length(network$genes))
					
					if (length(fixedPositions) != 0)
						# fill fixed positions with the corresponding values
						mat[,fixedPositions] <- sapply(fixedPositions,function(x)
										rep(network$fixed[x],startStates))
					
					# generate other positions randomly
					mat[,nonFixedPositions] <- round(runif(n=startStates*length(nonFixedPositions)))
					
					# eliminate duplicates
					mat <- unique(mat)
					
					while (nrow(mat) != startStates)
					# if duplicates were removed, generate new states until the
					# desired number of states is reached
					{
						vec <- rep(0,length(network$genes))
						if (length(fixedPositions) != 0)
							# fill fixed positions
							vec[fixedPositions] <- sapply(fixedPositions,
										      function(x)network$fixed[x])
						# generate other positions randomly
						vec[nonFixedPositions] <- round(runif(n=length(nonFixedPositions)))
						mat <- unique(rbind(matrix,vec))
					}
					cat("Using initial states:\n")
					
					# print states and form a list
					res <- lapply(1:nrow(mat),function(i)
						{
							cat(paste(mat[i,],collapse=""),"\n",sep="")
							mat[i,]
						})
					cat("\n")
					res
				}
			 },
		chosen = {
				if (!is.list(startStates) | length(startStates) == 0)
					stop("No start states supplied!")
				if (!all(sapply(startStates,function(x)(length(x) == length(network$genes)))))
					stop(paste("Please provide binary vectors with",length(network$genes),
						"elements in startStates!"))
				startStates
			 }
	)
	
	specialInitialization <- NULL
	fixedPositions <- which(network$fixed != -1)
	
	convertedStartStates <- NULL
	
	if (length(startStates) > 0)
		convertedStartStates <- sapply(startStates,function(x)bin2dec(x,length(network$genes)))
	
	# build list of special initializations from genesON and genesOFF
	if (length(genesON) > 0) 
	{
		onPositions <- match(tolower(genesON), network$genes)
		# correct the indices for fixed positions (these genes are omitted)
		sapply(fixedPositions, function(x) 
		{
			idx <- (onPositions > x)
			onPositions[idx] <<- onPositions[idx] - 1
		})
		onPositions <- onPositions - 1
		# build a two-row matrix with the gene index in the first row
		# and the value in the second row
		specialInitialization <- rbind(onPositions,rep(1, length(genesON)))
	}
	if (length(genesOFF) > 0) 
	{
		offPositions <- match(tolower(genesOFF), network$genes)
		# correct the indices for fixed positions (these genes are omitted)
		sapply(fixedPositions, function(x) 
		{
		    idx <- (offPositions > x)
		    offPositions[idx] <<- offPositions[idx] - 1
		})
		offPositions <- offPositions - 1
		specialInitialization <- cbind(specialInitialization,
					  rbind(offPositions,rep(0,length(genesOFF))))
	}

	# the C code requires all interactions to be coded into one vector:
	# Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
	inputGenes <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$input)))
	inputGenePositions <- as.integer(cumsum(c(0,sapply(network$interactions,
					 function(interaction)length(interaction$input)))))

	# Do the same for the transition functions.
	transitionFunctions <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$func)))
	transitionFunctionPositions <- as.integer(cumsum(c(0,sapply(network$interactions,
						  function(interaction)length(interaction$func)))))

	# Call the C code
	result <- .Call("getAttractors_R",inputGenes,inputGenePositions,
		    transitionFunctions,transitionFunctionPositions,
		    as.integer(network$fixed),as.integer(specialInitialization),
		    as.integer(convertedStartStates),
		    PACKAGE="BoolNet")
	
	if (is.null(result))
		stop("An error occurred in external C code!")
	
	
	if (length(network$genes) %% 32 == 0)
		numElementsPerEntry <- as.integer(length(network$genes) / 32)
	else
		numElementsPerEntry <- as.integer(length(network$genes) / 32  + 1)
	
	result$stateInfo$table <- matrix(result$stateInfo$table,nrow=numElementsPerEntry)
	
	if (!is.null(result$stateInfo$initialStates))
		result$stateInfo$initialStates <- matrix(result$stateInfo$initialStates,nrow=numElementsPerEntry)
	
	for (i in 1:length(result$attractors))
	{
		result$attractors[[i]]$involvedStates <- matrix(result$attractors[[i]]$involvedStates,nrow=numElementsPerEntry)
		if (canonical)
		# reorder states
			result$attractors[[i]]$involvedStates <- canonicalStateOrder(result$attractors[[i]]$involvedStates)
	}
	
	# extend the resulting structure by additional information, and assign a class	    
	result$stateInfo$genes <- network$genes
	result$stateInfo$fixedGenes <- network$fixed

	class(result$stateInfo) <- "BooleanStateInfo"
	class(result) <- "AttractorInfo"
	return(result)
}
