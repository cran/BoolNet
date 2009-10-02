
# Encode a vector of binary values <bin> with <len> bits
# to a decimal number
bin2dec <- function(bin,len)
{
	if (len %% 32 == 0)
		numElts <- as.integer(len / 32)
	else
		numElts <- as.integer(len / 32) + 1

	dec = rep(0,numElts)
	
	dec = .C("bin2dec",as.integer(dec),as.integer(bin),as.integer(len))[[1]]
}

# Decode the <len> low-order bits of <dec> to a vector of binary values,
# where the first entry is the least significant bit
dec2bin <- function(dec,len)
{
	bin = rep(0,len)
	
	bin = .C("dec2bin",as.integer(bin),as.integer(dec),as.integer(len))[[1]]
}

# Generate a list of all assignments of n variables with N possible values
allcombn <- function(N,n)
{
	rownum = N^n
	sapply(n:1,function(i)
		{
      			rep(1:N,each=N^(i-1),len=rownum)
      		})
}

# Get DNF representation of a truth table <truthTable>
# using the gene names in <genes>
getDNF <- function(truthTable,genes)
{
	# check for constant functions
	if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
		return("0")
	else
	if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
		return("1")
	
	# build canonical DNF from 1 entries in truth table	
	entries <- allcombn(2,length(genes)) - 1
	conjunctions <- sapply(1:nrow(entries),function(row)
		{
			conj <- entries[row,]
			if (truthTable[row])
			{
				paste("(",paste(sapply(1:length(conj),function(lit)
					{
						if (conj[lit])
							genes[lit]
						else
							paste("!",genes[lit],sep="")
					}),collapse=" & "),")",sep="")
			}
			else
				""
		})
		paste(conjunctions[conjunctions != ""],collapse = " | ")		
}

# Retrieves a string representation of an interaction function by either
# building a DNF (if <readableFunction> is false)
# or returning an unspecific function description (if <readableFunction> is true).
# <truthTable> contains the result column of the interaction's truth table, and
# <genes> contains the names of the involved genes.
getInteractionString <- function(readableFunctions,truthTable,genes)
{
	if (readableFunctions)
		getDNF(truthTable,genes)
	else
	{
		if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
			return("0")
		else
		if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
			return("1")
		else
			paste("<f(",
				paste(genes,collapse=","),"){",
				paste(truthTable,collapse=""),"}>", sep="")
	}		
}

# Retrieves an internal graph representation of a state table
# in <attractorInfo>
# The vertices are the states, and the edges are the transitions.
# Returns a list with a $vertices vector and a 2-column $edges matrix.
getStateGraphStructure <- function(attractorInfo)
{	
	numRows <- ncol(attractorInfo$stateInfo$table)	
	fixedGenes <- which(attractorInfo$stateInfo$fixedGenes != -1)
	nonFixedGenes <- which(attractorInfo$stateInfo$fixedGenes == -1)
	
	if (is.null(attractorInfo$stateInfo$initialStates))
	{
		inputStates <- matrix(ncol=length(attractorInfo$stateInfo$genes),nrow=2^length(nonFixedGenes))
	
		# encode the initial states
		# first, encode the changing part by calculating all combinations
		temp <- allcombn(2,length(nonFixedGenes)) - 1
		inputStates[,nonFixedGenes] <- temp[,ncol(temp):1]
	
		if (length(fixedGenes) > 0)
		# if there are fixed genes, encode their values
			inputStates[,fixedGenes] <- sapply(attractorInfo$stateInfo$fixedGenes[fixedGenes],
							     function(value)rep(value,2^length(nonFixedGenes)))	
	}
	else
	{
		inputStates <- t(apply(attractorInfo$stateInfo$initialStates,2,dec2bin,length(attractorInfo$stateInfo$genes)))
	}
	
	inputStates <- apply(inputStates,1,function(state)paste(state,collapse=""))
	
	# build "hash table" of state numbers
	stateHash <- 1:length(inputStates)
	names(stateHash) <- inputStates
	
	
	# encode output states
	outputStates <- apply(attractorInfo$stateInfo$table,2,function(state)
				{
					bin <- paste(dec2bin(state,length(attractorInfo$stateInfo$genes)),collapse="")
				})
	return(list(vertices=inputStates,edges=cbind(1:length(inputStates),
			sapply(outputStates,function(state)stateHash[state]))))
}

# Reorders a matrix of states <stateMatrix> (with each column being one state)
# in a canonical way such that the "smallest" state is the first.
# This makes attractor representations unique.
canonicalStateOrder <- function(stateMatrix)
{
	smallestIndex <- -1
	smallestVal <- rep(Inf,nrow(stateMatrix))
	for (i in 1:ncol(stateMatrix))
	# iterate over states
	{
		for (j in 1:nrow(stateMatrix))
		# iterate over elements of encoded state
		{
			if (stateMatrix[j,i] < smallestVal[j])
			# determine new minimum
			{
				smallestVal <- stateMatrix[,i]
				smallestIndex <- i
				break
			}
			else
			{
				if (stateMatrix[j,i] > smallestVal[j])
					break
			}
		}
	}
	if (smallestIndex != 1)
		# rearrange matrix
		return(cbind(stateMatrix[,smallestIndex:ncol(stateMatrix),drop=FALSE],
			     stateMatrix[,(1:(smallestIndex-1)),drop=FALSE]))
	else
		return(stateMatrix)
}
