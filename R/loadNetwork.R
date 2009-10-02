# Load a network in a specified rule description language
# from file <file>.
# <bodySeparator> is the character that separates targets and factors
loadNetwork <- function(file,bodySeperator=",") 
{
	op <- c("!","&", "\\|", "\\(", "\\)")

	func <- tolower(readLines(file,-1)[-1])

	# / in a gene name disturbs parsing
	func <- gsub("/","_",func) 

	tmp <-  unname(sapply(func,function(x){strsplit(x,bodySeperator)[[1]]}))
	targets <- tmp[1,]
	factors <- tmp[2,]

	factors.tmp <- sapply(factors,function(x)
	# extract gene names from Boolean expressions
	{
	sapply(op,function(y)
	{
	x <<- gsub(y," ",x)
	})
	tmp <- strsplit(x," ")[[1]]

	tmp <- unique(tmp[tmp != ""])
	})

	# create list of all gene names in both sides of the functions
	genes <- unique(c(targets,unname(unlist(factors.tmp))))

	# extract "real" gene names from the list, drop constants
	suppressWarnings(genes <- genes[is.na(as.integer(genes))])

	fixed <- rep(-1,length(genes))

	interactions <- NULL

	for(i in 1:length(targets))
	{
		inputGenes <- factors.tmp[[i]]
		if(suppressWarnings(is.na(as.integer(inputGenes[1]))))
		# the input is not a number
		{
			inputIndices = match(inputGenes,genes)
			exp = as.matrix(allcombn(2,length(inputIndices)) - 1)
			for(j in 1:length(inputIndices))
		      	{
				tmp1 = paste(exp[,j],collapse=",")
				eval(parse(text = paste(inputGenes[j],"=c(",tmp1,")")))
		      	}
			interactions[[i]] = list(input = inputIndices,
						 func = as.numeric(eval(parse(text=factors[i]))),
						 expression = factors[i])
		  }
		  else
		  # this is a constant gene
		  {
		  	fixed[i] = as.integer(inputGenes)
		  	interactions[[i]] = list(input=0,func=fixed[i],expression = fixed[i])
		  }        
	}

	if(length(targets) < length(genes))
	# some genes are only used as inputs, but are not involved in the network
	# -> create dummy input and function
	{
		for(i in (length(targets)+1):length(genes))
		{
	  		interactions[[i]] = list(input = 0,func=-1)
		}
	}
	
	names(interactions) <- genes
	
	res <- list(interactions = interactions,
		    genes = genes,
		    fixed = fixed)
	class(res) <- "BooleanNetwork"
	return(res)
}

