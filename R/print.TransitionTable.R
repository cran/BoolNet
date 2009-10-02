# Custom print function for class TransitionTable
print.TransitionTable <- function(x, quote = FALSE, max.levels = NULL,width = getOption("width"), ...)
{
	numGenes <- (ncol(x) - 2) / 2
	cat(format("State",width=max(7,numGenes),justify="right"),"   ",
	    format("Next state",width=max(12,numGenes + 3),justify="right"),
	    format("Attr. basin",width=13,justify="right"),
	    format("# trans. to attr.",width=19,justify="right"),"\n",sep="")
	colIndices <- c(1,numGenes,numGenes + 1, 2*numGenes, 
			2*numGenes + 1, 2*numGenes + 2)
	apply(x,1,function(row)
	{
		# paste all states of input and output into one string, and put out all columns of the table in a
		# formatted way
		cat(format(paste(row[colIndices[1]:colIndices[2]],collapse=""),width=max(7,numGenes),justify="right"),
		    " =>",
		    format(paste(row[colIndices[3]:colIndices[4]],collapse=""),width=max(12,numGenes + 3),justify="right"),
		    format(row[colIndices[5]],width=13,justify="right"),
    		    format(row[colIndices[6]],width=19,justify="right"),"\n",sep="")
	})	
	return(invisible(x))
}
