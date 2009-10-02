# Custom print function for class BooleanStateInfo
print.BooleanStateInfo <- function(x, quote = FALSE, max.levels = NULL,width = getOption("width"), ...)
{
	# Create a TransitionTable object and print it
	cat("Transition table of Boolean network\n\n")
	cat("Genes are encoded in the following order: ",paste(x$genes,collapse=" "),"\n\n",sep="")
	
	transitionTable <- .internal.getTransitionTable(x)
	
	print(transitionTable)
}
