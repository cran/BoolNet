# Custom print function for class BooleanNetwork
print.BooleanNetworkCollection <- function(x, quote = FALSE, max.levels = NULL,width = getOption("width"), ...)
{
	cat("Multiple solutions for a Boolean network with",length(x$genes),"genes\n\n")
	cat("Involved genes:\n",paste(x$genes,collapse=" "),"\n\n",sep="")
	cat("Transition functions:\n")
		
	mapply(function(gene,interaction)
		{
			cat("\nAlternative transition functions for gene ",gene,":\n",sep="")
			# print original expressions read from the files (if available)
			lapply(interaction,function(func)
			{
				cat(gene," = ",func$expression," (error: ",func$error,")\n",sep="")
			})
		},	
		x$genes,x$interactions)

	if (sum(x$fixed != -1) > 0)
	{
		cat("\nKnocked-out and over-expressed genes:\n")
		mapply(function(gene,fixed)
			{
				if (fixed != -1)
					cat(gene," = ",fixed,"\n",sep="")
			},
			x$genes,x$fixed)
	}	
		
	return(invisible(x))
}
