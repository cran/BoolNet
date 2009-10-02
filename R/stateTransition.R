stateTransition <- function(network,state)
{
	nonFixedIndices = (network$fixed == -1)

	res = rep(-1,length(network$genes))

	res[nonFixedIndices] <- sapply(which(nonFixedIndices),function(i)
		{    
			if(network$interactions[[i]]$input[1] == 0)
			# this is a constant gene with no transition function
				return(state[i])
			input = state[network$interactions[[i]]$input]
			return(network$interactions[[i]]$func[bin2dec(rev(input),length(input)) + 1])
		})

	res[!nonFixedIndices] = network$fixed[!nonFixedIndices]
	return(res)
}
