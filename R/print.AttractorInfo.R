# Custom print function for class AttractorInfo
print.AttractorInfo <- function(x, activeOnly = FALSE, ...)
{
  numGenes <- length(x$stateInfo$genes)
  attractors <- x$attractors
  
  lapply(seq_along(attractors),function(i)
  {
    if (is.null(attractors[[i]]$initialStates))
    # simple attractor
    {
      # decode binary representation of states involved in the attractors
      binMatrix <- matrix(apply(attractors[[i]]$involvedStates,2,function(state)
                          dec2bin(state,numGenes)),nrow=numGenes)
    
      # print general information on the attractor
      if (is.na(attractors[[i]]$basinSize))
        cat("Attractor ",i," is a simple attractor consisting of ",ncol(attractors[[i]]$involvedStates),
            " state(s)",sep="")
      else
        cat("Attractor ",i," is a simple attractor consisting of ",ncol(attractors[[i]]$involvedStates),
            " state(s) and has a basin of ",attractors[[i]]$basinSize, " state(s)",sep="")
      
      # print a graphical representation of the attractor cycle
      if (activeOnly)
      {
        cat(".\nActive genes in the attractor state(s):\n")
        
        i <- 0
        apply(binMatrix,2,function(col)
        {
          i <<- i + 1
          state <- paste(x$stateInfo$genes[which(col == 1)],collapse=", ")
          if (state == "")
            state <- "--"
          cat("State ",i,": ",state,"\n",sep="")
        })
        cat("\n")
      }
      else
      {
        cat(":\n\n")
        cat(" |--<",paste(rep("-",numGenes-1),collapse=""),"|\n",sep="")
        cat(" V ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
        apply(binMatrix,2,function(col)
        {
          cat(" ",col,"   |\n",sep="")
          cat(" | ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
        
        })
        cat(" V ",paste(rep(" ",numGenes-1),collapse=""),"  |\n",sep="")
        cat(" |-->",paste(rep("-",numGenes-1),collapse=""),"|\n\n",sep="")
      }
     }
     else
     {
           # print general information on the attractor
        cat("Attractor ",i," is a complex/loose attractor consisting of ",ncol(attractors[[i]]$involvedStates),
          " state(s) and ",ncol(attractors[[i]]$initialStates), " transition(s)",sep="")
        
        if (activeOnly)
        {
          cat(".\nActive genes in the state transitions: \n")
          initialStates <- t(apply(attractors[[i]]$initialStates,2,function(state)
                                 dec2bin(state,numGenes)))
          nextStates <- t(apply(attractors[[i]]$nextStates,2,function(state)
                                 dec2bin(state,numGenes)))
                                 
          binMatrix <- data.frame(initialStates,nextStates)                                 

          apply(binMatrix,1,function(row)
          {
            state1 <- paste(x$stateInfo$genes[which(row[seq_len(numGenes)] == 1)],collapse=", ")
            
            if (state1 == "")
              state1 <- "--"
            
            state2 <- paste(x$stateInfo$genes[which(row[seq_len(numGenes) + numGenes] == 1)],collapse=", ")
            
            if (state2 == "")
              state2 <- "--"
            
            cat(state1," => ",state2,"\n",sep="")
          })

        }
        else
        { 
          cat(":\n\n")
          initialStates <- apply(attractors[[i]]$initialStates,2,function(state)
                                 paste(dec2bin(state,numGenes),collapse=""))
          nextStates <- apply(attractors[[i]]$nextStates,2,function(state)
                                 paste(dec2bin(state,numGenes),collapse=""))

          binMatrix <- data.frame(initialStates,nextStates)

          apply(binMatrix,1,function(row)
          {         
              cat(row[1]," => ",row[2],"\n",sep="")
          })
        }
     }
     if (!activeOnly)
       cat("\nGenes are encoded in the following order: ",paste(x$stateInfo$genes,collapse=" "),"\n\n",sep="")
  })  
  return(invisible(x))
}
