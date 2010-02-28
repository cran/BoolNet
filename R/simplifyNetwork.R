# Simplify a Boolean network by dropping all input genes that are not relevant 
# for the output
simplifyNetwork <- function(network,readableFunctions=FALSE)
{
  stopifnot(inherits(network,"BooleanNetwork") | inherits(network,"ProbabilisticBooleanNetwork")
          | inherits(network,"BooleanNetworkCollection"))
  
  if (inherits(network,"BooleanNetwork"))
  {
    j <- 0
    network$interactions <- lapply(network$interactions,function(interaction)
          {
            if (interaction$input[1] != 0)
            # no constant gene
            {
              j <<- j + 1
        
              table <- allcombn(2,length(interaction$input)) - 1
        
              dropGenes <- apply(table,2,function(gene)
              # drop all genes that have no influence on the results,
              # i.e. the result column is equal for 0 values and 1 values
                  {
                    (identical(interaction$func[gene==1],
                         interaction$func[gene==0]))
                  })
              if (sum(!dropGenes) == 0)
              {
                network$fixed[j] <<- interaction$func[1]
                interaction$input <- 0
                interaction$func <- interaction$func[1]
              }
              else
              if (sum(dropGenes) > 0)
              {
                # update network
                network$fixed[j] <<- -1
                dropFunctionIndices <- unlist(sapply(1:length(dropGenes),function(i)
                      {
                        if(dropGenes[i])
                          which(table[,i] == 0)
                        else
                          NULL
                      }))
                interaction$input <- interaction$input[!dropGenes]
                interaction$func <- interaction$func[-dropFunctionIndices]
                interaction$expression <- 
                  getInteractionString(readableFunctions,
                           interaction$func,
                           network$genes[interaction$input])
              }
            }
            interaction
          })
  }
  else
  {
    network$interactions <- lapply(network$interactions,function(gene)
                                   lapply(gene,function(interaction)
      {
        if (interaction$input[1] != 0)
        # no constant gene
        {

          table <- allcombn(2,length(interaction$input)) - 1
  
          dropGenes <- apply(table,2,function(gene)
          # drop all genes that have no influence on the results,
          # i.e. the result column is equal for 0 values and 1 values
              {
                (identical(interaction$func[gene==1],
                     interaction$func[gene==0]))
              })
          if (sum(!dropGenes) == 0)
          {
            interaction$input <- 0
            interaction$func <- interaction$func[1]
          }
          else
          if (sum(dropGenes) > 0)
          {
            # update network
            dropFunctionIndices <- unlist(sapply(1:length(dropGenes),function(i)
                  {
                    if(dropGenes[i])
                      which(table[,i] == 0)
                    else
                      NULL
                  }))
            interaction$input <- interaction$input[!dropGenes]
            interaction$func <- interaction$func[-dropFunctionIndices]
            interaction$expression <- 
              getInteractionString(readableFunctions,
                       interaction$func,
                       network$genes[interaction$input])
          }
        }
        interaction
      }))
  }
  return(network)
}