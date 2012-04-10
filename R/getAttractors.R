# Identify attractors in a Boolean network.
# <network> is a BooleanNetwork structure specifying the network.
# <method> specifies what kind of search is conducted: "exhaustive" performs a full search over all states,
# "random" generates <startStates> random initial states, and "chosen" uses the states supplied in <startStates>.
# <genesON> and <genesOFF> are lists of genes to be set to ON/1 or OFF/0 respectively.
# If <canonical> is true, states in the attractors are reordered such that the "smallest" state is the first
getAttractors <- function (network, type=c("synchronous","asynchronous"), 
         method=c("exhaustive","random","chosen"), startStates=list(),
         genesON = c(), genesOFF = c(), canonical=TRUE,
         randomChainLength = 10000, avoidSelfLoops = TRUE, 
         geneProbabilities = NULL, 
         returnTable=TRUE) 
{
  stopifnot(inherits(network,"BooleanNetwork"))
  
  if (match.arg(type) == "asynchronous" & length(method) == 1 & match.arg(method) == "exhaustive")
    stop("Asynchronous attractor search cannot be performed in exhaustive search mode!")
  
  if (match.arg(type) == "asynchronous" & length(geneProbabilities) > 0 )
  {
    if (length(geneProbabilities) != length(network$genes))
      stop("Please supply exactly one probability for each gene!")
    if (abs(1.0 - sum(geneProbabilities)) > 0.0001)
      stop("The supplied gene probabilities do not sum up to 1!")
  }

  if (length(method) == 3)
  # if no method is supplied, infer method from the type of <startStates>
  {
    if (match.arg(type) == "asynchronous" & length(startStates) == 0)
    {
      startStates <- max(round(2 ^ sum(network$fixed == -1) / 20), 5)
      method = "random"
    }
    else
    if (is.numeric(startStates))
    {
      if (length(startStates) > 1)
        stop("Please supply either the number of start states or a list of start states in startStates!")
      else
        method <- "random"
    }
    else
    if (is.list(startStates) & (length(startStates) > 0))
      method = "chosen"
    else
      method = "exhaustive"
  }

      
  if ((length(network$genes) > 29) & (match.arg(method) == "exhaustive") & (match.arg(type) == "synchronous"))
    stop(paste("An exhaustive attractor search is not supported for networks with more than 29 genes!",
          "Please use a heuristic method or provide start states!"))
  
  # fix genes according to genesON and genesOFF
  if (length(genesON) > 0) 
  {
    network <- fixGenes(network,genesON,1)
  }
  if (length(genesOFF) > 0) 
  {
    network <- fixGenes(network,genesOFF,0)
  }
  
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
          if (match.arg(type) == "synchronous")
          {
            list()
            warning(paste("The number of random states is set larger than the total",
                          "number of states. Performing an exhaustive search!"))
          }
          else
          {
            warning(paste("The number of random states is set larger than the total ",
                        "number of states! The maximum number of different states is ",2 ^ length(nonFixedPositions)),"!",sep="")
            startStates = 2 ^ length(nonFixedPositions)            
          }     
        }
        if (startStates <= (2 ^ length(nonFixedPositions)))
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
            mat <- unique(rbind(mat,vec))
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
          fixedPositions <- which(network$fixed != -1)
  
          fixedGenes <- which(network$fixed != -1)
          statesValid <- sapply(startStates,function(state)
                                {
                                  isTRUE(all(state[fixedGenes] == network$fixed[fixedGenes]))
                                })
      startStates <- startStates[statesValid]
          if (!isTRUE(all(statesValid)))
            warning("Some of the supplied start states did not match the restrictions of the fixed genes and were removed!")    
            
        startStates
       }
  )
  
  specialInitialization <- NULL
  
  convertedStartStates <- NULL

  if (length(startStates) > 0)
    convertedStartStates <- sapply(startStates,function(x)bin2dec(x,length(network$genes)))

  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputGenes <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$input)))
  inputGenePositions <- as.integer(cumsum(c(0,sapply(network$interactions,
           function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  transitionFunctions <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$func)))
  transitionFunctionPositions <- as.integer(cumsum(c(0,sapply(network$interactions,
              function(interaction)length(interaction$func)))))

  networkType <- switch(match.arg(type,c("synchronous","asynchronous")),
    synchronous = 0,
    asynchronous = 1)
  
  on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))
  # Call the C code
  result <- .Call("getAttractors_R",inputGenes,inputGenePositions,
        transitionFunctions,transitionFunctionPositions,
        as.integer(network$fixed),
        as.integer(convertedStartStates),
        as.integer(networkType),
        as.double(geneProbabilities),
        as.integer(randomChainLength),
        as.integer(avoidSelfLoops),
        as.integer(returnTable),
        PACKAGE="BoolNet")
  
  if (is.null(result))
    stop("An error occurred in external C code!")
  
  if (length(result$attractors) == 0)
    stop("getAttractors() was not able to identify any attractors! Please check the supplied parameters and restart!")
  
  if (length(network$genes) %% 32 == 0)
    numElementsPerEntry <- as.integer(length(network$genes) / 32)
  else
    numElementsPerEntry <- as.integer(length(network$genes) / 32  + 1)
  
  if (!is.null(result$stateInfo))
  {
    result$stateInfo$table <- matrix(result$stateInfo$table,nrow=numElementsPerEntry)
  
    if (!is.null(result$stateInfo$initialStates))
      result$stateInfo$initialStates <- matrix(result$stateInfo$initialStates,nrow=numElementsPerEntry)
  }
  
  for (i in 1:length(result$attractors))
  {
    result$attractors[[i]]$involvedStates <- matrix(result$attractors[[i]]$involvedStates,nrow=numElementsPerEntry)
    if (canonical)
    # reorder states
      result$attractors[[i]]$involvedStates <- canonicalStateOrder(result$attractors[[i]]$involvedStates)
      
    if (!is.null(result$attractors[[i]]$initialStates))
      result$attractors[[i]]$initialStates <- matrix(result$attractors[[i]]$initialStates,nrow=numElementsPerEntry)
    if (!is.null(result$attractors[[i]]$nextStates))
      result$attractors[[i]]$nextStates <- matrix(result$attractors[[i]]$nextStates,nrow=numElementsPerEntry)
      
    if (result$attractors[[i]]$basinSize == 0)
      result$attractors[[i]]$basinSize <- NA
  }
  
  # order attractors according to their lengths
  attractorLengths <- sapply(result$attractors,function(attractor)ncol(attractor$involvedStates))  
  reordering <- order(attractorLengths)
  result$attractors <- result$attractors[reordering]
  
  if (!is.null(result$stateInfo))
  {
    inverseOrder <- sapply(1:length(reordering),function(x)which(reordering == x))    
    result$stateInfo$attractorAssignment <- inverseOrder[result$stateInfo$attractorAssignment]
  }
  
  # extend the resulting structure by additional information, and assign a class      
  result$stateInfo$genes <- network$genes
  result$stateInfo$fixedGenes <- network$fixed

  if (!is.null(result$stateInfo$table))
    class(result$stateInfo) <- "BooleanStateInfo"
  class(result) <- "AttractorInfo"
  return(result)
}
