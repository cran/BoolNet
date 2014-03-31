
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
  
  bin = .C("dec2bin",as.integer(bin),as.integer(dec),as.integer(len),NAOK=TRUE)[[1]]
}

# Generate a list of all assignments of n variables with N possible values
allcombn <- function(N,n)
{
  rownum = N^n
  sapply(n:1,function(i)
    {
            rep(seq_len(N),each=N^(i-1),len=rownum)
          })
}

# Get DNF representation of a truth table <truthTable>
# using the gene names in <genes>. 
# If <mode> is "canonical", build a canonical DNF.
# If <mode> is "short", join terms to reduce the DNF
getDNF <- function(truthTable, genes, mode = c("short","canonical"))
{
  if (mode[1] == TRUE)
    mode <- (if (length(genes) <= 12) "short" else "canonical")
    
  mode <- match.arg(mode, c("short","canonical"))
  # check for constant functions
  if (isTRUE(all.equal(truthTable,rep(0,length(truthTable)))))
    return("0")
  else
  if (isTRUE(all.equal(truthTable,rep(1,length(truthTable)))))
    return("1")
  
  # generate truth table
  entries <- allcombn(2,length(genes)) - 1
  colnames(entries) <- genes

 
  if (mode == "short")
  {
    # heuristic minimization
    
    # the 1 terms that need to be covered
    uncoveredEntries <- which(truthTable == 1)

    # current conjunction list
    conjunctions <- list()  
    
    while (length(uncoveredEntries) > 0)
    {
      # take an uncovered entry and expand it
      currentEntry <- entries[uncoveredEntries[1],]
          
      for (gene in genes)
      # test for each gene whether it can be eliminated from the term
      {
        geneIdx <- which(names(currentEntry) == gene)
        candidate <- currentEntry[-geneIdx]
        condition <- rep(TRUE,length(truthTable))
        for (i in seq_along(candidate))
        {
          condition <- condition & (entries[,names(candidate)[i]] == candidate[i])
        }
        
        if (length(unique(truthTable[condition])) == 1)
          # eliminate gene
          currentEntry <- currentEntry[-geneIdx]
      }
      
      # determine which truth table result entries are now covered
      eliminatedEntries <- rep(TRUE,length(truthTable))
      for (i in seq_along(currentEntry))
      {
        eliminatedEntries <- eliminatedEntries & 
                             (entries[,names(currentEntry)[i]] == currentEntry[i])
      }
      uncoveredEntries <- setdiff(uncoveredEntries, which(eliminatedEntries))
      
      # remember conjunction
      conjunctions <- c(conjunctions, list(currentEntry))
    }
    return(paste(paste("(",sapply(conjunctions, function(conj)
    {
      paste(mapply(function(gene, val)
      {
        if (val == 1)
          return(gene)
        else
          return(paste("!",gene,sep=""))
      }, names(conj), conj), collapse=" & ")
    }), ")", sep=""), collapse=" | "))
  }
  else
  {
    # canonical DNF
    conjunctions <- apply(entries[truthTable==1,,drop=FALSE],1,function(conj)
    {
      paste("(",paste(sapply(seq_along(conj),function(lit)
      {
        if (conj[lit])
          genes[lit]
        else
          paste("!",genes[lit],sep="")
      }),collapse=" & "),")",sep="")
    })
    return(paste(conjunctions[conjunctions != ""],collapse = " | "))
  }
}

# Retrieves a string representation of an interaction function by either
# building a DNF (if <readableFunction> is false)
# or returning an unspecific function description 
# (if <readableFunction> is true or "canonical" or "short").
# <truthTable> contains the result column of the interaction's truth table, and
# <genes> contains the names of the involved genes.
getInteractionString <- function(readableFunctions,truthTable,genes)
{
  if (readableFunctions != FALSE)
    getDNF(truthTable,genes, readableFunctions)
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
  stateHash <- seq_len(length(inputStates))
  names(stateHash) <- inputStates
  
  
  # encode output states
  outputStates <- apply(attractorInfo$stateInfo$table,2,function(state)
        {
          bin <- paste(dec2bin(state,length(attractorInfo$stateInfo$genes)),collapse="")
        })
  return(list(vertices=inputStates,edges=cbind(seq_along(inputStates),
      sapply(outputStates,function(state)stateHash[state]))))
}

# Reorders a matrix of states <stateMatrix> (with each column being one state)
# in a canonical way such that the "smallest" state is the first.
# This makes attractor representations unique.
canonicalStateOrder <- function(stateMatrix)
{
  smallestIndex <- -1
  smallestVal <- rep(Inf,nrow(stateMatrix))
  for (i in seq_len(ncol(stateMatrix)))
  # iterate over states
  {
    for (j in seq_len(nrow(stateMatrix)))
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
           stateMatrix[,seq_len(smallestIndex-1),drop=FALSE]))
  else
    return(stateMatrix)
}

# Find a child node named <name>
# of <node>, or return NULL if such a node
# does not exist.
# If <throwError>, throw an error if the node does not exist.
xmlFindNode <- function(node,name,throwError=FALSE)
{
    indices <- which(xmlSApply(node, xmlName) ==  name)
    if (length(indices) == 0)
    {
      if (throwError)
        stop(paste("Node \"",name,"\" is required, but missing!", sep=""))
      else
        return(NULL)
    }
    return(xmlChildren(node)[indices])
}

# Remove leading and trailing whitespace characters from <string>
trim <- function(string)
{
  string <- gsub("^[ ]+", "", string)
  string <- gsub("[ ]+$", "", string)
  return(string)
}

# Generate an interaction by parsing <expressionString>
# and building the corresponding truth table.
# Here, <inputGenes> are the names of the inputs, 
# and <allGenes> is a list of all genes in the network.
# Returns an interaction as used in the BooleanNetwork class.
generateInteraction <- function(expressionString, inputGenes, allGenes)
{

  # strip leading and trailing whitespace characters
  expressionString <- trim(expressionString)
  
  if (length(inputGenes) > 0)
  {
    inputIndices <- match(inputGenes, allGenes)

    names(inputIndices) <- inputGenes
    inputIndices <- sort(inputIndices)
  
    truthTable <- as.matrix(allcombn(2,length(inputIndices)) - 1)
  
    for (gene in inputGenes)
    {
      assign(gene, truthTable[,which(names(inputIndices) == gene)])
    }
  }
  else
    inputIndices <- 0
    
  tryCatch(
    interaction <- list(input = unname(inputIndices),
                     func = as.numeric(eval(parse(text=expressionString))),
                     expression = expressionString),
    error = function(err)
            {
              stop(paste("An error was detected while parsing an expression: ",
                          err$message, "\nExpression: \"",
                         expressionString,"\"",
                         sep=""))
            })
  if (length(interaction$func) == 1 && length(inputGenes) > 0)
  {
    warning("There seems to be a constant function that has specified inputs!")
    interaction$func <- rep(interaction$func, nrow(truthTable))
  }                         
  return(interaction)                         
}

