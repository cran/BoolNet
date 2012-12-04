# Load a network in a specified rule description language
# from file <file>.
# <bodySeparator> is the character that separates targets and factors
# <lowercaseGenes> specifies whether gene names are converted to lower case
loadNetwork <- function(file, bodySeparator=",", lowercaseGenes=FALSE) 
{
  op <- c("!","&", "\\|", "\\(", "\\)")

  func <- readLines(file,-1)[-1]
  
  func <- func[sapply(func, function(str)substr(str,1,1) != "#")]
  
  if (lowercaseGenes)
    func <- tolower(func)
  
  func <- func[nchar(func) > 0]

  # / in a gene name disturbs parsing
  func <- gsub("/","_",func) 

  tmp <-  unname(lapply(func,function(x){strsplit(x,bodySeparator)[[1]]}))
  targets <- sapply(tmp,function(rule)rule[1])
  factors <- sapply(tmp,function(rule)rule[2])
  probabilities <- sapply(tmp,function(rule)
                          {
                            if (length(rule) >= 3)
                              as.numeric(rule[3])
                            else
                              1.0
                          })

  factors.tmp <- lapply(factors,function(x)
  # extract gene names from Boolean expressions
  {
    # replace operators by spaces
    sapply(op,function(y)
    {
      x <<- gsub(y," ",x)
    })
    
    # create a list of involved factors
    tmp <- strsplit(x," ")[[1]]
    tmp <- unique(tmp[tmp != ""])
  })

  # create list of all gene names in both sides of the functions
  genes <- unique(c(targets,unname(unlist(factors.tmp))))
  
  isProbabilistic <- (length(unique(targets)) < length(targets))

  # extract "real" gene names from the list, drop constants
  suppressWarnings(genes <- genes[is.na(as.integer(genes))])

  fixed <- rep(-1,length(genes))
  names(fixed) <- genes

  interactions <- list()

  for(i in 1:length(targets))
  {
    target <- targets[i]
    inputGenes <- factors.tmp[[i]]
    interaction <- list()
    if(suppressWarnings(is.na(as.integer(inputGenes[1]))))
    # the input is not a number
    {
      interaction <- generateInteraction(factors[i], inputGenes, genes) 
    }
    else
    # this is a constant gene
    {
      if (!isProbabilistic)
        fixed[target] <- as.integer(inputGenes)
      interaction <- list(input=0,func=as.integer(inputGenes),expression = inputGenes)
    }
    if (isProbabilistic)
    {
      interaction$probability <- probabilities[i]
      interactions[[target]][[length(interactions[[target]]) + 1]] <- interaction
    }
    else
      interactions[[target]] <- interaction
    
  }

  onlyInputs <- setdiff(genes,targets)
  if(length(onlyInputs) > 0)
  # some genes are only used as inputs, but are not involved in the network
  # -> create dummy input and function
  {
    for(gene in onlyInputs)
    {
      warning(paste("There is no transition function for gene \"",
                     gene,"\"! Assuming an input!",sep=""))
      if (isProbabilistic)
        interactions[[gene]][[1]] = list(list(input = length(interactions)+1,func=c(0,1),
        									  expression = gene))
      else
        interactions[[gene]] = list(input = length(interactions)+1,func=c(0,1),
        							expression = gene)
    }
  }
  
  if (isProbabilistic)
  {
    wrongProb <- sapply(interactions,function(interaction)
                                    abs(1.0-sum(sapply(interaction,function(func)func$probability))) > 0.0001)
    if (any(wrongProb))
      stop(paste("The probabilities of gene(s) ",paste(genes[wrongProb],collapse=", ")," do not sum up to 1!",sep=""))
  }  

  res <- list(interactions = interactions,
              genes = genes,
              fixed = fixed)
  
  if (isProbabilistic)
    class(res) <- c("ProbabilisticBooleanNetwork","BooleanNetworkCollection")
  else
    class(res) <- "BooleanNetwork"
  return(res)
}

