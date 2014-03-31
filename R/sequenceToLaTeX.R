
sequenceToLaTeX <- function(network, startState, includeAttractorStates = c("all","first","none"),
                            sequence, title="", grouping = list(), plotFixed = TRUE,
                            onColor="[gray]{0.9}",offColor="[gray]{0.6}", file="sequence.tex")
{
  if (!missing(network))
  {
    stopifnot(inherits(network,"BooleanNetwork"))
    if (missing(startState) || !missing(sequence))
      stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
    
    sequence <- getPathToAttractor(network = network,
                                   state = startState, 
                                   includeAttractorStates = includeAttractorStates)

    numGenes <- length(network$genes)                                   
    whichFixed <- which(network$fixed != -1)
    if (plotFixed | (length(whichFixed) == 0))
      plotIndices <- seq_len(numGenes)
    else
      plotIndices <- seq_len(numGenes)[-whichFixed]
      
    sequence <- sequence[,plotIndices]                                   
  }
  else
  {
    if (missing(sequence) || !missing(startState))
        stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
  }
  
  # escape "_" in LaTeX
  genes = gsub("_", "\\_", colnames(sequence))
  
  # determine list of genes to be plotted
  
  # Open output file, and print header
  sink(file)
  cat("% Please include packages tabularx and colortbl in your master document:\n",
      "% \\usepackage{tabularx,colortbl}\n\n\n",sep="")
      
  totalMatrix <- t(sequence)
  colnames(totalMatrix) <- seq_len(ncol(totalMatrix))
  
  if(length(grouping)>0)
  {
     # reorder genes according to the supplied groups
    totalMatrix <- totalMatrix[unlist(grouping$index),]
    separationPositions <- c(1,cumsum(sapply(grouping$index,length)+1))
  }
  else
    separationPositions <- c()

  # output table header
  cat("\\begin{table}[ht]\n",
       "\\begin{center}\n",
       "\\caption{",title,"}\n",
       "\\begin{tabularx}{\\linewidth}{l", 
     paste(rep(paste(rep(">{\\centering\\arraybackslash}X", 
          ncol(totalMatrix), collapse = " "))),collapse=""), 
  "}\\hline\n", sep="")
       
  cat("\\textbf{Time}\t&\t",paste(seq_len(ncol(totalMatrix)),collapse="\t&\t"),"\\\\")
  
   if(length(grouping) == 0)
     cat("\\hline\n")
   else
     cat("\n")  

  # output active and inactive states
  for(j in seq_len(nrow(totalMatrix)))
  {
    separator <- which(separationPositions==j)
    if (length(separator) != 0)
    {
      cat("\\hline \\multicolumn{",ncol(totalMatrix) + 1,"}{c}{",grouping$class[separator],"}\\\\ \\hline \n",sep="")
    }
    cat("\\textbf{",rownames(totalMatrix)[j],"}\t&\t",sep="")
    for(i in seq_len(ncol(totalMatrix)))
    {
      if(totalMatrix[j,i] == 1)
        cat("\\cellcolor",onColor,"1",sep="")
      else
        cat("\\cellcolor",offColor,"0",sep="")
      if (i < ncol(totalMatrix))
        cat("\t&\t")
    }
    cat("\\\\\n")
  }

  cat("\\hline\\end{tabularx}\n\\end{center}\n",
      "\\end{table}\n\n",sep="")

  # return a list of accumulated matrices
  sink()
  return(totalMatrix)
}

