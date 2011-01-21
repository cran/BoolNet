# Create LaTeX state tables of all attractors in <attractorInfo>.
# Genes are grouped according to <grouping>.
# An additional title can be supplied in <title>.
# If <plotFixed> is set, fixed variables are included in the plot.
# <onColor> and <offColor> specify the colors of ON/1 and OFF/0 states.
# <file> is the name of the output LaTeX document.
attractorsToLaTeX <- function (attractorInfo, subset, title = "", grouping = list(), plotFixed = TRUE,
        onColor="[gray]{0.9}",offColor="[gray]{0.6}", file="attractors.tex")  
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))
  
  if (missing(subset))
      subset <- 1:length(attractorInfo$attractors)
  else
    if (any(subset > length(attractorInfo$attractors)))
      stop("You specified an attractor index that is greater than the total number of attractors in 'subset'!")

  
  numGenes <- length(attractorInfo$stateInfo$genes)
  
  # escape "_" in LaTeX
  genes = gsub("_", "\\_", attractorInfo$stateInfo$genes)
  
  # determine list of genes to be plotted
  whichFixed <- which(attractorInfo$stateInfo$fixedGenes != -1)
  if (plotFixed | (length(whichFixed) == 0))
    plotIndices <- 1:numGenes
  else
    plotIndices <- (1:numGenes)[-whichFixed]
  
      # convert decimal state numbers to binary state matrices (one for each attractor)
  binMatrices <- lapply(attractorInfo$attractors,function(attractor)
          {
            res <- matrix(apply(attractor$involvedStates,2,function(state)
              dec2bin(state,numGenes)[plotIndices]),nrow=length(plotIndices))
          })

  # count the numbers of attractors with equal lengths
  attractorLengths <- sapply(attractorInfo$attractors,function(attractor)
                             {
                                if (is.null(attractor$initialStates))
                                # simple attractor
                                  ncol(attractor$involvedStates)
                                else
                                # complex attractor => extra treatment
                                  -1
                             })  
  lengthTable <- table(attractorLengths)
  lengthTable <- lengthTable[as.integer(names(lengthTable)) != -1]
  
  # Open output file, and print header
  sink(file)
  cat("% Please include packages tabularx and colortbl in your master document:\n",
      "% \\usepackage{tabularx,colortbl}\n\n\n",sep="")
      
  res <- lapply(1:length(lengthTable),function(i)
  # accumulate all attractors with equal length in one matrix and plot them
  {
    len <- as.integer(names(lengthTable)[i])
    
     if (length(intersect(which(attractorLengths == len),subset)) > 0)
     {
      cnt <- lengthTable[i]

      # build accumulated matrix     
      totalMatrix <- c()
      for (mat in binMatrices[intersect(which(attractorLengths == len),subset)])
      {
        totalMatrix <- cbind(totalMatrix,mat)
      }
      rownames(totalMatrix) <- attractorInfo$stateInfo$genes[plotIndices]
      colnames(totalMatrix) <- sapply(intersect(which(attractorLengths == len),subset),function(i)paste("Attr",i,".",1:len,sep=""))
    
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
               "\\caption{",
               title, "Attractors with ",len," state(s)}\n",
               "\\begin{tabularx}{\\linewidth}{l", 
             paste(rep(paste(rep(">{\\centering\\arraybackslash}X", 
                  len), collapse = " "),length(intersect(which(attractorLengths == len),subset))),collapse="|"), 
          "}\\hline\n",
               sep="")
    
      cat("\t&\t",paste(paste("\\multicolumn{",len,"}{c}{Attr. ",intersect(which(attractorLengths == len),subset),"}",
            sep=""),collapse="\t&\t"),"\\\\\n")    
    
      # output active and inactive states
      for(j in 1:nrow(totalMatrix))
      {
        separator <- which(separationPositions==j)
        if (length(separator) != 0)
        {
          cat("\\hline \\multicolumn{",ncol(totalMatrix) + 1,"}{c}{",grouping$class[separator],"}\\\\ \\hline \n",sep="")
        }
        cat("\\textbf{",rownames(totalMatrix)[j],"}\t&\t",sep="")
        for(i in 1:ncol(totalMatrix))
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
    
          # output frequency of attractor (basin size / number of states)
      freq <- round(sapply(attractorInfo$attractors[intersect(which(attractorLengths == len),subset)],
          function(attractor)attractor$basinSize/ncol(attractorInfo$stateInfo$table)) * 100,2)

      if (!isTRUE(all(is.na(freq))))
      {
        cat("\\hline Freq.\t&\t",paste(paste("\\multicolumn{",len,"}{c}{",freq,"\\%}",
              sep=""),collapse="\t&\t"),"\\\\\n")
      }

      cat("\\hline\\end{tabularx}\n\\end{center}\n",
          "\\end{table}\n\n",sep="")

      totalMatrix
    }
  })
  
  # return a list of accumulated matrices
  sink()
  names(res) <- names(lengthTable)
  return(res)
}
