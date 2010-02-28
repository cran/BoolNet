# Plot state tables of all attractors in <attractorInfo>.
# Genes are grouped according to <grouping>.
# An additional title can be supplied in <title>.
# If <plotFixed> is set, fixed variables are included in the plot.
# <onColor> and <offColor> specify the colors of ON/1 and OFF/0 states.
plotAttractors <- function (attractorInfo, subset, title = "", mode=c("table","graph"),
                            grouping = list(), plotFixed = TRUE, onColor="green",offColor="red",
                            layout=layout.circle, drawLabels=TRUE,...) 
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))
  numGenes <- length(attractorInfo$stateInfo$genes)
  
  if (missing(subset))
      subset <- 1:length(attractorInfo$attractors)
  else
    if (any(subset > length(attractorInfo$attractors)))
      stop("You specified an attractor index that is greater than the total number of attractors in 'subset'!")
  
  switch(match.arg(mode,c("table","graph")),
  table =
  {
  
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
  
    res <- lapply(1:length(lengthTable),function(i)
    # accumulate all attractors with equal length in one matrix and plot them
    {
      len <- as.integer(names(lengthTable)[i])
      if (length(intersect(which(attractorLengths == len),subset)) > 0)
      {
        cnt <- lengthTable[i]

        # initialize with empty plot
        plot(c(),c(),xlim=c(0,len*cnt),ylim=c(-2,length(plotIndices)+1),xlab="",ylab="",
             axes=FALSE,main=paste(title, "Attractors with ",len," state(s)",sep=""))
    
        # build accumulated matrix     
        totalMatrix <- c()
        for (mat in binMatrices[intersect(which(attractorLengths == len),subset)])
        {
          totalMatrix <- cbind(totalMatrix,mat)
        }
        rownames(totalMatrix) <- attractorInfo$stateInfo$genes[plotIndices]
        colnames(totalMatrix) <- sapply(intersect(which(attractorLengths == len),subset),function(i)paste("Attr",i,".",1:len,sep=""))
    
        if(length(grouping)>0)
           # reorder genes according to the supplied groups
          totalMatrix = totalMatrix[unlist(grouping$index),]
    
        par(yaxt='s',las=2)
            axis(2,(1:length(plotIndices))-0.5,rownames(totalMatrix))

        # plot active and inactive states
        for(i in 1:ncol(totalMatrix))
          for(j in 1:nrow(totalMatrix))
          {
            if(totalMatrix[j,i] == 1)
              rect(i-1,j-1,i,j,col=onColor,border="gold")
            else
              rect(i-1,j-1,i,j,col=offColor,border="gold")
          }
    
        # draw vertical separators between attractors  
            sep = seq(0,ncol(totalMatrix),by=len)    
            abline(v = sep[-1],col="white",lwd=3)
            
        # output frequency of attractor (basin size / number of states)
        freq <- sapply(attractorInfo$attractors[intersect(which(attractorLengths == len),subset)],
            function(attractor)attractor$basinSize/length(attractorInfo$stateInfo$table)) * 100

        if (!isTRUE(all(is.na(freq))))
        {
          text(sep[1:(length(sep)-1)] + len/2, rep(0.4+nrow(totalMatrix),ncol(totalMatrix)),
            paste(freq,"%",sep=""),cex=.75,font=3)
        }
      
        if(length(grouping)>0)
        # draw separators between groups, and print group names
        {
          sepPos = cumsum(sapply(grouping$index,length))
          abline(h=sepPos[-length(sepPos)],col="black",lwd=3)
          text(ncol(totalMatrix)/2,sepPos-0.5,grouping$class,cex=0.9)
              }

          legend(x="bottomright",pch=c(15,15),col=c(onColor,offColor),legend = c("active","inactive"),cex=0.7,horiz=T)
        totalMatrix
       }
    })
  
    # return a list of accumulated matrices
    names(res) <- names(lengthTable)
    return(res)
  },
  
  graph = 
  {
    if (!require(igraph))
      stop("Please install the igraph package before using the \"graph\" mode of this function!")
    
     args <- list(...)
      
    if (is.null(args$vertex.size))
      args$vertex.size <- 2
    
    if (is.null(args$edge.arrow.mode))
      args$edge.arrow.mode <- 0
    
    if (is.null(args$vertex.label.cex))
      args$vertex.label.cex <- 0.75
    
    if (is.null(args$vertex.label.dist))
      args$vertex.label.dist <- 1  
      
     if (is.null(args$vertex.color))
      args$vertex.color <- "grey"
    
    if (is.null(args$edge.arrow.size))
      args$edge.arrow.size <- 0.5  
      
    lapply(attractorInfo$attractors[subset],function(attractor)
    {
      nodes <- data.frame(apply(attractor$involvedStates,2,function(state)
                          paste(dec2bin(state,numGenes),collapse="")))
    
      if (!is.null(attractor$initialStates))
      # asynchronous complex attractor
      {
        initialStates <- apply(attractor$initialStates,2,function(state)
                                paste(dec2bin(state,numGenes),collapse=""))
        nextStates <- apply(attractor$nextStates,2,function(state)
                                paste(dec2bin(state,numGenes),collapse=""))
      }
      else
      {
        initialStates <- apply(attractor$involvedStates,2,function(state)
                                paste(dec2bin(state,numGenes),collapse=""))
        if (length(initialStates) == 1)
        # steady state
          nextStates <- initialStates
        else
        # synchronous attractor with more than one state
          nextStates <- initialStates[c(2:length(initialStates),1)]                                
      }
      edgeMatrix <- data.frame(initialStates,nextStates)
      graph <- graph.data.frame(edgeMatrix,vertices=nodes,directed=TRUE)
      
      if (drawLabels)
        labels <- nodes[,1]
      else
        labels <- NA
      
       plot(graph,layout=layout,vertex.label=labels,vertex.label.cex=args$vertex.label.cex,
         vertex.size=args$vertex.size,vertex.color=args$vertex.color,
         vertex.label.dist=args$vertex.label.dist,
         edge.arrow.size=args$edge.arrow.size,
         main=title,...)
      graph
    })
  })
}
