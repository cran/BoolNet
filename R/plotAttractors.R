# Plot state tables of all attractors in <attractorInfo>.
# Genes are grouped according to <grouping>.
# An additional title can be supplied in <title>.
# If <plotFixed> is set, fixed variables are included in the plot.
# <onColor> and <offColor> specify the colors of ON/1 and OFF/0 states.
# TODO: remove old definition from comments
#plotAttractors <- function (attractorInfo, subset, title = "", mode=c("table","graph"),
#                            grouping = list(), plotFixed = TRUE, onColor="green",offColor="red",
#                            layout=layout.circle, drawLabels=TRUE, drawLegend=TRUE, ask=TRUE, 
#                            reverse=FALSE, ...) 
#{}
plotAttractors <- function (attractorInfo, subset, 
                            title = "", 
                            mode = c("table", "graph"),
                            grouping = list(), 
                            plotFixed = TRUE,
                            onColor = "#4daf4a", offColor = "#e41a1c",
                            layout = layout.circle,
                            drawLabels = TRUE, drawLegend = TRUE, ask = TRUE, reverse = FALSE, 
                            # new parameters
                            borderColor = "black",
                            eps = 0.1, allInOnePlot = FALSE, ...)
{
  stopifnot(inherits(attractorInfo, "AttractorInfo") || inherits(attractorInfo,
                                                                 "SymbolicSimulation"))
  if (inherits(attractorInfo, "AttractorInfo")) {
    numGenes <- length(attractorInfo$stateInfo$genes)
    geneNames <- attractorInfo$stateInfo$genes
  }
  else {
    numGenes <- ncol(attractorInfo$attractors[[1]])
    geneNames <- colnames(attractorInfo$attractors[[1]])
  }
  if (missing(subset))
    subset <- seq_along(attractorInfo$attractors)
  else if (any(subset > length(attractorInfo$attractors)))
    stop("You specified an attractor index that is greater than the total number of attractors in 'subset'!")
  res <- switch(match.arg(mode, c("table", "graph")), table = {
    whichFixed <- which(attractorInfo$stateInfo$fixedGenes != -1)
    if (plotFixed | (length(whichFixed) == 0)) {
      plotIndices <- seq_len(numGenes)
    } else {
      plotIndices <- seq_len(numGenes)[-whichFixed]
    }
    if (inherits(attractorInfo, "AttractorInfo")) {
      binMatrices <- lapply(attractorInfo$attractors, function(attractor) {
        res <- matrix(apply(attractor$involvedStates, 2, function(state) dec2bin(state, numGenes)[plotIndices]), nrow = length(plotIndices))
      })
      attractorLengths <- sapply(attractorInfo$attractors, function(attractor) {
        ifelse (is.null(attractor$initialStates), ncol(attractor$involvedStates), -1)
      })
    } else {
      binMatrices <- lapply(attractorInfo$attractors, t)
      attractorLengths <- sapply(binMatrices, ncol)
    }
    lengthTable <- table(attractorLengths)
    lengthTable <- lengthTable[as.integer(names(lengthTable)) != -1]
    oldAsk <- par("ask")
    oldmfrow <- NULL
    if (allInOnePlot) {
      oldmfrow <- par()$mfrow
      #par(mfrow=c(1,length(table(sapply(sapply(attractorInfo$attractors, "[[", "involvedStates", simplify=F), ncol)))))
      par(mfrow=c(1,length(lengthTable)))
    }
    res <- lapply(seq_along(lengthTable), function(i) {
      len <- as.integer(names(lengthTable)[i])
      attractorIndices <- intersect(which(attractorLengths == len), subset)
      if (length(attractorIndices) > 0) {
        cnt <- length(attractorIndices)
        totalMatrix <- c()
        for (mat in binMatrices[attractorIndices]) {
          totalMatrix <- cbind(totalMatrix, mat)
        }
        rownames(totalMatrix) <- geneNames[plotIndices]
        colnames(totalMatrix) <- sapply(attractorIndices, function(i) paste("Attr", i, ".", seq_len(len), sep = ""))
        if (length(grouping) > 0) totalMatrix <- totalMatrix[unlist(grouping$index), , drop = FALSE]
        par(ask = ask && i > 1 && dev.interactive())
        if (reverse) {
          plot(c(), c(),
               xlim = c(0, len *cnt), ylim = c(-2, nrow(totalMatrix) + 1),
               xlab = "", ylab = "", axes = FALSE,
               main = paste(title, "\nAttractors with ", len, " state(s)", sep = ""))
        } else {
          plot(c(), c(),
               xlim = c(0, len * cnt), ylim = c(nrow(totalMatrix) + 1, -2),
               xlab = "", ylab = "", axes = FALSE,
               main = paste(title, "\nAttractors with ", len," state(s)", sep = ""))
        }
        #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
        axis(2, seq_len(nrow(totalMatrix)) - 0.5, rownames(totalMatrix), yaxt = "s", las = 2)
        xStart <- 1
        unitFactor <- (len - (2 * eps))/len
        startX <- eps
        for (i in seq_len(ncol(totalMatrix))) {
          for (j in seq_len(nrow(totalMatrix))) {
            rectCol <- ifelse(totalMatrix[j, i], onColor, offColor)
            rect(startX, j - 1, startX + unitFactor, j, col = rectCol, border = borderColor, lwd = 2)
          }
          startX <- startX + unitFactor
          if (i %% len == 0) {
            startX <- startX + 2*eps
          }
        }
        sep = seq(0, ncol(totalMatrix), by = len)
        #abline(v = sep[-1], col = par()$bg, lwd = 3)
        if (inherits(attractorInfo, "AttractorInfo")) {
          if (is.null(attractorInfo$stateInfo$table)) {
            freq <- rep(NA, length(attractorIndices))
          }
          else {
            freq <- round(100 * sapply(attractorInfo$attractors[attractorIndices], function(attractor) {
              attractor$basinSize/ncol(attractorInfo$stateInfo$table)}), 2)
          }
        } else {
          if (!is.null(attractorInfo$graph)) {
            freq <- round(sapply(attractorIndices, function(i) sum(attractorInfo$graph$attractorAssignment == i)/nrow(attractorInfo$graph)) * 100, 2)
          } else {
            freq <- rep(NA, length(attractorIndices))
          }
        }
        if (!isTRUE(all(is.na(freq)))) {
          text(sep[seq_along(sep) - 1] + len/2, ifelse(reverse, -nrow(totalMatrix)-1.5, 1) + rep(1 + nrow(totalMatrix), ncol(totalMatrix)), paste(freq, "%", sep = ""), cex = 0.75, font = 3)
        }
        if (length(grouping) > 0) {
          if(!is.null(grouping$class)) {
            sepPos = cumsum(sapply(grouping$index, length))
            abline(h = sepPos[-length(sepPos)], col = "black", lwd = 3)
            text(ncol(totalMatrix)/2, sepPos - 0.5, grouping$class, cex = 0.9)
          }
        }
        #if (drawLegend) legend(x = "bottomright", pch = c(15, 15), col = c(onColor, offColor), legend = c("active", "inactive"), cex = 0.7, horiz = T, xpd = T)
        if (drawLegend) legend(x = 0, y = ifelse(reverse, -nrow(totalMatrix)-1, 2) + nrow(totalMatrix), pch = c(15, 15), col = c(onColor, offColor), legend = c("active", "inactive"), cex = 0.7, horiz = T, xpd=T)
        totalMatrix
      }
    })
    par(ask = oldAsk)
    if(!is.null(oldmfrow)) par(mfrow = oldmfrow)
    names(res) <- names(lengthTable)
    return(res)
    
  }, graph = {
    args <- list(...)
    if (is.null(args$vertex.size)) args$vertex.size <- 2
    if (is.null(args$edge.arrow.mode)) args$edge.arrow.mode <- 0
    if (is.null(args$vertex.label.cex)) args$vertex.label.cex <- 0.75
    if (is.null(args$vertex.label.dist)) args$vertex.label.dist <- 1
    if (is.null(args$vertex.color)) args$vertex.color <- "grey"
    if (is.null(args$edge.arrow.size)) args$edge.arrow.size <- 0.5
    lapply(attractorInfo$attractors[subset], function(attractor) {
      if (inherits(attractorInfo, "AttractorInfo")) {
        nodes <- data.frame(apply(attractor$involvedStates, 2, function(state) paste(dec2bin(state, numGenes), collapse = "")), stringsAsFactors = FALSE)
        if (!is.null(attractor$initialStates)) {
          initialStates <- apply(attractor$initialStates, 2, function(state) paste(dec2bin(state, numGenes), collapse = ""))
          nextStates <- apply(attractor$nextStates, 2, function(state) paste(dec2bin(state, numGenes), collapse = ""))
        } else {
          initialStates <- apply(attractor$involvedStates, 2, function(state) paste(dec2bin(state, numGenes), collapse = ""))
          if (length(initialStates) == 1) nextStates <- initialStates else nextStates <- initialStates[c(2:length(initialStates), 1)]
        }
      } else {
        initialStates <- apply(attractor, 1, paste, collapse = "")
        nextStates <- apply(attractor[c(2:nrow(attractor), 1), , drop = FALSE], 1, paste, collapse = "")
        nodes <- data.frame(unique(c(initialStates, nextStates)), stringsAsFactors = FALSE)
      }
      edgeMatrix <- data.frame(initialStates, nextStates)
      graph <- graph.data.frame(edgeMatrix, vertices = nodes, directed = TRUE)
      if (drawLabels) labels <- nodes[, 1] else labels <- NA
      args$layout <- layout
      args$x <- graph
      args$vertex.label <- labels
      args$main <- title
      do.call("plot", args)
      graph
    })
  })
}
