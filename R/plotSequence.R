plotSequence <- function (network, startState, includeAttractorStates = c("all", "first", "none"),
                          sequence, title = "", mode = c("table", "graph"),
                          plotFixed = TRUE, grouping = list(),
                          onColor = "#4daf4a", offColor = "#e41a1c",
                          layout, drawLabels = TRUE, drawLegend = TRUE,
                          highlightAttractor = TRUE, reverse = FALSE,
                          ### new parameters
                          borderColor = "black", eps=0.1,
                          attractor.sep.lwd = 2, attractor.sep.col = "blue", ...)
{
  if (!missing(network)) {
    stopifnot(inherits(network, "BooleanNetwork") || inherits(network,
                                                              "SymbolicBooleanNetwork"))
    if (missing(startState) || !missing(sequence))
      stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
    sequence <- getPathToAttractor(network = network, state = startState, includeAttractorStates = includeAttractorStates)
    numGenes <- length(network$genes)
    if (match.arg(mode, c("table", "graph")) == "table") {
      whichFixed <- which(network$fixed != -1)
      if (plotFixed | (length(whichFixed) == 0))
        plotIndices <- seq_len(numGenes)
      else plotIndices <- seq_len(numGenes)[-whichFixed]
      attractor <- attributes(sequence)$attractor
      sequence <- sequence[, plotIndices, drop = FALSE]
      attributes(sequence)$attractor <- attractor
    }
  }
  else {
    if (missing(sequence) || !missing(startState))
      stop("Either \"network\" and \"startState\" or \"sequence\" must be provided!")
  }
  switch(match.arg(mode, c("table", "graph")), table = {
    totalMatrix <- t(sequence)
    if (length(grouping) > 0) 
      totalMatrix = totalMatrix[unlist(grouping$index), , drop = FALSE]
    if (is.null(colnames(totalMatrix))) {
      colnames(totalMatrix) <- seq_len(ncol(totalMatrix)) 
    }
    else if (length(grep("t = ", colnames(totalMatrix)) == ncol(totalMatrix))) {
      colnames(totalMatrix) <- sapply(colnames(totalMatrix),
                                      gsub, 
                                      pattern = "t = ", replacement = "", fixed = TRUE)
    }
    if (!reverse) 
      totalMatrix <- totalMatrix[nrow(totalMatrix):1, , drop = F]
    plot(0, type="n", 
         xlim = c(0, ncol(totalMatrix)), ylim = c(-2, nrow(totalMatrix)+1), 
         xlab = "", ylab = "", axes = FALSE, main = title, ...)
    
    len <- ncol(totalMatrix)
    xStart <- 1
    unitFactor <- (len - (2 * eps))/len
    
    axis(3, eps, "t = ", 
         lty = "blank", yaxt = "s", xaxt = "s", xaxs = "i", pos = nrow(totalMatrix)+1)
    axis(3, eps + (1:ncol(totalMatrix))*unitFactor - unitFactor/2, ((1:ncol(totalMatrix)) - 1),
         lty = "blank", yaxt = "s", xaxt = "s", xaxs = "i", pos = nrow(totalMatrix)+1)
    axis(2, seq_len(nrow(totalMatrix)) - 0.5, rownames(totalMatrix),     yaxt = "s", xaxt = "s", xaxs = "i", las = 2)
    
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

    if (length(grouping) > 0) {
      if(!is.null(grouping$class)) {
        sepPos = cumsum(sapply(grouping$index, length))
        abline(h = sepPos[-length(sepPos)], col = "black", lwd = 3)
        text(ncol(totalMatrix)/2, sepPos - 0.5, grouping$class, cex = 0.9)
      }
    }
    if (!is.null(attributes(sequence)$attractor) && highlightAttractor) {
      attrStart <- min(attributes(sequence)$attractor) - 1
      #lines(x = c(attrStart, attrStart), y = c(-1, nrow(totalMatrix)) + 0.5)
      #abline(v =  eps+attrStart*unitFactor, lwd = attractor.sep.lwd, col = attractor.sep.col)
      rect(xleft = eps+attrStart*unitFactor, xright = eps+attrStart*unitFactor,
           ybottom = -0.5, ytop = nrow(totalMatrix)+0.75, lwd = attractor.sep.lwd, border = attractor.sep.col)
      arrows(x0 =  eps+attrStart*unitFactor, y0 = nrow(totalMatrix)+0.5,
             x1 = ncol(totalMatrix) - eps, y1 = nrow(totalMatrix)+0.5,
             length = 0.1, angle = 20, code = 3, col = attractor.sep.col)
      xpd.prev <- par()$xpd
      par(xpd=T)
      text(eps+attrStart*unitFactor, nrow(totalMatrix)+1, "Attractor", pos = 4)
      par(xpd=xpd.prev)
    }
    #if (drawLegend) legend(x = "bottomright", pch = c(15, 15), col = c(onColor, offColor), legend = c("active", "inactive"), cex = 0.7, horiz = T)
    if (drawLegend) legend(x = 0, y = -2, pch = c(15, 15), col = c(onColor, offColor), legend = c("active", "inactive"), cex = 0.7, horiz = T, xpd=T)
    #return(totalMatrix[network$genes,])
    return(totalMatrix)
  }, graph = {
    if (installed.packages()["igraph", "Version"] < package_version("0.6")) bias <- 1 else bias <- 0
    args <- list(...)
    if (is.null(args$vertex.size)) args$vertex.size <- 2
    if (is.null(args$edge.arrow.mode)) args$edge.arrow.mode <- 0
    if (is.null(args$rescale)) args$rescale <- !missing(layout)
    if (is.null(args$vertex.label.cex)) args$vertex.label.cex <- 0.75
    if (is.null(args$vertex.label.dist)) args$vertex.label.dist <- 0.25
    if (is.null(args$vertex.label.degree)) args$vertex.label.degree <- -pi/2
    if (is.null(args$vertex.color)) args$vertex.color <- "grey"
    if (is.null(args$edge.arrow.size)) args$edge.arrow.size <- 0.5
    if (missing(layout)) layout <- matrix(c(seq(-1, 1, length.out = nrow(sequence)),
                                            rep(0, nrow(sequence))), ncol = 2)
    states <- apply(sequence, 1, paste, collapse = "")
    nodes <- data.frame(seq_along(states), stringsAsFactors = FALSE)
    lastAttractorEdgeIndex <- NULL
    if (length(states) > 1) {
      initialStates <- 1:(length(states) - 1)
      nextStates <- 2:length(states)
      edgeMatrix <- data.frame(initialStates, nextStates)
      if (!is.null(attributes(sequence)$attractor)) {
        attractorEdge <- c(max(attributes(sequence)$attractor),
                           min(attributes(sequence)$attractor))
        if (all(attractorEdge <= length(states))) {
          edgeMatrix <- rbind(edgeMatrix, attractorEdge)
          lastAttractorEdgeIndex <- nrow(edgeMatrix)
        }
      }
    } else {
      edgeMatrix <- data.frame(matrix(nrow = 0, ncol = 2))
    }
    graph <- graph.data.frame(edgeMatrix, vertices = nodes,
                              directed = TRUE)
    if (drawLabels) labels <- states else labels <- NA
    if (highlightAttractor && !is.null(attributes(sequence)$attractor)) {
      attractorEdgeIndices <- intersect(seq_len(nrow(edgeMatrix)),
                                        c(attributes(sequence)$attractor, nrow(edgeMatrix))) -
        bias
      graph <- set.edge.attribute(graph, "width", index = attractorEdgeIndices,
                                  value = 3)
    }
    if (!is.null(lastAttractorEdgeIndex)) {
      graph <- set.edge.attribute(graph, "curved", value = 0)
      graph <- set.edge.attribute(graph, "curved", 
                                  index = lastAttractorEdgeIndex - bias, value = 0.5)
    }
    plot(graph, layout = layout, vertex.label = labels,
         vertex.label.cex = args$vertex.label.cex, vertex.size = args$vertex.size,
         vertex.color = args$vertex.color, vertex.label.dist = args$vertex.label.dist,
         vertex.label.degree = args$vertex.label.degree,
         edge.arrow.size = args$edge.arrow.size, rescale = args$rescale,
         main = title, ...)
    return(graph)
  })
}

