# Export a state table in <attractorInfo> to a Pajek graph
toPajek <- function (attractorInfo, file="boolean.net", includeLabels=FALSE) 
{
  stopifnot(inherits(attractorInfo,"AttractorInfo"))
  
  if (is.null(attractorInfo$stateInfo$table))
    stop(paste("This AttractorInfo structure does not contain transition table information.",
           "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))

  graphStruct <- getStateGraphStructure(attractorInfo)
  
  sink(file)
  cat("*Vertices ", length(graphStruct$vertices), "\r\n", sep = "")
  if (includeLabels)
  {
    lapply(seq_along(graphStruct$vertices),function(i)
      cat(i," \"",graphStruct$vertices[i],"\"\r\n",sep=""))
  }
  cat("*Arcs\r\n")
  apply(graphStruct$edges,1,function(edge)
    cat(edge[1]," ",edge[2]," 1\r\n",sep=""))
  sink()
}
