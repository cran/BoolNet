\name{plotPBNTransitions}
\alias{plotPBNTransitions}
\title{
Visualize the transitions in a probabilistic Boolean network
}
\description{
Visualizes the state transitions and their probabilities in a probabilistic Boolean network. This takes the transition table information calculated by the \code{\link{markovSimulation}} method. Only transitions with non-zero probability are included in the plot. The function requires the \pkg{igraph} package.
}
\usage{
plotPBNTransitions(markovSimulation, 
                   stateSubset, 
                   drawProbabilities = TRUE, 
                   drawStateLabels = TRUE, 
                   layout = layout.fruchterman.reingold, 
                   plotIt = TRUE, ...)
}
\arguments{
  \item{markovSimulation}{
  An object of class \code{MarkovSimulation}, as returned by \code{\link{markovSimulation}}. As the transition table information in this structure is required, \code{\link{markovSimulation}} must be called with \code{returnTable} set to TRUE.
}
  \item{stateSubset}{
An optional list of states, where each element of the list must be a vector with a 0/1 entry for each gene. If this argument is supplied, the graph only contains the specified states and transitions between these states.
}
  \item{drawProbabilities}{
If set to true, the edges of the graph are annotated with the probabilities of the corresponding transitions. Default is TRUE.
}
  \item{drawStateLabels}{
If set to true, the vertices of the graph are annotated with the gene values of the corresponding states. Defaults to TRUE.
}
  \item{layout}{
A layouting function that determines the placement of the nodes in the graph. Please refer to the \code{\link[igraph]{layout}} manual entry in the \pkg{igraph} package for further details. By default, the Fruchterman-Reingold algorithm is used.
}
  \item{plotIt}{
If this is true, a plot is generated. Otherwise, only an object of class \code{igraph} is returned, but no plot is drawn.
}
  \item{\dots}{
Further graphical parameters to be passed to \code{\link[igraph:plot.graph]{plot.igraph}}.
}
}

\value{
Returns an invisible object of class \code{igraph} containing the wiring graph.
}

\details{
This function uses the \code{\link[igraph:plot.graph]{plot.igraph}} function from the \pkg{igraph} package. The plots are customizeable using the \code{\dots} argument. For details on possible parameters, please refer to \code{\link[igraph:plot.common]{igraph.plotting}}.
}

\seealso{
\code{\link{markovSimulation}}
}
\examples{
# load example network
data(examplePBN)

# perform a Markov chain simulation
sim <- markovSimulation(examplePBN)

# plot the transitions and their probabilities
\donttest{plotPBNTransitions(sim)}
}

\keyword{probabilistic Boolean network
        PBN
        Boolean network
        plot
        transition
        graph}
