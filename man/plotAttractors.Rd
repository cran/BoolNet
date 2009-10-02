\name{plotAttractors}
\Rdversion{1.1}
\alias{plotAttractors}

\title{
Plot state table of attractors
}
\description{
Visualizes attractors by drawing a table of the involved states in two colors.
}
\usage{
plotAttractors(attractorInfo, 
               grouping = list(), 
               title = "", 
               plotFixed = TRUE, 
               onColor = "green", offColor = "red")
}
\arguments{
  \item{attractorInfo}{
	An object of class \code{AttractorInfo}, as returned by \code{\link{getAttractors}}
}
  \item{grouping}{An optional structure to form groups of genes in the plot. This is a list with the following elements: 
 \describe{
  	\item{class}{A vector of names for the groups. These names will be printed in the region belonging to the group in the plot.}
  	\item{index}{A list with the same length as \code{class}. Each element is a vector of gene indices belonging to the group.}}
}
 
  \item{title}{
	An optional title for the plot
}
  \item{plotFixed}{
	If this is true, genes with fixed values are included in the plot. Otherwise, these genes are not drawn.
}
  \item{onColor}{
	An optional color value for the 1/ON values in the table. Defaults to green.
}
  \item{offColor}{
	An optional color value for the 0/OFF values in the table. Defaults to red.
}
}
\details{
Attractors in \code{attractorInfo} are first grouped by length. Then, a figure is plotted to the currently selected device for each attractor length (i.e. one plot with all attractors consisting of 1 state, one plot with all attractors consisting of 2 states, etc.). This means that on the standard X11 output device, only the last plot is displayed unless you set \code{par(mfrow=c(...))} accordingly.
The figure is a table with the genes in the rows and the states of the attractors in the columns. Cells of the table are (by default) red for 0/OFF values and green for 1/ON values. If \code{grouping} is set, the genes are rearranged according to the indices in the group, horizontal separation lines are plotted between the groups, and the group names are printed.
}
\value{
A list of matrices corresponding to the plots is returned. Each of these matrices has the genes in the rows and the states of the attractors in the columns.
}
\seealso{
\code{\link{getAttractors}}, \code{\link{attractorsToLaTeX}}
}
\examples{
library(BoolNet)

# load example data
data(cellcycle)

# get attractors
attractors <- getAttractors(net)

# plot attractors side by side
par(mfrow=c(1,length(attractors$attractors)))
plotAttractors(attractors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LaTeX
	Boolean network
	attractor
	cycle
	basin
	transition
	state
	graph
	plot}
