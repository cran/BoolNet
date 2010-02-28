\name{yeastTimeSeries}
\alias{yeastTimeSeries}
\docType{data}
\title{
Yeast cell cycle time series data
}
\description{
Preprocessed time series measurements of four genes from the yeast cell cycle data by Spellman et al.
}
\usage{data(yeastTimeSeries)}
\format{
  A matrix with 14 measurements for the genes Fhk2, Swi5, Sic1, and Clb1. Each gene is a row of the matrix, and each column is a measurement.
}
\details{
The data were obtained from the web site of the yeast cell cycle analysis project at \url{http://genome-www.stanford.edu/cellcycle}. The time series synchronized with the elutriation method were extracted for the genes Fhk2, Swi5, SIC1, and Clb1. In a preprocessing step, missing values were imputed by taking the means of the measurements of the same genes at neighbouring time points.
}
\source{
P. T. Spellman, G. Sherlock, M. Q. Zhang, V. R. Iyer, K. Anders, M. B. Eisen, P. O. Brown, D. Botstein, B. Futcher (1998), Comprehensive Identification of Cell Cycle-regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray Hybridization. Molecular Biology of the Cell 9(12):3273--3297.

Yeast cell cycle analysis project web site: \url{http://genome-www.stanford.edu/cellcycle}
}
\examples{
library(BoolNet)
data(yeastTimeSeries)

# the data set is stored in a variable called 'yeastTimeSeries'
print(yeastTimeSeries)
}
\keyword{datasets
        yeast cell cycle
        time series
        binarization
        binarize}
