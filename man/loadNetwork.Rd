\name{loadNetwork}
\Rdversion{1.1}
\alias{loadNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load a Boolean network from a file
}
\description{
Loads a Boolean network or probabilistic Boolean network from a file and converts it to an internal transition table representation.
}
\usage{
loadNetwork(file, 
            bodySeparator = ",",
            lowercaseGenes = FALSE)
}

\arguments{
  \item{file}{
	The name of the file to be read
}

  \item{bodySeparator}{
	An optional separation character to divide the target factors and the formulas. Default is ",".
}

  \item{lowercaseGenes}{
    If set to \code{TRUE}, all gene names are converted to lower case, i.e. the gene names are case-insensitive. This corresponds to the behaviour of \pkg{BoolNet} versions prior to 1.5. Defaults to \code{FALSE}.
  }
}
\details{
The routine supports a description language based on the Boolean operators AND (&), or (|), and NOT (!). 
The first line contains a header. In case of a Boolean network with only one function per gene, the header is "targets, functions"; in a probabilistic network, there is an optional third column "probabilities". All subsequent lines contain Boolean rules.
A rule consists of a target gene, a separator, a Boolean expression to calculate a transition step for the target gene, and an optional probability for the rule (for probabilistic Boolean networks only -- see below). The EBNF description of the rule format is as follows:

\preformatted{Rule = GeneName Separator BooleanExpression [Separator Probability];
BooleanExpression = GeneName 
                    | "!" BooleanExpression 
                    | "(" BooleanExpression ")" 
                    | BooleanExpression " & " BooleanExpression  
                    | BooleanExpression " | " BooleanExpression;
GeneName = ? A gene name from the list of involved genes ?;
Separator = ? The field <bodySeparator>, usually "," ?;
Probability = ? A floating-point number ?;
}

If there is exactly one rule for each gene, a Boolean network of class \code{BooleanNetwork} is created. In these networks, constant genes are automatically fixed (e.g. knocked-out or over-expressed). This means that they are always set to the constant value, and states with the complementary value are not considered in transition tables etc. If you would like to change this behaviour, use \code{\link{fixGenes}} to reset the fixing.

If two or more rules exist for the same gene, the function returns a probabilistic network of class \code{ProbabilisticBooleanNetwork}. In this case, alternative rules may be annotated with probabilities, which must sum up to 1 for all rules that belong to the same gene. If no probabilities are supplied, uniform distribution is assumed.
}
\value{
If only one function per gene is specified, a structure of class \code{BooleanNetwork} representing the network is returned. It has the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the transition function for the i-th gene. Each element has the following sub-components:
  \describe{
  	\item{input}{A vector of indices of the genes that serve as the input of the Boolean transition function. If the function has no input (i.e. the gene is constant), the vector consists of a zero element.}
  	\item{func}{The transition function in truth table representation. This vector has \if{latex}{\cr}\code{2^length(input)} entries, one for each combination of input variables. If the gene is constant, the function is 1 or 0.}
  	\item{expression}{A string representation of the Boolean expression from which the truth table was generated}
  }}
  \item{fixed}{A vector specifying which genes are knocked-out or over-expressed. For each gene, there is one element which is set to 0 if the gene is knocked-out, to 1 if the gene is over-expressed, and to -1 if the gene is not fixed at all, i. e. can change its value according to the supplied transition function. Constant genes are automatically set to fixed values.}

If there is at least one gene with two or more alternative transition functions, a structure of class \code{ProbabilisticBooleanNetwork} is returned. This structure is similar to \code{BooleanNetwork}, but allows for storing more than one function in an interaction. It consists of the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the alternative transition functions for the i-th gene. Each element is a list of transition functions. In this second-level list, each element has the the following sub-components:
  \describe{
	  \item{input}{A vector of indices of the genes that serve as the input of the Boolean transition function. If the function has no input (i.e. the gene is constant), the vector consists of a zero element.}
	  \item{func}{The transition function in truth table representation. This vector has \if{latex}{\cr}\code{2^length(input)} entries, one for each combination of input variables. If the gene is constant, the function is -1.}
	  \item{expression}{A string representation of the underlying Boolean expression}
	  \item{probability}{The probability that the corresponding transition function is chosen}	  
  }}
  \item{fixed}{A vector specifying which genes are knocked-out or over-expressed. For each gene, there is one element which is set to 0 if the gene is knocked-out, to 1 if the gene is over-expressed, and to -1 if the gene is not fixed at all, i. e. can change its value according to the supplied transition function. You can knock-out and over-express genes using \code{\link{fixGenes}}.}
}
\seealso{
\code{\link{getAttractors}}, \code{\link{markovSimulation}}, \code{\link{toPajek}}, \code{\link{stateTransition}}, \code{\link{fixGenes}}
}
\examples{
# write example network to file
sink("testNet.bn")
cat("targets, functions\n")
cat("Gene1, !Gene2 | !Gene3\n")
cat("Gene2, Gene3 & Gene4\n")
cat("Gene3, Gene2 & !Gene1\n")
cat("Gene4, 1\n")
sink()

# read file
net <- loadNetwork("testNet.bn")
print(net)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Boolean network
  probabilistic Boolean network
  PBN
	file
	logic
	parse}
