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
            lowercaseGenes = FALSE,
            symbolic = FALSE)
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
  
  \item{symbolic}{
If set to \code{TRUE}, a symbolic representation of class \code{SymbolicBooleanNetwork} is returned. This is not available for asynchronous or probabilistic Boolean networks, but is required for the simulation of networks with extended temporal predicates and time delays (see \code{\link{simulateSymbolicModel}}). If such predicates are detected, the switch is activated by default. 
}
  
}
\details{

Depending on whether the network is loaded in truth table representation or not, the supported network file formats differ slightly.

For the truth table representation (\code{symbolic=FALSE}), the language basically consists of expressions based on the Boolean operators AND (&), or (|), and NOT (!). In addition, some convenience operators are included (see EBNF and operator description below).
The first line contains a header. In case of a Boolean network with only one function per gene, the header is "targets, functions"; in a probabilistic network, there is an optional third column "probabilities". All subsequent lines contain Boolean rules or comment lines that are omitted by the parser.
A rule consists of a target gene, a separator, a Boolean expression to calculate a transition step for the target gene, and an optional probability for the rule (for probabilistic Boolean networks only -- see below). 

The EBNF description of the network file format is as follows:

\preformatted{
Network           = Header Newline {Rule Newline | Comment Newline};
Header            = "targets" Separator "factors";
Rule              = GeneName Separator BooleanExpression [Separator Probability];
Comment           = "#" String;
BooleanExpression = GeneName 
                    | "!" BooleanExpression 
                    | "(" BooleanExpression ")" 
                    | BooleanExpression " & " BooleanExpression  
                    | BooleanExpression " | " BooleanExpression;
                    | "all(" BooleanExpression {"," BooleanExpression} ")"
                    | "any(" BooleanExpression {"," BooleanExpression} ")"
                    | "maj(" BooleanExpression {"," BooleanExpression} ")"
                    | "sumgt(" BooleanExpression {"," BooleanExpression} "," Integer ")"
                    | "sumlt(" BooleanExpression {"," BooleanExpression} "," Integer ")";
GeneName          = ? A gene name from the list of involved genes ?;
Separator         = ",";
Integer           = ? An integer value?;
Probability       = ? A floating-point number ?;
String            = ? Any sequence of characters (except a line break) ?;
Newline           = ? A line break character ?;
}

The extended format for Boolean networks with temporal elements that can be loaded if \code{symbolic=TRUE} additionally allows for a specification of time steps. Furthermore, the operators can be extended with iterators that evaluate their arguments over multiple time steps.

\preformatted{
Network               = Header Newline
                        {Function Newline | Comment Newline};
Header                = "targets" Separator "factors";
Function              = GeneName Separator BooleanExpression;
Comment               = "#" String;
BooleanExpression     = GeneName | GeneName TemporalSpecification | BooleanOperator | TemporalOperator
BooleanOperator       =   BooleanExpression 
                        | "!" BooleanExpression 
                        | "(" BooleanExpression ")" 
                        | BooleanExpression " & " BooleanExpression  
                        | BooleanExpression " | " BooleanExpression;
TemporalOperator      =   "all" [TemporalIteratorDef] 
                                "(" BooleanExpression {"," BooleanExpression} ")"
                        | "any" [TemporalIteratorDef] 
                                "(" BooleanExpression {"," BooleanExpression} ")"
                        | "maj" [TemporalIteratorDef] 
                                "(" BooleanExpression {"," BooleanExpression} ")"
                        | "sumgt" [TemporalIteratorDef] 
                                  "(" BooleanExpression {"," BooleanExpression} "," Integer ")"
                        | "sumlt" [TemporalIteratorDef] 
                                  "(" BooleanExpression {"," BooleanExpression} "," Integer ")"
                        | "timeis" "(" Integer ")"
                        | "timegt" "(" Integer ")"
                        | "timelt" "(" Integer ")";
TemporalIteratorDef   = "[" TemporalIterator "=" Integer ".." Integer "]";
TemporalSpecification = "[" TemporalOperand {"+" TemporalOperand | "-" TemporalOperand} "]";
TemporalOperand       = TemporalIterator | Integer
TemporalIterator      = ? An alphanumeric string ?;
GeneName              = ? A gene name from the list of involved genes ?;
Separator             = ",";
Integer               = ? An integer value?;
String                = ? Any sequence of characters (except a line break) ?;
Newline               = ? A line break character ?;
}

The meaning of the operators is as follows:
\describe{

\item{\code{all}}{Equivalent to a conjunction of all arguments. For symbolic networks, the operator can have a time range, in which case the arguments are evaluated for each time point specified in the iterator.}

\item{\code{any}}{Equivalent to a disjunction of all arguments. For symbolic networks, the operator can have a time range, in which case the arguments are evaluated for each time point specified in the iterator.}

\item{\code{maj}}{Evaluates to true if the majority of the arguments evaluate to true. For symbolic networks, the operator can have a time range, in which case the arguments are evaluated for each time point specified in the iterator.}

\item{\code{sumgt}}{Evaluates to true if the number of arguments (except the last) that evaluate to true is greater than the number specified in the last argument. For symbolic networks, the operator can have a time range, in which case the arguments are evaluated for each time point specified in the iterator.}

\item{\code{sumlt}}{Evaluates to true if the number of arguments (except the last) that evaluate to true is less than the number specified in the last argument. For symbolic networks, the operator can have a time range, in which case the arguments are evaluated for each time point specified in the iterator.}

\item{\code{timeis}}{Evaluates to true if the current absolute time step (i.e. number of state transitions performed from the current start state) is the same as the argument.}

\item{\code{timelt}}{Evaluates to true if the current absolute time step (i.e. number of state transitions performed from the current start state) is the less than the argument.}

\item{\code{timegt}}{Evaluates to true if the current absolute time step (i.e. number of state transitions performed from the current start state) is greater than the argument.}

}

If \code{symbolic=FALSE} and there is exactly one rule for each gene, a Boolean network of class \code{BooleanNetwork} is created. In these networks, constant genes are automatically fixed (e.g. knocked-out or over-expressed). This means that they are always set to the constant value, and states with the complementary value are not considered in transition tables etc. If you would like to change this behaviour, use \code{\link{fixGenes}} to reset the fixing.

If \code{symbolic=FALSE} and two or more rules exist for the same gene, the function returns a probabilistic network of class \code{ProbabilisticBooleanNetwork}. In this case, alternative rules may be annotated with probabilities, which must sum up to 1 for all rules that belong to the same gene. If no probabilities are supplied, uniform distribution is assumed.

If \code{symbolic=TRUE}, a symbolic representation of a (possibly temporal) Boolean network of class \code{SymbolicBooleanNetwork} is created.
}
\value{
If \code{symbolic=FALSE} and only one function per gene is specified, a structure of class \code{BooleanNetwork} representing the network is returned. It has the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the transition function for the i-th gene. Each element has the following sub-components:
  \describe{
  	\item{input}{A vector of indices of the genes that serve as the input of the Boolean transition function. If the function has no input (i.e. the gene is constant), the vector consists of a zero element.}
  	\item{func}{The transition function in truth table representation. This vector has \if{latex}{\cr}\code{2^length(input)} entries, one for each combination of input variables. If the gene is constant, the function is 1 or 0.}
  	\item{expression}{A string representation of the Boolean expression from which the truth table was generated}
  }}
  \item{fixed}{A vector specifying which genes are knocked-out or over-expressed. For each gene, there is one element which is set to 0 if the gene is knocked-out, to 1 if the gene is over-expressed, and to -1 if the gene is not fixed at all, i. e. can change its value according to the supplied transition function. Constant genes are automatically set to fixed values.}

If \code{symbolic=FALSE} and there is at least one gene with two or more alternative transition functions, a structure of class \code{ProbabilisticBooleanNetwork} is returned. This structure is similar to \code{BooleanNetwork}, but allows for storing more than one function in an interaction. It consists of the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the alternative transition functions for the i-th gene. Each element is a list of transition functions. In this second-level list, each element has the the following sub-components:
  \describe{
	  \item{input}{A vector of indices of the genes that serve as the input of the Boolean transition function. If the function has no input (i.e. the gene is constant), the vector consists of a zero element.}
	  \item{func}{The transition function in truth table representation. This vector has \if{latex}{\cr}\code{2^length(input)} entries, one for each combination of input variables. If the gene is constant, the function is -1.}
	  \item{expression}{A string representation of the underlying Boolean expression}
	  \item{probability}{The probability that the corresponding transition function is chosen}	  
  }}
  \item{fixed}{A vector specifying which genes are knocked-out or over-expressed. For each gene, there is one element which is set to 0 if the gene is knocked-out, to 1 if the gene is over-expressed, and to -1 if the gene is not fixed at all, i. e. can change its value according to the supplied transition function. You can knock-out and over-express genes using \code{\link{fixGenes}}.}
  
If \code{symbolic=TRUE}, a structure of class \code{SymbolicBooleanNetwork} that represents the network as expression trees is returned. It has the following components:
  \item{genes}{A vector of gene names involved in the network. This list determines the indices of genes in inputs of functions or in state bit vectors.}
  \item{interactions}{A list with \code{length(genes)} elements, where the i-th element describes the transition function for the i-th gene in a symbolic representation. Each such element is a list that represents a recursive expression tree, possibly consisting of sub-elements (operands) that are expression trees themselves. Each element in an expression tree can be a Boolean/temporal operator, a literal ("atom") or a numeric constant.}
  \item{internalStructs}{A pointer referencing an internal representation of the expression trees as raw C objects. This is used for simulations and must be set to NULL if \code{interactions} are changed to force a refreshment.
  }
  \item{timeDelays}{An integer vector storing the temporal memory sizes required for each of the genes in the network. That is, the vector stores the minimum number of predecessor states of each gene that need to be saved to determine the successor state of the network.}
  \item{fixed}{A vector specifying which genes are knocked-out or over-expressed. For each gene, there is one element which is set to 0 if the gene is knocked-out, to 1 if the gene is over-expressed, and to -1 if the gene is not fixed at all, i. e. can change its value according to the supplied transition function. Constant genes are automatically set to fixed values.}  
}
\seealso{
\code{\link{getAttractors}}, \code{\link{simulateSymbolicModel}}, \code{\link{markovSimulation}}, \code{\link{stateTransition}}, \code{\link{fixGenes}}, \code{\link{loadSBML}}, \code{\link{loadBioTapestry}}
}
\examples{
# write example network to file
sink("testNet.bn")
cat("targets, factors\n")
cat("Gene1, !Gene2 | !Gene3\n")
cat("Gene2, Gene3 & Gene4\n")
cat("Gene3, Gene2 & !Gene1\n")
cat("Gene4, 1\n")
sink()

# read file
net <- loadNetwork("testNet.bn")
print(net)
}
\keyword{Boolean network
  probabilistic Boolean network
  PBN
	file
	logic
	parse}
