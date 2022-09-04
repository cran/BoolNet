\name{generateRandomNKNetwork}
\Rdversion{1.1}
\alias{generateRandomNKNetwork}
\title{
Generate a random N-K Boolean network
}
\description{
Generates a random N-K Boolean network (see Kauffman, 1969) using different configurations for the topology, the linkage, and the functions.
}
\usage{
generateRandomNKNetwork(n, k, 
                        topology = c("fixed", "homogeneous", "scale_free"), 
                        linkage = c("uniform", "lattice"), 
                        functionGeneration = c("uniform", "biased"), 
                        validationFunction, failureIterations=10000,
                        simplify = FALSE, noIrrelevantGenes=TRUE, 
                        readableFunctions = FALSE, 
                        d_lattice = 1, zeroBias = 0.5, 
                        gamma = 2.5, approx_cutoff = 100)
}
\arguments{
  \item{n}{
	The total number of genes in the network
}
  \item{k}{
	If this is a single number, this is either the maximum number of genes in the input of a transition function (for \code{topology="fixed"} and \code{topology="scale_free"}) or the mean number of genes in the input of a function (for \code{topology="homogeneous"}). If \code{topology="fixed"}, this can also be a vector with \code{n} elements specifying the number of input genes for each gene separately.
}
  \item{topology}{
	If set to "fixed", all transition functions of the network depend on exactly \code{k} input genes (unless there are irrelevant input genes to be removed if \code{simplify=TRUE} and \code{noIrrelevantGenes=FALSE}).
	
	If set to "homogeneous", the number of input genes is drawn independently at random from a Poisson distribution with lambda = k.
	
	If set to "scale_free", the number of input genes of each function is drawn from a Zeta distribution with parameter \code{gamma}.
}
  \item{linkage}{
	If this parameter is "uniform", the actual input genes are drawn uniformly at random from the total \code{k} genes.
	
	If set to "lattice", only genes from the neighbourhood \if{latex}{\cr}\code{(i - d_lattice * k_i):(i + d_lattice * k_i)} are taken, which means that all genes are dependent from other genes in the direct neighbourhood.
}
  \item{functionGeneration}{
  This parameter specifies how the truth tables of the transition functions are generated. 
	If set to "uniform", the truth table result column of the function is filled uniformly at random with 0 and 1. If set to "biased", a bias is introduced, where the probability of drawing a 0 is determined by the parameter \code{zeroBias}.
	
As a third option, \code{functionGeneration} can be set to a user-defined function that generates the truth tables. This function must have a single parameter \code{input} that is supplied with a vector of input gene indices. It must return a binary vector of size \code{2^length(input)} corresponding to the result column of the truth table. For the generation of canalyzing and nested canalyzing functions that are often assumed to be biologically plausible, the generation functions \code{\link{generateCanalyzing}} and \code{\link{generateNestedCanalyzing}} are included in \pkg{BoolNet}.
}

  \item{validationFunction}{
An optional function that restricts the generated Boolean functions to certain classes. This can be used if no explicit generation function can be specified in \code{functionGeneration}, but it is nevertheless possible to check whether a generated function belongs to that class or not. The function should have two input parameter \code{input} and \code{func} that receive a candidate function. \code{input} is a matrix of 0/1 integer values specifying the input part of the truth table of the candidate function, with the input genes in the columns. Each of the 2^k rows of \code{input}  (where k is the number of input genes) corresponds to one entry of \code{func}, which is an integer vector of 0/1 values corresponding to the output of the candidate function. The validation function should return \code{TRUE} if the candidate function is accepted or \code{FALSE} if it is rejected.
}

  \item{failureIterations}{
The maximum number of iterations the generator tries to generate a function that is accepted by \code{validationFunction} before it gives up and throws an error. Defaults to 10000.  
}

  \item{simplify}{
	If this is true, \code{\link{simplifyNetwork}} is called to simplify the gene transition functions after the perturbation. This removes irrelevant input genes. Should not be used together with \code{noIrrelevantGenes=TRUE}, as this automatically generates a network that cannot be simplified any further. Defaults to FALSE.
}
  \item{noIrrelevantGenes}{
  If set to true, gene transition functions are not allowed to contain irrelevant genes, i.e. the functions have exactly the number of input genes determined by the \code{topology} method. This means that the network cannot be simplified any further, and \code{simplify} should be turned off. The default value is TRUE.
}

  \item{readableFunctions}{
This parameter specifies if readable DNF representations of the transition function truth tables are generated and displayed when the network is printed. If set to FALSE, the truth table result column is displayed. If set to "canonical", a canonical Disjunctive Normal Form is generated from each truth table. If set to "short", the canonical DNF is minimized by joining terms (which can be time-consuming for functions with many inputs). If set to TRUE, a short DNF is generated for functions with up to 12 inputs, and a canonical DNF is generated for functions with more than 12 inputs.
}

  \item{d_lattice}{
	The dimension parameter for the lattice if \code{linkage="lattice"}. Defaults to 1.
}
  \item{zeroBias}{
	The bias parameter for biased functions for \code{functionGeneration="biased"}. Defaults to 0.5 (no bias).
}
  \item{gamma}{
	The Gamma parameter of the Zeta distribution for \code{topology="scale_free"}. Default is 2.5.
}
  \item{approx_cutoff}{
	This parameter is only used with \code{topology="scale_free"}. It sets the number of iterations in the sum used to approximate the Riemann Zeta function. Defaults to 100.
}
}
\details{
The function supports a high number of different configurations to generate random networks. Several of the parameters are only needed for special configurations. The generated networks have different structural properties. Refer to the literature for more details.

Constant genes are automatically fixed (e.g. knocked-out or over-expressed). This means that they are always set to the constant value, and states with the complementary value are not considered in transition tables etc. If you would like to change this behaviour, use \code{\link{fixGenes}} to reset the fixing.
}
\value{
An object of class \code{BooleanNetwork} containing the generated random network. The class \if{latex}{\cr}\code{BooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\references{
S. A. Kauffman (1969), Metabolic stability and epigenesis in randomly constructed nets. J. Theor. Biol. 22:437--467.

S. A. Kauffman (1993), The Origins of Order. Oxford University Press.

M. Aldana (2003), Boolean dynamics of networks with scale-free topology. Physica D 185: 45--66.

M. Aldana and S. Coppersmith and L. P. Kadanoff (2003), Boolean dynamics with random coupling. In E. Kaplan, J. E. Marsden  and K. R. Sreenivasan (editors): Perspectives and Problems in Nonlinear Science, Springer.
}
\seealso{
\code{\link{perturbNetwork}},\code{\link{loadNetwork}}, \code{\link{simplifyNetwork}}, \code{\link{fixGenes}}
}
\examples{
\dontrun{
# generate different random networks
net1 <- generateRandomNKNetwork(n=10, k=10, 
                                topology="scale_free", 
                                linkage="uniform", 
                                functionGeneration="uniform",
                                noIrrelevantGenes=FALSE, 
                                simplify=TRUE)
                                
net2 <- generateRandomNKNetwork(n=10, k=3, 
                                topology="homogeneous", 
                                linkage="lattice", 
                                functionGeneration="uniform", 
                                d_lattice=1.5, 
                                simplify=TRUE)
                                
net3 <- generateRandomNKNetwork(n=10, k=2, 
                                topology="fixed", 
                                linkage="uniform", 
                                functionGeneration="biased", 
                                noIrrelevantGenes=FALSE,                                
                                zeroBias=0.6)

# get attractors
print(getAttractors(net1))
print(getAttractors(net2))
print(getAttractors(net3))
}
}