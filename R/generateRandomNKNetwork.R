generateRandomNKNetwork <- function(n,k,topology=c("fixed","homogeneous","scale_free"),
					linkage=c("uniform","lattice"),
					functionGeneration=c("uniform","biased"),
					simplify=FALSE, noIrrelevantGenes=TRUE, readableFunctions=FALSE,
					d_lattice=1, zeroBias=0.5, gamma=2.5, approx_cutoff=100)
{
	k_i_vec <- switch(match.arg(topology),
			fixed = {
					if (length(k) == n)
						k
					else
					if (length(k) == 1)
						rep(k,n)
					else
						stop("k must have 1 or n element(s)!")
				},
			homogeneous = round(rpois(n,k)),
			scale_free = rzeta(n,k,gamma=gamma,approx_cutoff=approx_cutoff),
			stop("'topology' must be one of \"fixed\",\"homogeneous\"")
			)
			
	k_i_vec[k_i_vec > n] <- n
	
	geneNames <- paste("Gene",1:n)
			
	interactions <- mapply(function(i,k_i)
			{
				if (k_i == 0)
				{
					genes <- 0
					func <- round(runif(min=0,max=1,n=1))
				}
				else
				{
					genes <- switch(match.arg(linkage,c("uniform","lattice")),
						uniform = sample(1:n,k_i,replace=FALSE),
						lattice = {
								region <- c(max(1,round(i - k_i*d_lattice)):max(1,i-1),
									    min(n,i+1):min(n,round(i + k_i*d_lattice)))
								
								sample(region,k_i,replace=FALSE)
							},
						stop("'linkage' must be one of \"uniform\",\"lattice\""))

					  containsIrrelevant <- TRUE
					  while(containsIrrelevant)
					  {
							  func <- switch(match.arg(functionGeneration,c("uniform","biased")),
									  uniform = round(runif(min=0,max=1,n=2^k_i)),
									  biased = as.integer(runif(min=0,max=1,n=2^k_i) > zeroBias),
									  stop("'functionGeneration' must be one of \"uniform\",\"biased\""))

								if (noIrrelevantGenes)
								{
									table <- allcombn(2,k_i) - 1

									dropGenes <- apply(table,2,function(gene)
									# determine all genes that have no influence on the results,
									# i.e. the result column is equal for 0 values and 1 values
												  {
													  (identical(func[gene==1],
															 	 func[gene==0]))
												  })
																		
									if (sum(dropGenes) == 0)
										containsIrrelevant <- FALSE
									
								}
								else
									containsIrrelevant <- FALSE
						}


				}
				return(list(input=genes,func=func,
					expression=getInteractionString(readableFunctions,
							    	     func,geneNames[genes])))
			},1:length(k_i_vec),k_i_vec,SIMPLIFY=FALSE)
	fixed <- sapply(interactions,function(i)
			{
				if (i$input[1] == 0)
					i$func[1]
				else
					-1
			})
	net <- list(genes=geneNames,interactions=interactions,fixed=fixed)
	class(net) <- "BooleanNetwork"
	if (simplify)
		net <- simplifyNetwork(net,readableFunctions)
	return(net)	
}
