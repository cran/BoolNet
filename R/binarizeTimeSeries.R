# A simple binarization for time series. <measurements> is a list of matrices with the
# genes in the rows, each specifying one time series.
# <nstart> and <iter.max> are the corresponding parameters for k-means. This clustering
# is employed for binarization. If <dropInsignificant> is true, a Wilcoxon test on the
# binarized genes is performed, and insignificant genes (with p < <alpha>) are removed.
# Returns a list containing the binarized matrix list, a vector of p-values, and a vector
# of thresholds.
binarizeTimeSeries <- function(measurements,nstart=100,iter.max=1000,dropInsignificant=TRUE,alpha=0.05)
{
	if (!is.null(dim(measurements)))
		fullData <- measurements
	else
	# in case of list, paste all matrices before clustering
	{
		fullData <- measurements[[1]]
		for (m in measurements[-1])
		{
			fullData <- cbind(fullData,m)
		}
	}

	cluster <- apply(fullData,1,function(gene)
	# cluster data using k-means
	{
		cl_res <- kmeans(gene, 2, nstart=nstart,iter.max=iter.max)
		
		if (cl_res$centers[1] > cl_res$centers[2])
		# exchange clusters if necessary, so that smaller numbers
		# are binarized to 0, and larger numbers are binarized to 1
			group <- abs(cl_res$cluster-2)
		else
			group <- cl_res$cluster-1
		
		# check significance of partitioning	
		pval <- wilcox.test(gene[group==0], gene[group==1], exact=FALSE)$p.value
		
		# calculate the binarization threshold
		threshold <-  min(cl_res$centers) + dist(cl_res$centers)[1]/2
		list(bin=group, threshold=threshold, p.value=pval)
	})
			
	# determine significant genes
	significant <- sapply(cluster,function(cl)(cl$p.value<alpha))
		
	if (any(!significant))
	{
		if (dropInsignificant)
		# remove insignificant genes
			cluster <- cluster[significant]
		else
		# throw warning
			warning("The following genes show a uniform behaviour and should possibly be excluded from binarization:", paste(rownames(fullData)[!significant],collapse=" "))
	}
	
	if (is.null(dim(measurements)))
	# split up the collated binarized measurements into a list of matrices of the original size
	{
		startIndex <- 0
		binarizedTimeSeries <- lapply(measurements,function(m)
					{
						currentSplit <- (startIndex+1):(startIndex+ncol(m))
						startIndex <<- startIndex + ncol(m)
						t(sapply(cluster,function(cl)
							cl$bin[currentSplit]))					
					})
	}
	else
		binarizedTimeSeries <- 	t(sapply(cluster,function(cl)
							cl$bin))

	return(list(binarizedMeasurements=binarizedTimeSeries,
		    pvals=sapply(cluster,function(cl)cl$p.value),
		    thresholds=sapply(cluster,function(cl)cl$threshold)))
}
