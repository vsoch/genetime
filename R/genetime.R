#' genetime
#'
#' This function will predict "important" time periods for a list of gene probes based on RNA-seq expression.
#' @param genelist A list of gene symbol identifiers
#' @param thresh The quantile of expression values to take per gene.  Default is 0.95
#' @param iters Number of iterations for permutations to create null distributions,  Default is 1000
#' @param bf bonferonni correction to use (Default is 0.05)
#' @param out.pdf output pdf of temporal enrichment patterns, defaults to FALSE.  Default output name (set with out.name) is "genetime.pdf"
#' @param out.name name of output pdf of temporal enrichment patterns, default is "genetime.pdf"
#' @keywords temporal gene enrichment
#' @export
#' @examples
#' mygenes = "ADNP" "AGAP2" "ANK2" "APH1A"
#' thresh = 0.99
#' genetime(genelist=mygenes,thresh=thresh,out.pdf=TRUE,out.name="/home/myfile.pdf")

genetime = function(genelist,thresh=0.95,out.pdf=FALSE,out.name="genetime.pdf",iters=1000,bf=0.05){

  library(pheatmap)
  library(fields)
  library(MASS)
  
  # So we can plot outside of main plot area
  par(xpd=TRUE,new=FALSE)
  
  # First get raw times, regions, and genes for the set
  raw = probePeaks(genelist,thresh=thresh,out.pdf=out.pdf,out.name=out.name)

  # Get our region colors
  colors = getColors()
  
  # Here is our x axis
  timepoints = raw$timepoints

  # Genetime matrix
  genetimes = raw$matrix
  
  # We will weight the timepoint by the summed RNA expression values
  weights = array(0,length(timepoints))
  for (tt in 1:length(timepoints)){
    t = timepoints[tt]
    subset = genetimes[which(genetimes$timepoint == t),]
    weights[tt] = sum(subset$value)
  }

  names(weights) = timepoints
  par(mfrow=c(1,1))
  barplot(weights,main="weighted probes per timepoint, all regions, (sum of expression values)",col="tomato")

  cat ("Press [enter] to continue")
  line = readline()
  
  # Now make a matrix of scores for each region
  uniqueregions = unique(genetimes$region)

  # Let's make a matrix of these counts for each region
  countmatrix = array(0,dim=c(length(uniqueregions),length(timepoints)))
  rownames(countmatrix) = uniqueregions
  colnames(countmatrix) = timepoints
  for (rr in 1:length(uniqueregions)){
    region = uniqueregions[rr]
    subset = genetimes[which(genetimes$region==region),]
    vector = array(0,length(timepoints))
    for (tt in 1:length(timepoints)){
      t = timepoints[tt]
      subsett = subset[which(subset$timepoint == t),]
      vector[tt] = sum(subsett$value)  
    }
    countmatrix[rr,] = vector  
  }

  # Here is a non corrected "enrichment" graphic for all regions and timepoints
  ph = pheatmap(countmatrix,cluster_cols=FALSE,main="Uncorrected temporal enrichment map")
  cat ("Press [enter] to continue")
  line = readline()

  # This function will rearrange matrix to be same as heatmap (regions clustered, not timepoints)
  formatMatrix = function(matrix){
    matrix = matrix[ph$tree_row$order,]
    rownames(matrix) = ph$tree_row$labels
    colnames(matrix) = timepoints
    return(matrix)
  }
  
  # We will save a matrix of pvalues
  pvalue_matrix = array(dim=c(length(uniqueregions),length(timepoints)))
  qvalues = array(dim=c(length(uniqueregions),length(timepoints)))
  
  cat("Running permutations...\n")
  
  # We are going to save p values for how well the null distributions fit negative binomial
  negative_binomial_p = c()
  
  # Create a null distribution for each region
  for (rr in 1:length(uniqueregions)){
    cat("region",rr,"of",length(uniqueregions),"\n")
    region = uniqueregions[rr]
    subset = genetimes[which(genetimes$region==region),]
    nulldist = array(0,dim=c(iters,length(timepoints)))
    colnames(nulldist) = timepoints
    for (i in 1:iters) {
      # Get a random sample of timepoints
      shuffled = sample(timepoints,nrow(subset),replace=TRUE)
      subsett = subset    
      subsett$timepoint = shuffled
      for (tt in 1:length(timepoints)) {
        t = timepoints[tt]
        shuffled_subset = subsett[which(subsett$timepoint == t),]
        nulldist[i,tt] = sum(shuffled_subset$value)    
      }
    }
  
    # Fit null to a negative binomial, and calculate a p-value
    # Save pvalues in a vector to correct for 16 tests
    pvalues = c()
    for (tt in 1:length(timepoints)){
      single_null_dist = nulldist[,tt]
      #hist(single_null_dist,col="purple",main=paste("Null Distribution timepoint",timepoints[tt]))
      #lines(rep(confidence_intervals[1],100),seq(1,100),col="red",lwd=2)
      #lines(rep(confidence_intervals[2],100),seq(1,100),col="red",lwd=2)
      # Use glm.nb to fit a negative binomial to one of our nulls
      # Round so we are using counts
      nully = round(single_null_dist)
      out = glm.nb(nully~1,link=identity)
      # Here is the intercept
      # coef(out)
      # theta: this is the "dispersion parameter" which must control how wide it is?
      theta = out$theta
      summary = summary(out)
      
      # I *think* we can use this pvalue to see how well the model fits the data?
      # (this confirms the null distrubution is modeled well by the binomial?)
      # Let's save this value
      negative_binomial_p = c(negative_binomial_p,summary$coefficients[4])
      
      # Now create a density with the      
      bin = dnbinom(1:max(single_null_dist), mu=mean(single_null_dist), size=theta)
      #par(mfrow=c(1,2))
      #hist(nully,col="tomato",main="Actual null distribution")
      
      # Now calculate where the actual value is
      actual = countmatrix[rr,tt]
      
      # I am SO close - I know I need to use the pnbinom function to get a pvalue, but I'm not sure what the input is.
      # Get p value, single sided test
      pvalue = pnbinom(actual, mu=mean(single_null_dist), size=theta,lower.tail=FALSE)      
      pvalues = c(pvalues,pvalue)
      #plot(bin,main=paste("NB Model w/ Actual Score, pvalue",pvalue),type="l")
      #lines(rep(actual,3),c(0,.001,.002),col="red",lwd=5)
      #cat ("Press [enter] to continue")
      #line = readline()
      pvalue_matrix[rr,tt] = pvalue  
    }
    # Now let's correct for multiple comparisons - we do for each row (region)
    
    qvalues[rr,] = p.adjust(pvalues,method="bonferroni")
  }
  
  # Now let's correct for multiple comparisons - we do for each row (region)
  # qvalues = apply(pvalue_matrix,1,p.adjust,method="fdr")

  # Plot p vs qvalues, all
  par(mfrow=c(1,2))
  hist(pvalue_matrix,main="pvalues, all regions and timepoints",col="red",xlab="pvalue",breaks=25)
  hist(qvalues,main="qvalues, all regions and timepoints",col="turquoise",xlab="qvalue",breaks=25)

  cat ("Press [enter] to continue")
  line = readline()
  
  # This will be a binary matrix of significant qvalues
  bin = qvalues
  bin[is.nan(bin)] = -99
  bin[bin>=bf] = -99
  bin[bin!=-99] = 1
  bin[bin==-99] = 0

  # Format sig matrix same as heatmap
  # Rearrange sig to be in same clustering as the heatmap
  bin = formatMatrix(bin)
  pheatmap(bin,cluster_rows=FALSE,cluster_cols=FALSE,main="Corrected temporal enrichment map")

  # Calculate probability scores for each timepoint
  probscores = colSums(bin) / nrow(bin)

  # Format the other matrices to return to user
  pvalue_matrix = formatMatrix(pvalue_matrix)
  qvalues = formatMatrix(qvalues)
  
  par(mfrow=c(1,1))
  image.plot(cbind(probscores,probscores),yaxt="n",xaxt="n",lab.breaks=timepoints,main="Temporal Probability Scores for Probe Set Influence",sub="Percentage of Brain Areas With Significant Enrichment")
  axis(1,labels=timepoints,seq(0,1,by =((1 - 0)/(length(timepoints) - 1))),pos=.1)
  cat ("Press [enter] to continue")
  line = readline()
  
  # Finally, filter original data to timepoints and regions with significant p values, plot in one place
  sigpairs = which(bin==1,arr.ind=TRUE)
  finallist = c()
  for (i in 1:nrow(sigpairs)){
    coord = sigpairs[i,]
    rg = rownames(sigpairs)[i]
    tp = coord[[2]]
    rowstoadd = intersect(which(genetimes$timepoint == timepoints[tp]),which(genetimes$region==rg))
    finallist = rbind(finallist,genetimes[rowstoadd,])
  }

  # Now plot!
  x = match(finallist$timepoint,colnames(interpolated$timeseries))
  plot(x,finallist$value,col=colors,pch=19,ylab="rna-seq expression",xaxt="n",main=paste("asd genes for significantly enriched timepoints and regions"),xlab="age",ylim=c(min(genetimes$value),max(genetimes$value)))
  # Add the x axis
  axis(1,labels=timepoints,at=seq(1,length(timepoints)))
  # Now let's add text for genes
  text(x,finallist$value,label=finallist$gene,cex=0.6,pos=4) 
  
  tmp = list(uncorrected=genetimes,corrected=finallist,pvalues=pvalue_matrix,qvalues=qvalues,binary=bin,timepoints=raw$timepoints,genes=raw$genes,regions=raw$regions,colors=colors)
  return(tmp)
}
