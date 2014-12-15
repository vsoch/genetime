#' genetime
#'
#' This function will predict "important" time periods for a list of gene probes based on RNA-seq expression.
#' @param genelist A list of gene symbol identifiers
#' @param threshold The quantile of expression values to take per gene.  Default is 0.95
#' @param iters Number of iterations for permutations to create null distributions,  Default is 1000
#' @param out.pdf output pdf of temporal enrichment patterns, defaults to FALSE.  Default output name (set with out.name) is "genetime.pdf"
#' @param out.name name of output pdf of temporal enrichment patterns, default is "genetime.pdf"
#' @keywords temporal gene enrichment
#' @export
#' @examples
#' mygenes = "ADNP" "AGAP2" "ANK2" "APH1A"
#' thresh = 0.99
#' genelist(genelist=mygenes,threshold=thresh,out.pdf=TRUE,out.name="/home/myfile.pdf")

genetime = function(genelist,thresh=0.95,out.pdf=FALSE,out.name="genetime.pdf",iters=1000){

  library(ade4)
  library(pheatmap)
  library(fields)
  
  # So we can plot outside of main plot area
  par(xpd=TRUE)
  
  # First get raw times, regions, and genes for the set
  raw = probePeaks(devgenes,thresh=thresh,out.pdf=out.pdf,out.name=out.name)

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

  cat("Running permutations...\n")
  
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
    
    # Do shapiro test to test for normality
    #sha = apply(nulldist, 2, shapiro.test)
    #nor=0
    #for (i in sha){
    #  if ((i$p.value>0.05)==T){nor=nor+1}
    #}    
  
    # Save pvalues in a vector to correct for 16 tests
    for (tt in 1:length(timepoints)){
      single_null_dist = nulldist[,tt]
      confidence_intervals = mean(single_null_dist) + c(-1,1)*qnorm(0.95)*sd(single_null_dist)
      #hist(single_null_dist,col="purple",main=paste("Null Distribution timepoint",timepoints[tt]))
      #lines(rep(confidence_intervals[1],100),seq(1,100),col="red",lwd=2)
      #lines(rep(confidence_intervals[2],100),seq(1,100),col="red",lwd=2)
      # Calculate Z score for our timepoint score
      score = countmatrix[rr,tt]
      z = (mean(single_null_dist)-score)/sd(single_null_dist)
      # Get p value, single sided Z test
      pvalue = 1-(pnorm(-abs(z),lower.tail=FALSE))
      pvalue_matrix[rr,tt] = pvalue
    }
  }
  
  cat("Correcting...\n")  
  # Now let's correct for multiple comparisons
  qvalues = p.adjust(pvalue_matrix,method="fdr")

  # Plot p vs qvalues, all
  par(mfrow=c(1,2))
  hist(pvalue_matrix,main="pvalues, all regions and timepoints",col="red",xlab="pvalue",breaks=25)
  hist(qvalues,main="qvalues, all regions and timepoints",col="turquoise",xlab="qvalue",breaks=25)

  cat ("Press [enter] to continue")
  line = readline()
  
  # Which ones are significant?
  qvalues = array(qvalues,dim=c(length(uniqueregions),length(timepoints)))

  # This will be a binary matrix
  bin = qvalues
  bin[is.nan(bin)] = -99
  bin[bin>=0.05] = -99
  bin[bin!=-99] = 1
  bin[bin==-99] = 0

  # Format sig matrix same as heatmap
  # Rearrange sig to be in same clustering as the heatmap
  bin = formatMatrix(bin)
  pheatmap(bin,cluster_rows=FALSE,cluster_cols=FALSE,main="Corrected temporal enrichment map")

  # Calculate probability scores for each timepoint
  probscores = colSums(sig) / nrow(sig)

  # Format the other matrices to return to user
  pvalue_matrix = formatMatrix(pvalue_matrix)
  qvalues = formatMatrix(qvalues)
  
  par(mfrow=c(1,1))
  image.plot(cbind(probscores,probscores),yaxt="n",xaxt="n",lab.breaks=timepoints,main="Temporal probability scores for probes being active")
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
  plot(x,finallist$value,col=regioncolor,pch=19,ylab="rna-seq expression",xaxt="n",main=paste("asd genes for significantly enriched timepoints and regions"),xlab="age",ylim=c(min(genetimes$value),max(genetimes$value)))
  # Add the x axis
  axis(1,labels=timepoints,at=seq(1,length(timepoints)))
  # Now let's add text for genes
  text(x,finallist$value,label=finallist$gene,cex=0.6,pos=4) 
  
  tmp = list(uncorrected=genetimes,corrected=finallist,pvalues=pvalue_matrix,qvalues=qvalues,binary=bin,timepoints=raw$timepoints,genes=raw$genes,regions=raw$regions,colors=colors)
  return(tmp)
}