#' probePeaks
#'
#' This function will extract temporal features for a gene set of interest to predict times of temporal enrichment.
#' @param genelist A list of gene symbol identifiers
#' @param thresh The quantile of expression values to take per gene.  Default is 0.95
#' @param out.pdf output pdf of temporal enrichment patterns, defaults to FALSE.  Default output name (set with out.name) is "genetime.pdf"
#' @param out.name name of output pdf of temporal enrichment patterns, default is "genetime.pdf"
#' @keywords temporal gene enrichment
#' @export
#' @examples
#' mygenes = "ADNP" "AGAP2" "ANK2" "APH1A"
#' thresh = 0.99
#' probePeaks(genelist=mygenes,threshold=thresh,out.pdf=TRUE,out.name="/home/myfile.pdf")

probePeaks = function(genelist,thresh=0.95,out.pdf=FALSE,out.name="genetime.pdf"){
  
  # This is so legend can print outside plot area
  par(xpd=TRUE)
  
  # Get colors for the brain regions
  colors = getColors()
      
  # A helper function to return threshold for each row
  get_top_quantile = function(row,thresh) {
    qpos = quantile(row,thresh)
    return(qpos)
  }

  # Get genes in developmental enrichment database
  devgenes = probeLookup(genelist)
  
  # Exit if there is no gene overlap
  stopifnot(length(devgenes) != 0)
  
  # Load temporal enrichment interpolated timeseries
  interpolated = getData()
  
  # Here are the x axis labels
  xlabels = colnames(interpolated$timeseries)
  
  # If the user wants to print genes to file
  if (out.pdf==TRUE){
    wd = getwd()
    outfile = paste(wd,"/",out.name,sep="")
    pdf(outfile,onefile=TRUE,width=12)
  }
  
  # We will save data (genes,times,regions,expression) to "genetimes"
  genetimes = c()
  
  # Iterate through genes
  for (gene in devgenes) {
  
    # Get the subset of the data
    subset = interpolated$timeseries[which(interpolated$genes==gene),]
    regions = gsub(paste(gene,"_",sep=""),"",rownames(subset))
    
    # If out.pdf is true, make plots of each gene
    if (out.pdf == TRUE){ 
    
      # Print temporal plot of each gene to pdf file  
      plot(subset[1,],col=colors[1],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("rna-seq expression for ",gene),xlab="age",ylim=c(min(subset),max(subset)))
      axis(1,labels=colnames(subset),at=seq(1,length(colnames(subset)))) 
      for (o in 2:nrow(subset)) {
        points(subset[o,], col = colors[o],pch=19) 
      }
      legend(25,max(subset),regions,col=colors,cex=.6,pch=19)
  
    }
    
    # Method: Max values in top thresh (quartile) of distribution
    
    # Try getting quantiles for top values
    quantiles = apply(subset,1,get_top_quantile,thresh)
  
    # If any quantiles are == 0, set to impossibly large number
    quantiles[which(quantiles == 0)] = 99999
    
    # Get values above threshold for each row
    maxvalues = c()
    timepoints = c()
    reg = c()
    for (r in 1:nrow(subset)) {
      region = r
      idx = which(subset[r,]>=quantiles[r])
      maxvalues = c(maxvalues,subset[r,idx])
      timepoints = c(timepoints,colnames(subset)[idx])
      reg = c(reg,rep(region,length(idx)))
    }

    # Now let's make into a data frame, save to result
    result = data.frame(gene=rep(gene,length(maxvalues)),timepoint = timepoints,region=regions[reg],value=maxvalues,thresh=rep(thresh,length(maxvalues)))    
    genetimes = rbind(genetimes,result)

    # What are the unique times?
    uniquetimes = unique(timepoints)
    
    # If the user wants to plot
    if (out.pdf==TRUE){   
      time = uniquetimes[1]
      x = which(xlabels == time)
      values =  maxvalues[which(timepoints == time)]  
      x = rep(x,length(values))
      r = reg[which(timepoints==time)]
      plot(x,as.numeric(values),col=colors[r],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("timepoints of top",1-thresh,"of values for each region for",gene),xlab="age",ylim=c(min(subset),max(subset)),xlim=c(1,ncol(subset)))
      axis(1,labels=colnames(subset),at=seq(1,length(colnames(subset)))) 
  
      # Add the other points
      for (m in 2:length(uniquetimes)) {
        # Get the x coordinate
        time = uniquetimes[m]
        x = which(xlabels == time)
        # The color is just the value of count, corresponding to row (region)
        values =  maxvalues[which(timepoints == time)]
        x = rep(x,length(values))
        r = reg[which(timepoints==time)]
        points(x,as.numeric(values), col = colors[r],pch=19)
      }
      legend(25,max(subset)/1.5,pch=19,regions,col=colors,cex=.6)
    }  
  }

  # Prepare final output object  
  genetimes = list(matrix=genetimes,timepoints=xlabels,regions=regions,genes=devgenes)
  
  # Now make plots by region
  if (out.pdf==TRUE){
    for (r in 1:length(regions)) {
      
      region = regions[r]
      subset = genetimes[which(genetimes$region == region),]
      regioncolor = colors[r]
      
      # First plot all the genes for a single region
      idx = grep(region,rownames(interpolated$timeseries))
      singleregion = interpolated$timeseries[idx,]
      
      # Plot the gene times in the region
      # Here are the x values
      x = match(subset$timepoint,colnames(interpolated$timeseries))
      plot(x,subset$value,col=regioncolor,pch=19,ylab="rna-seq expression",xaxt="n",main=paste("important timepoints and genes for",region),xlab="age",ylim=c(min(subset$value),max(subset$value)))
      axis(1,labels=colnames(interpolated$timeseries),at=seq(1,length(labels)))
      # Now let's add text for genes
      text(x,subset$value,label=subset$gene,cex=0.6,pos=4) 
    }
    dev.off()
  }
  return(genetimes)
}