# Genetime
Predict "important" time periods for a list of gene probes based on RNA-seq expression.

###Usage

genetime(genelist, thresh = 0.95, out.pdf = FALSE, out.name = "genetime.pdf", iters = 1000)

###Arguments

genelist	
A list of gene symbol identifiers

out.pdf	
output pdf of temporal enrichment patterns, defaults to FALSE. Default output name (set with out.name) is "genetime.pdf"

out.name	
name of output pdf of temporal enrichment patterns, default is "genetime.pdf"

iters	
Number of iterations for permutations to create null distributions, Default is 1000

threshold	
The quantile of expression values to take per gene. Default is 0.95

###Examples

mygenes = "ADNP" "AGAP2" "ANK2" "APH1A"
thresh = 0.99
genelist(genelist=mygenes,threshold=thresh,out.pdf=TRUE,out.name="/home/myfile.pdf")

###Installation
Since there is a data file over 100MB, it cannot be stored on github, so you must download them separately.  First clone the repository
  
  git clone https://github.com/vsoch/genetime
  cd genetime/data
  
Download the data file

  wget http://www.vbmis.com/bmi/data/interpolated_timeseries.rda
  
Change to the directory above the "genetime" folder, start R, and install from source
  
  R
  install.packages("genetime",rep=NULL,type="source")
  library("genetime")


