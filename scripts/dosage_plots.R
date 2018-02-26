# Careful, running this requires a local install of R packages on any cluster you intend to work on. I did not forsee this.
# The snakemake approach: specify appropriate versions of R and associated packages in environment.yaml, assmuing the conda versions of these packages are compatible

# USAGE: RScript dosage_plots.R bedfiles.txt <controlfile>

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
# setwd("") # only if needed

header.temp <- c("chrom", "start", "end", "readcount", "bases_covered", "binsize", "breadth")
filenames <- read.csv(args[1],header=F) # hardcoded path!
input <- apply(filenames, 1, function(x) read.csv(x, sep='\t', comment.char="#", header=F, col.names=header.temp)) # need generic column names and specific list names
# potential memory fault, with a very large number of bed files. each bed file is ~ dozens of KB.

names(input) <- gsub("_coverage.bed", "",filenames$V1) # hardcoded, but should be ok as long as snakemake pipeline naming is maintained.

# specify control file here, need a way to make this more command line friendly
control <- read.csv(args[2], sep='\t', comment.char="#", header=F, col.names=header.temp) # idea, snakemake this pipeline and use config.yaml to specify control file

# If statement catches first data frame only, should be sufficient.
if ("normcov" %in% colnames(input[[1]])) {} else {
  input <- lapply(input, function(x) cbind(x, "normcov"=2*(x$readcount/sum(x$readcount)) / (control$readcount / sum(control$readcount)))) # did manual math on representative bins, came out fine
} 

# plot dosage for each entry of data frame.

mid.dosage <- function(x,y) { # x, y is title of plot
  n=floor(nrow(x)/2)
  return(x[n,9])
}

plot.dosage <- function(x,y) {
  numblanks.dosage <- 15 # is aesthetically pleasing for 1Mb bins, but can customize to accommodate any bin size or number of chromosomes
  stuf.d <- c(rep(NA, numblanks.dosage))
  dosage.stuffer <- data.frame("chrom"=stuf.d, "start"=stuf.d, "end"=stuf.d,
                               "readcount"=stuf.d, "bases_covered"=stuf.d,
                               'binsize'=stuf.d,"breadth"=stuf.d,"normcov"=stuf.d,
                               "bin"=stuf.d)
  x$bin <- seq(1,nrow(x))
  chr.list.dosage <- split(x, f=x$chrom)
  chr.list.dosage.stuffed <- lapply(chr.list.dosage[2:12], function(x) rbind(x,dosage.stuffer))
  chr.list.dosage.stuffed <- dplyr::bind_rows(chr.list.dosage.stuffed, chr.list.dosage[13])
  chr.list.dosage.stuffed$bin2 <- seq(1:nrow(chr.list.dosage.stuffed))
  chr.list.dosage.stuffed$normcov <- as.numeric(chr.list.dosage.stuffed$normcov) # force as numeric if not
  
  midpoints.dosage <- sapply(chr.list.dosage[2:13],mid.dosage)
  
  plt.dosage <- ggplot(chr.list.dosage.stuffed, aes(x=bin2,y=normcov)) +
    labs(x="",y="Copy Number") +
    geom_line(color="#008080",fill="white",size=1) + # color currently hardcoded
    ggtitle(names(y)) +
    geom_point(aes(x=bin2,y=normcov), size=1.0, color="black") +
    guides(fill=F,color=F) +
    scale_x_continuous(breaks=which(chr.list.dosage.stuffed$bin %in% midpoints.dosage),
                       labels=names(midpoints.dosage)) +
    scale_y_continuous(limits=c(0.75,4.25), breaks=seq(1,4), labels=seq(1,4)) +
    theme(panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(color="black",linetype="dashed"),
          panel.background=element_rect(fill="white",color="black"),
          axis.text.x=element_text(size=12,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          axis.title.y=element_text(size=18,angle=90,vjust=0.5),
          axis.ticks=element_blank(),
          plot.title=element_text(size=24,face="bold",hjust=0))
  write.table(chr.list.dosage.stuffed, paste(names(y),"_plot_table.tsv",sep=""), quote=F,eol='\n',
              row.names=F)
  tosave <- paste(names(y),"_dosage_plot.pdf",sep="")
  ggsave(tosave, width=10,height=6,units="in",plot=plt.dosage,device="pdf")
}

# Loop for making plots
for (i in seq(1:length(input))) {
  plot.dosage(input[[i]],input[i])
}