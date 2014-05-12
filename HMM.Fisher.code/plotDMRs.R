# plotDMRs.R
#
# R script to plot specified DMRs
#
# Usage: 
# R CMD BATCH '--args Input1 Input2  index  extend  test  control header output'   HMM.Fisher.code/plotDMR.R
#
############################################################################################################
# 1) Input1: the mC.matrix.txt results from HMM-Fisher program
# pos     test_1     test_2      test_3       test_4  control_1  control_2  control_3  control_4       
# 497  0.988701   0.886364    0.886598     0.977778   0.602339   0.956522   0.979899   0.936508   
# 525  0.971591   1.000000    0.964286     0.956522   0.970930   0.949640   0.959799   0.984375   
# 542  0.944056   1.000000    0.978495     0.932584   0.942149   0.992647   0.946524   0.909091   
###########################################################################################################
#
############################################################################################################
# 2) Input2: the DMRs.txt results from HMM-Fisher program
# chr     start   end     len     DM      num.CG  total.CG        meanCov.test    meanCov.control meanDiff
# chr1    795361  795361  1       hyper   1       1       69      70.25   0.4652
# chr1    848868  848873  6       hyper   2       2       25.12   30.62   0.4272
# chr1    851081  851123  43      hyper   5       5       5.45    8.1     0.6126
# chr1    858338  858368  31      hyper   2       2       22.25   12      0.5629
###########################################################################################################
#
# 3) index: a vector of the index indicating which DMRs in DMR.txt file that the user what to plot.
#           e.g, c(1,4:6) means to plot the 1st, and 4th to 6th DMRs in the DMR.txt file.
# 4) extend: number of bps extended from both ends of the region to plot
# 5) test: number of test samples
# 6) control: number of control samples
# 7) header: whether the input2(DMR.txt) has header line, T or F. Note: The DMRs.txt from HMM-DM program has a header line.
# 8) outpput: the name of output .ps file.
#
# example command line: R CMD BATCH '--args   mC.matrix.txt   DMRs.txt  c(26:27:59)  100   4  4  T  example.DMR.plot'   HMM.Fisher.code/plotDMRs.R

args<-commandArgs(trailingOnly = TRUE);
headerL<-as.logical(args[7]);
# read in the mC.matrix.txt and DMRs.txt files
mC.input <- read.table(file=args[1], header=F);
DMR.input <- read.table(file=args[2], header=headerL); 

# define the other parameters
index <- eval( parse(text=args[3]) );
extend <-as.integer(args[4]);
control <-as.integer(args[5]);
test <-as.integer(args[6]);
output <-as.character(args[8]);

output.title <- paste(output, ".ps",sep="");

# plot the region
postscript(output.title,horizontal=T, paper="letter")
par(mfrow=c(1,1))
  for (i in index)
  {   
    # define the start and end position of the region
    start <- DMR.input[i,2]
    end <- DMR.input[i,3]
    
    # define the region after extension
    plot.start <- start-extend
    plot.end <- end+extend
    region <- mC.input[mC.input[,1]>=plot.start & mC.input[,1]<=plot.end,]
    
    title <- paste("DMR ", DMR.input[i,1], " ",start, ":", end, sep="")
    plot( c(plot.start, plot.end), c(0,1), type="n", xlab="Position",  ylab="Methylation Level",  main=title,cex.main=1.5, cex.axis=1.5,cex.lab=1.5)
    rect(xleft=start-1, xright=end+1, ybottom=-0.02, ytop=1.02,  density=NA, col="lightgray")
    for (j in 1:test)
     { points(region[,1],region[,1+j], pch=15, col="red", cex=1.5)}
    for (j in (test+1):(control+test))
     { points(region[,1],region[,1+j], pch=17, col="blue", cex=1.5)}
    legend ("topleft", c("test", "control"), pch=c(15,17), , cex=1.5, col=c("red", "blue"))
  }
dev.off()


