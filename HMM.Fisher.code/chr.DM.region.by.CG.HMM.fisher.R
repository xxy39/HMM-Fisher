chr.DM.region.by.CG.HMM.fisher<-function(chr.DM.status.matrix, raw.CG, chr, distance.threshold, report.singleCG, empty.CG)
{
   # This function is used to get large regions that is either hyper or hypo methylated from a whole chromosome. 
   # The basic idea of this function is: 
   # First, get the index (i.e., 1, 2, .., n) and then group those "indexs" into regions based on their distance. 
   # ---------------------------------------------------------------------------
   # Note 1: in "chr.DM.status.matrix", usually results from HMM.fisher. It has the following 22 columns: 
   #    column 1: chr
   #    column 2: position of ith CG
   #    column 3: position of (i+1)th CG. (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 4: p-value of fisher-exact test
   #    column 5: estimated states for ith CG sites in control samples, seperated by ":"
   #    column 6: estimated states for ith CG sites in test samples, seperated by ":"
   #    column 7: estimated states for (i+1)th CG sites in control samples, seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 8: estimated states for (i+1)th CG sites in test sampels seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 9: raw mC levels for ith CG sites in control samples, seperated by ":"
   #    column 10: raw mC levels for ith CG sites in test samples, seperated by ":"
   #    column 11: raw mC levels for (i+1)th CG sites in control samples, seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 12: raw mC levels for (i+1)th CG sites in test samples, seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 13: posterior prob for ith CG sites in control samples, seperated by ":"
   #    column 14: posterior prob for ith CG sites in test samples, seperated by ":"
   #    column 15: posterior prob for (i+1)th CG sites in control samples, seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 16: posterior prob for (i+1)th CG sites in test samples, seperated by ":". (if ith CG is not combined with (i+1)th CG, this column is "NA" )
   #    column 17: meanDiff: (mean mC of control) - (mean mC of test)
   #    column 18: raw DM state: 1: hyper (ER+ has higher mC); 0: EM (ER+ and ER- have similar mC); -1: hypo (ER- has higher mC).   
   #    column 19: DM state (with meanDiff): 1: hyper (control has higher mC); 0: EM (two groups have similar mC); -1: hypo (test has higher mC).   
   #    column 20: index of CG sites
   #    column 21: mean coverage of control
   #    column 22: mean coverage of test
   #---------------------------------------------------------------------------
   # Note 2: "raw.CG", is a vector of all CG positions on that chr 
   # ---------------------------------------------------------------------------
   # Note 3: "chr", a numeric shows the chr number 
   # ---------------------------------------------------------------------------
   # Note 4: "distance", a numeric shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # Note 5: "empty.CG", a numeric shows the threshold of number of CGs without coverage between consecutive CG sites to combine together.
   # Note 6: The summary output file include a lot of singlar points and also some DM region with low coverage, 
   #         Here we use a flag: report.singleCG=FALSE to remove it at the end. 
   # ---------------------------------------------------------------------------
   # The output file has 10 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG    ER+.ave.cov ER-.ave.cov   meandiff
   # "chr1" "4252119" "4252128" "10" "hyper" "2"       "2"             "49.5"      "97.75"      
   # "chr1" "4252143" "4252143" "1"  "hyper" "1"       "1"            "46.5"      "90.5"       
   # "chr1" "4248700" "4248702" "3"  "hypo"  "2"       "2"           "7.125"     "13.25"       
 

   chr.DM.status.matrix<-chr.DM.status.matrix[order(chr.DM.status.matrix[,2]),]

   # get hyper regions
   hyper.regions<-get.DM.region(DM.type="hyper",chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)

   # get hypo regions
   hypo.regions<-get.DM.region(DM.type="hypo",chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)

   DM.regions<-rbind(hyper.regions, hypo.regions)
   DM.regions<-as.data.frame(DM.regions, stringsAsFactors = FALSE)
   DM.regions[,c(2:4,6:10)]<-apply(DM.regions[,c(2:4,6:10)], 2, function(x) as.numeric(x))

   row.names(DM.regions)<-NULL
   if (report.singleCG==TRUE)
    { return(DM.regions) }
   else 
    { 
      singular.index<-(1:dim(DM.regions)[1])[as.numeric(DM.regions[,4])==1]
      DM.regions.NOsingular<-DM.regions[-singular.index, ]
      return(DM.regions.mat.NOsingular)
    }
} 
