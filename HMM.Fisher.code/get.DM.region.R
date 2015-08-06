get.DM.region<-function(DM.type,chr.DM.status.matrix, raw.CG, chr, distance.threshold, empty.CG)
{
   # This function is used to get hyper or hypo DM regions. It is called by chr.DM.region.by.CG.HMM.fisher()
   # The basic idea of this function is: 
   # First,get the index (i.e., 1, 2, .., n) and then group those "indexs" into regions based on their distance. 
   #
   # Note 1: DM.type, the type of DM regions to summarize, 1 (hyper) or -1 (hypo).
   # Note 2: chr.DM.status.matrix is the result from HMM.fisher. It has 22 columns:
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
   # Note 3: "raw.CG", is a vector of all CG positions on that chr 
   # Note 4: "chr", a numeric value shows the chr number 
   # Note 5: "distance.threshold", a numeric value shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # Note 6: "empty.CG", a numeric value shows the threshold of number of CGs without coverage between consecutive CG sites to combine together.
   #
   # The output file has 10 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG    ER+.ave.cov ER-.ave.cov   meandiff
   # "chr1" "4252119" "4252128" "10" "DM" "2"       "2"             "49.5"      "97.75"      
    
   chr.DM.region.mat<-matrix(NA, nrow=0, ncol=10)
   colnames(chr.DM.region.mat)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "pos.ave.cov", "neg.ave.cov", "meanDiff")

   DM.index<-(1:dim(chr.DM.status.matrix)[1])[as.character(chr.DM.status.matrix[, 19])== as.character(DM.type)]
   
   ####################
   #  Get DM region 
   ####################
   
   if ( length(DM.index)<=1) 
   {  cat("There are:", length(DM.index), " ", DM.type, " DM CG sites, we do not need to summarize \n") }  
   if ( length(DM.index)>1) 
   {  
      # get the index for DM CG sites (using the 20th column of chr.DM.status.matrix)
      DM.index.dis<- chr.DM.status.matrix[DM.index[-1],20] - chr.DM.status.matrix[DM.index[-length(DM.index)],20] 
      physical.dis<-chr.DM.status.matrix[DM.index[-1],2] - chr.DM.status.matrix[DM.index[-length(DM.index)],2]
      start.index<-DM.index[1]; end.index<-DM.index[1];
      empty<-rep(0, (length(DM.index)-1))
      for (j in 1:(length(DM.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(chr.DM.status.matrix[DM.index[j],2]) & as.numeric(raw.CG)< as.numeric(chr.DM.status.matrix[DM.index[j+1],2]) ])
	}

      i<-1
      while(i<=(length(DM.index)))
      {   # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
        
         if ( i < length(DM.index) )
         {   if ( DM.index.dis[i]==1 && physical.dis[i]<= distance.threshold && empty[i]<= empty.CG )
             { 
                # That is, the two consecutive index are two successive number, that is, x and x+1, and the physical distance between them is <= distance, and the number of CGs without coverage <=empty.CG 
                # We should add more to this region.
                end.index<-DM.index[i+1]
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
                i<-i+1
             }
             else 
             { 
                # That is, the two consecutive index are NOT two successive number, that is, they are NOT x and x+1, rather than x and x+n. Or the distance is > distance.
                # We should end this region and start a new one.
                chrom<-paste("chr",chr,sep="")  ; start<-chr.DM.status.matrix[start.index, 2]
                end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
                DM<-DM.type
                num.CG<- end.index - start.index + 1 
                group1.ave.cov<-round (mean(chr.DM.status.matrix[start.index:end.index, 21]),2)
                group2.ave.cov<-round( mean(chr.DM.status.matrix[start.index:end.index, 22]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.index:end.index, 17]),4)
                # smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])


                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group2.ave.cov, meanDiff.mC)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
                i<-i+1
                if ( i<=length(DM.index) ) {  start.index<-DM.index[i]; end.index<-DM.index[i] }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the DM index if ( i==length(DM.index) )
           # here we print out everything we got:

              chrom<-paste("chr",chr,sep="")    ; start<-chr.DM.status.matrix[start.index, 2]
              end<- chr.DM.status.matrix[end.index, 2]; leng<-end-start +1
              DM<-DM.type
              # Then if there are many "0, 1", we will report 1, that is, as long as there are at least one 
              num.CG<- end.index - start.index + 1 
                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.index:end.index, 21]),2)
                group2.ave.cov<round(-mean(chr.DM.status.matrix[start.index:end.index, 22]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.index:end.index, 17]),4)
               #  smooth.meanDiff<-mean(chr.DM.status.matrix[start.index:end.index, 13])


		total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions


              new.vec<-c(chrom, start, end, leng, DM,  num.CG, total.CG,  group1.ave.cov, group2.ave.cov, meanDiff.mC)
              chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
              # cat("when i is:", i, "start and end index:", c(start.index, end.index), "start and end pos:", c(start, end), "\n") 
             
              break
         } 
      }
   }
  return (chr.DM.region.mat)
}
