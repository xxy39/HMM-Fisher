DMR.combine<-function( DM.type, regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, p.threshold, empty.CG)
{

   # This function is used to further combine small regions generated from chr.DM.region.by.CG.HMM.fisher. It is called by DM.region.combine.HMM.fisher().
   # 1) DM.type, the type of DM regions to summarize, "DM" or "hypo"
   # 2) regions: the output from chr.DM.region.by.CG.HMM.fisher
   # 3) chr.DM.status.matrix: the same input as in chr.DM.region.by.CG.HMM.fisher
   # 4) raw.CG, 
   # 5) distance.threshold: a numeric shows the threshold of physical distance. The CG sites with distance larger than this value won't be in the same region.
   # 6) num.CG.between: the max number of EM CG sites allowed between any two small DM regions
   # 7) posterior.threshold: the max posterio probability for the EM CGs inbetween
   # 8) empty.CG: a numeric shows the threshold of number of CGs without coverage between consecutive CG sites to combine together 
   # The output file has 10 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff 

   chr.DM.region.mat<-matrix(NA, nrow=0, ncol=10)
   colnames(chr.DM.region.mat)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "meanCov.control", "meanCov.test", "meanDiff")

   DM.index<-(1:dim(regions)[1])[as.character(regions[, 5])== as.character(DM.type)]
   
   if ( length(DM.index)<=1) 
   {  cat("There are:", length(DM.index), " ", DM.type," CG regions, we do not need to summarize \n") }  
   if ( length(DM.index)>1) 
   {       
      # get the distance between DM regions
      physical.dis<-regions[DM.index[-1],2] - regions[DM.index[-length(DM.index)],3]
      start.index<-1; end.index<-1;
      empty<-rep(0, (length(DM.index)-1))
      for (j in 1:(length(DM.index)-1))
        {
           empty[j]<-length((1:length(raw.CG))[as.numeric(raw.CG)> as.numeric(regions[DM.index[j],3]) & as.numeric(raw.CG)< as.numeric(regions[DM.index[j+1],2]) ])
	}

      i<-1
      while(i<=(length(DM.index)))
      {   # cat("when i is:", i, "start and end index:", c(start.index, end.index), "\n") 
         
	 add<-0
         if ( i < length(DM.index) )
         {   if (  physical.dis[i]<= distance.threshold)
             { 
                # the physical distance between them is <= distance

                # get the other CG sites inbetween
                between.CG.index<-(1:dim(chr.DM.status.matrix)[1])[chr.DM.status.matrix[,2] > regions[DM.index[i],3] & chr.DM.status.matrix[,2] < regions[DM.index[i+1],2] ]
              
                # if both regions have 1 DM CG site, only allow 1 CG site in between
		 if (regions[DM.index[i],6]==1 && regions[DM.index[i+1],6] == 1 )
		   { num.CG.between=1 }

		 if (length(between.CG.index)>0 && length(between.CG.index) <= num.CG.between && empty[i]<= empty.CG)
		     {
                              between.CG.state<-chr.DM.status.matrix[between.CG.index,19]
			      between.CG.raw.state<-chr.DM.status.matrix[between.CG.index,18]
                              between.CG.p<-chr.DM.status.matrix[between.CG.index,4]

                              vec<-rep(0,length(between.CG.index))
			      for (j in 1:length(between.CG.index))
			        {
                                     if (between.CG.state[j]=="EM" && between.CG.p[j]<=p.threshold || between.CG.raw.state[j]==as.character(DM.type) )
				      { vec[j]<-1}				     
				}

			      if (  sum(vec)== length(between.CG.index)  )
				       { add<-1 }
		     }			       
             }

             if (add==1)
	       {end.index<-i+1
	        i<-i+1 }
             else 
             { 
                chrom<-regions[1,1]  ; start<-regions[DM.index[start.index], 2]
                end<- regions[DM.index[end.index], 3]; leng<-regions[DM.index[end.index], 3] - regions[DM.index[start.index],2] + 1
                DM<-DM.type
                num.CG<- sum(regions[DM.index[start.index: end.index], 6])

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 20]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end,  20]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 21]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 22]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 17]),4)

                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions

                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group2.ave.cov, meanDiff.mC)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)
                i<-i+1
                if ( i<=length(DM.index) ) {  start.index<-i; end.index<-i }  # start a new search
            } 
        }

        else
        {  # Here is the case that we reach the end of the DM index if ( i==length(DM.index) )
                chrom<-regions[1,1]  ; start<-regions[DM.index[start.index], 2]
                end<- regions[DM.index[end.index], 3]; leng<-regions[DM.index[end.index], 3] - regions[DM.index[start.index],2] + 1
                DM<-DM.type
                num.CG<- sum(regions[DM.index[start.index: end.index], 6])

		start.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== start, 20]
		end.pos.index<- chr.DM.status.matrix[chr.DM.status.matrix[,2]== end, 20]

                group1.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 21]),2)
                group2.ave.cov<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 22]),2)
                meanDiff.mC<-round(mean(chr.DM.status.matrix[start.pos.index:end.pos.index, 17]),4)

                total.CG.index<-(1:length(raw.CG))[as.numeric(raw.CG)<=as.numeric(end) & as.numeric(raw.CG)>=as.numeric(start)]
		total.CG<-length(total.CG.index)  # total number of CG sits within this DM regions

                new.vec<-c(chrom, start, end, leng, DM, num.CG, total.CG, group1.ave.cov, group2.ave.cov, meanDiff.mC)
                chr.DM.region.mat<-rbind(chr.DM.region.mat, new.vec)             
                break
         } 
      }
   }
  return (chr.DM.region.mat)
}
