DM.region.combine.HMM.fisher<-function(regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, p.threshold, report.singleCG, empty.CG)
{

   # This function is used to further combine small regions generated from chr.DM.region.by.CG.HMM.fisher
   # 1)regions: the output from chr.DM.region.by.CG.HMM.fisher
   # 2)chr.DM.status.matrix: the same input as in chr.DM.region.by.CG.HMM.fisher
   # 3)raw.CG, distance.threshold,  report.singleCG, empty.CG are the same as in chr.DM.region.by.CG.ver3
   # 4) num.CG.between: the max number of EM CG sites allowed between any two small DM regions
   # 5) p.threshold: the max p-value for the EM CGs inbetween
   # 6) report.singleCG: report the single DM CG or not.
   # 7) empty.CG: a numeric shows the threshold of number of CGs without coverage between consecutive CG sites to combine together 
   # output:
   # The output file has 10 columns, for example: 
   #  chr     start     end       len  DM      num.CG total.CG     meanCov.control meanCov.test  meanDiff 


   ####################
   #  Get hyper region 
   ####################
   hyper.combine<-DMR.combine( DM.type="hyper", regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, p.threshold, empty.CG)
   ####################
   #  Get hypo region 
   ####################
   hypo.combine<-DMR.combine( DM.type="hypo", regions, chr.DM.status.matrix, raw.CG, distance.threshold, num.CG.between, p.threshold, empty.CG)

   combined<-rbind(hyper.combine, hypo.combine)
   combined<-as.data.frame(combined, stringsAsFactors = FALSE)
   combined[,c(2:4,6:10)]<-apply(combined[,c(2:4,6:10)], 2, function(x) as.numeric(x))


   row.names(combined)<-NULL
   if (report.singleCG==TRUE)
   {  return(combined)  }
   else 
   { 
      singular.index<-(1:dim(combined)[1])[as.numeric(combined[,4])==1]
      combined.NOsingular<-combined[-singular.index, ]
      return(combined.NOsingular)
   }
} 
