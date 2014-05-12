getMeth<-function(total.reads, meth.reads, n1, n2, min.percent)
  {
       # This function is used to get methylation ratio retults for CG sites that pass the quality control
       #
       # 1) totoal.reads: a matrix of the total reads covering each CG sites. 
       #    Columns: the position, the coverage for control sampels,  the coverage for test sampels
       # 2) meth.reads: a matrix of the methylated reads covering each CG sites. 
       #    Columns: the position, the coverage for control sampels,  the coverage for test sampels
       # 3) n1: number of control samples
       # 4) n2: number of test samples
       # 5) min.percent: the min percentage of samples with coverage in each group
       #
       # 6) This funtion RETURN a list with 2 elements:
       #     (1) element 1: the methylation ratio for each CG that passes the quality control
       #     (2) element 2: the index of CGs that pass the quality control
       
       meth.ratio<-meth.reads[,-1]/total.reads[,-1]
       meth.ratio[is.na(meth.ratio)]<-NA

       con.thres<-ceiling(n1 * min.percent)
       test.thres<-ceiling(n2 * min.percent)

       QC.index<- (1:dim(total.reads)[1])[ (apply(total.reads[,-1], 1, function(x) { return( sum((x[1:n1]>0)*1)>=con.thres & sum((x[(n1+1):(n1+n2)] >0)*1)>=test.thres)  }))*1 ==1]
       
       methy.level<-cbind(total.reads[,1], meth.ratio)
       QC.methy<-methy.level[QC.index,]
       return(list(QC.methy,QC.index))
  }

