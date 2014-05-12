get.H.max.string.ver2<-function(H1) 
{ 
  # This function is used to get the methylation "string" : 0 0 0 1/2 1/2 1/2 1 1 1, and also output the maximum probability.
  # 0, Not methylated; 1/2, partly methylated; 1, fully methylated
  #
  # The output file is the H.summary (the input file) and "max.p" and "mCStatus". 
  # It has the following columns: "No.m", "Par.m", "Ful.m", "max.p", "mCstatus"
  # ("No.m", "Par.m", "Ful.m") mean "not methylated", "partly methylated" and "fully methylated". 
  # Note, the input H1 should be a N * 3 matrix, not 3*N matrix. 
  #
  # This function convert hidden states coding from 1,2,3 to 0,1/2,1.
  
   methy.status<-c(0, 1/2, 1) # their index (1,2,3) corresponds to (0, 1/2, 1) state respectively.
   H.max.mat<-apply(H1, 1, function(v) 
   {  max.1<-methy.status[v==max(v)]  
      if (length(max.1)==1) { return(c(max(v), max.1)) }   # max.1 gives the mCstatus that has the max probability
      else { if (length(max.1) >1) 
      {
        xx<-sample(1:length(max.1), 1)
        return(c(max(v), max.1[xx]))}                      # max.1[xx] gives one of the mCstatus that have equal max probability
      }
   } ) 


   H2<-cbind(H1, t( H.max.mat) )
   colnames(H2)<-c("No.m", "Par.m", "Ful.m", "max.p", "mCstatus")
   return(H2)
}
