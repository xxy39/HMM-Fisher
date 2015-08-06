fisher.exact<-function(HMM.results, mC.matrix, distance.status, n1, n2)
  {
        #####################################################################################################################
	# Fisher's exact test step for HMM-Fisher
	#
	# 1. HMM.results: resutls from HMM step of HMM-Fisher.
	# 2. mC.matrix: methylation levels of all samples. This file is the output from function getMeth()
	# 3. distance.status: output from function get.distance.status.
        # 4. n1:  Numeric. number of testsamples
	# 5. n2:  Numeric. number of control samples
	#
	# output: matrix with 15 columns
	#         1-2: position for current and next CG site
	#         3: p-value of Fisher's exact test
	#         4: the states of current CG estimated by HMM for control samples. 0, Not methylaited; 0.5, Partly methylated; 1, fully methylated
	#         5: the states of current CG estimated by HMM for test samples. 
	#         6: the states of the next CG estimated by HMM for control samples
	#         7: the states of the next CG estimated by HMM for test samples
	#         8: the raw methylation level of current CG in control samples
	#         9: the raw methylation level of current CG in test samples
	#         10: the raw methylation level of next CG in control samples
	#         11: the raw methylation level of next CG in test samples
	#         12: the posterior probabilities of current CG in control samples
	#         13: the posterior probabilities of current CG in test samples
	#         14: the posterior probabilities of next CG in control samples
	#         15: the posterior probabilities of next CG in test samples
        #####################################################################################################################

	HMM.fisher.results<-matrix(NA, nrow=dim(HMM.results)[1], ncol=15)
        for (i in 1:dim(HMM.results)[1])
        {

	    if (distance.status[i]=="right" || distance.status[i]=="both" )
               {
		    
		    for.fisher<-matrix(NA,nrow=2, 3)
		    HMM.CG<-HMM.results[i,]
		    HMM.CG[is.na(HMM.CG)]<- 5 # assign NA any integer that is not confused with the analysis
        	    for.fisher[1,1:3]<-c(
		    as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==1))) + as.numeric(sum(1*(HMM.CG[(2+n1+n2): (2+n1+n2+n1-1)] ==1))), as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==0.5))) + as.numeric(sum(1*(HMM.CG[(2+n1+n2): (2+n1+n2+n1-1)] ==0.5))), as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==0))) + as.numeric(sum(1*(HMM.CG[(2+n1+n2): (2+n1+n2+n1-1)] ==0)))  )
        	    for.fisher[2,1:3]<-c(
		    as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==1)))+ as.numeric(sum(1*(HMM.CG[(2+2*n1+n2):(2+2*n1+2*n2-1)]==1))), as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==0.5)))+as.numeric(sum(1*(HMM.CG[(2+2*n1+n2):(2+2*n1+2*n2-1)]==0.5))), as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==0)))+ as.numeric(sum(1*(HMM.CG[(2+2*n1+n2):(2+2*n1+2*n2-1)]==0))) )

        	    HMM.fisher.results[i,1]<-as.numeric(mC.matrix[i,1])  # i position
        	    HMM.fisher.results[i,2]<-as.numeric(mC.matrix[i+1,1])  # i+1 position
        	    HMM.fisher.results[i,3]<-round(fisher.test(for.fisher)$p.value,4) # p-value
    
        	    # i position: estimated H for control and test
        	    HMM.fisher.results[i,4]<-paste(as.numeric(HMM.results[i, 2:(2+n1-1)]), collapse=":")    
        	    HMM.fisher.results[i,5]<-paste(as.numeric(HMM.results[i, (2+n1):(2+n1+n2-1)]), collapse=":")    

	            # i+1 position: estimated H for control and test
	            HMM.fisher.results[i,6]<-paste(as.numeric(HMM.results[i, (2+n1+n2):(2+2*n1+n2-1)]), collapse=":")    
	            HMM.fisher.results[i,7]<-paste(as.numeric(HMM.results[i, (2+2*n1+n2):(2+2*n1+2*n2-1)]), collapse=":")
	            # i position: raw methy level for control and test
	            HMM.fisher.results[i,8]<-paste(round(mC.matrix[i, 2:(2+n1-1)],2),collapse=":")
	            HMM.fisher.results[i,9]<-paste(round(mC.matrix[i,(2+n1):(2+n1+n2-1)],2),collapse=":")

	            # i+1 position: raw methy level for control and test
	            HMM.fisher.results[i,10]<-paste(round(mC.matrix[i+1, 2:(2+n1-1)],2),collapse=":")
	            HMM.fisher.results[i,11]<-paste(round(mC.matrix[i+1,(2+n1):(2+n1+n2-1)],2),collapse=":")
	            #i position: postrior prob for control and test
	            HMM.fisher.results[i,12]<-paste(round(as.numeric(HMM.results[i,(2+2*n1+2*n2):(2+3*n1+2*n2-1)]),2),collapse=":")
	            HMM.fisher.results[i,13]<-paste(round(as.numeric(HMM.results[i,(2+3*n1+2*n2):(2+3*n1+3*n2-1)]),2),collapse=":")
	            #i+1 position: postrior prob for control and test
	            HMM.fisher.results[i,14]<-paste(round(as.numeric(HMM.results[i,(2+3*n1+3*n2):(2+4*n1+3*n2-1)]),2),collapse=":")
	            HMM.fisher.results[i,15]<-paste(round(as.numeric(HMM.results[i,(2+4*n1+3*n2):(2+4*n1+4*n2-1)]),2),collapse=":")

                }

             else 
                {
	            for.fisher<-matrix(NA,nrow=2, 3)
        	    for.fisher[1,1:3]<-c(as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==1))), as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==0.5))), as.numeric(sum(1*(HMM.CG[2:(2+n1-1)]==0)))  )
        	    for.fisher[2,1:3]<-c(as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==1))), as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==0.5))), as.numeric(sum(1*(HMM.CG[(2+n1):(2+n1+n2-1)]==0))))
	            HMM.fisher.results[i,1]<-as.numeric(HMM.results[i,1])  # position
	            HMM.fisher.results[i,3]<-round(fisher.test(for.fisher)$p.value,4) # p-value
	    
	            #estimated H for control and test
        	    HMM.fisher.results[i,4]<-paste(as.numeric(HMM.results[i, 2:(2+n1-1)]), collapse=":")    
        	    HMM.fisher.results[i,5]<-paste(as.numeric(HMM.results[i, (2+n1):(2+n1+n2-1)]), collapse=":")    

	            #raw methy level for control and test
	            HMM.fisher.results[i,8]<-paste(round(mC.matrix[i, 2:(2+n1-1)],2),collapse=":")
	            HMM.fisher.results[i,9]<-paste(round(mC.matrix[i,(2+n1):(2+n1+n2-1)],2),collapse=":")

	            #postrior prob for control and test
	            HMM.fisher.results[i,12]<-paste(round(as.numeric(HMM.results[i,(2+2*n1+2*n2):(2+3*n1+2*n2-1)]),2),collapse=":")
	            HMM.fisher.results[i,13]<-paste(round(as.numeric(HMM.results[i,(2+3*n1+2*n2):(2+3*n1+3*n2-1)]),2),collapse=":")	   
                }
        }

        HMM.fisher.results<-cbind("chr1",HMM.fisher.results)
        HMM.fisher.results<-as.data.frame(HMM.fisher.results, stringsAsFactors = FALSE)
	HMM.fisher.results[,2:4]<-apply(HMM.fisher.results[,2:4], 2, function(x) as.numeric(x))
	return(HMM.fisher.results)

}
