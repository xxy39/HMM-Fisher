
HMM.Fisher<-function(total.reads, meth.reads,  n1,n2, chromosome, code.dir, output.dir, min.percent=0.8,  iterations=60, dist.combine=100,p.threshold=0.05, meanDiff.cut=0.3,  max.distance=100, max.empty.CG=3, max.EM=1, max.p=0.1, singleton=TRUE )
  {
  
        #####################################################################################################################
	# run HMM-DM with a single command
	#
	# 1. meth.reads: a matrix of methylation levels in two groups. col1, positions; col2-5, methylation levels for 4 samples in group 1, 
	#                col6-9, methylation  levels for 4 samples in group 2.
	# 2. total.reads: a matrix of covreage for all samples in two groups. col1, positions; col2-5, coverage for 4 samples in group 1, col6-9, 
	#                coverage for 4 samples in group 2. Note that the positions and order of samples should correspond to those in mC.matrix.
        # 3. n1:  Numeric. number of control samples
	# 4. n2:  Numeric. number of control samples
	# 5. chromosome:Character. The chromosome the users want to analyze, e.g., chromosome =1, or chromosome = 2.
	# 6. code.dir: String. The directory of the source code files of HMM-DM (e.g., /home/HMM.DM/HMM.DM.code).
	#                Note, there should be no "/" at the very end.
	# 7. output.dir: String. The directory for the output files (e.g., /home/HMM.DM.results). Note, there should be no "/" at the very end.
	#                 Two matrices will be generated from this function. See section 3.2 for more detail. 
	# 8. min.percent: Numeric. The CG sites should be covered in at least min.percent of the control samples AND of the test samples. 
	#                 Otherwise, the CG sites are dropped. Default = 0.8.
	# 9. iterations: Numeric. The number of iterations when running HMM-DM. Default = 60
	# 10. dist.combine: Numeric. COncecutive CG sites with distance > dist.combine will not be combined in Fisher exact step. Default = 100bp
	# 11. p.threshold: Numeric.the threshold of p-value for significance. CG sites with p-value <= p.threshold will be identified as DM.
	# 12. meanDiff.cut: Numeric. Minimum mean difference of methylation levels between the two groups to call a DM CG site. Default = 0.3
	# 13. chromosome:Character. The chromosome the users want to analyze, e.g., chromosome =1, or chromosome = 2.
	# 14. max.distance: Numeric. The maximum distance between any two DM CGs site within a DM region. Default = 100bp.
	# 15. max.empty.CG: Numeric. The maximum number of CGs without coverage between any two DM CG sites within a DM region. Default = 3.
	# 16. max.EM: Numeric. When combining two consecutive DM regions (with >=2 DM CGs), the maximum number of EM CG sites between these
	#            two DM regions. These EM CG sites can be 1) having a p-value >= p.threshold, but <= max.p; 
	#            or 2) identified as DM by HMM-DM but with small meanDiff. Default = 1. 
        #            Note: if either region is a singleton, only 1 EM CG is allowed.
	# 17. max.p: Numeric. The maximum p-value for the EM included in the combined DM region. Default = 0.1.
	# 18.singleton: Logical. Report the singletons or not in summarizing region step? If TRUE (default), the singletons will be 
	#                reported in the HMM.DM.results.txt. See section 3.2 for more detail.
        #####################################################################################################################

        library("msm")

	# The following scource code files are used in quality control
        source(paste(code.dir,"getMeth.R",sep="/"))

	# The following scource code files are used in HMM-fisher 
        source(paste(code.dir,"gibbs.sample.1.ID.R",sep="/"))
        source(paste(code.dir,"check.post.H.R",sep="/"))
        source(paste(code.dir,"joint.prob.R",sep="/"))
        source(paste(code.dir,"get.H.max.string.ver2.R",sep="/"))
        source(paste(code.dir,"count.T.R",sep="/"))
        source(paste(code.dir,"up.T.prob.R",sep="/"))
        source(paste(code.dir,"get.distance.status.R",sep="/"))
        source(paste(code.dir,"fisher.exact.R",sep="/"))

	# The following scource code files are used to form regions 
        source(paste(code.dir,"chr.DM.region.by.CG.HMM.fisher.R",sep="/"))
        source(paste(code.dir,"DM.region.combine.HMM.fisher.R",sep="/"))
        source(paste(code.dir, "get.DM.region.R",sep="/"))    
        source(paste(code.dir, "DMR.combine.R",sep="/"))    
	########################################################################################################################
	# 1. quality control
	########################################################################################################################
        # get the methy ratio for CG sites that pass the quality control
	QC.results<-getMeth(total.reads, meth.reads, n1, n2, min.percent )
        mC.matrix<-round(QC.results[[1]],6)
        mC.matrix.index<-QC.results[[2]]                
        write.table(mC.matrix, paste(output.dir, "mC.matrix.txt", sep="/"), quote=F, col.names=F, row.names=F, sep="\t")

	########################################################################################################################
	# 2. Identify DM CG sites
	########################################################################################################################
        #################################################################################
        # 2.1 HMM step
        #################################################################################
	set.seed(13874312)  
        # define prior
        mu<-c(0, 1/2, 1)  # the mean of the emission probability (truncated Normal distribution) 
        sigma2<-c(0.12^2, 0.15^2, 0.13^2) # variance of emission probability
        initial.pi<-c(1/3, 1/3, 1/3)  # intial transition probability 
        lower.TN<-c(0, 0, 0)
        upper.TN<-c(1, 1, 1) 
        T.prior<-matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow=TRUE)  # transition prior Dirichlet (1, 1, 1) 

        # define the number of iterations
        R<-iterations       
        G<-dim(mC.matrix)[1] # Number of CG sites we are concerned 

        # vecotr to sane transtion probabilities
        trans.p<-rep(0,9)

        ##########################################################
        # 2.1.1 seperate HMM chain for each sample
        # set the same seed for each sample
        # estimate a single transition matrix using all 8 samples 
        ##########################################################         
         sample.H.list<-list()
         H.list<-list()
         H.index.list<-list()
	 for (j in 1:(n1+n2))
	   {              
		  # the index of CG with coverage in each sample
		  H.index<- (1:dim(mC.matrix)[1])[!is.na(mC.matrix[,1+j])*1]
		  H.index.list<-c(H.index.list, list(H.index))

		  # initial states for each sample
		  H<-sample(c(1, 2, 3), length(H.index), replace=T)
                  H.list<-c(H.list, list(H))

                  # emtpy matrix to store udated H for each iteration            
                  sample.H.list<-c(sample.H.list, list(matrix(rep(NA), nrow=R, ncol=length(H.index))))                  
	   }
         # empty matrix to store joint probability
         joint.p.matrix<-matrix(rep(NA),nrow=(n1+n2), ncol=R)

        ########################
        # 2.1.2 use the estimated transition matrix to update states for each sample
        #########################
	date(); 
        for ( i in 1: R)
        {
           cat("--------------- when iteration is", i," , time is", date(), "\n") 
           # H<-cbind(H.control1, H.control2,H.control3,H.control4,  H.test1, H.test2,  H.test3, H.test4)  # combine H from all 8 samples
           trans.list<-lapply(1:(n1+n2), function(x) { T.count<-count.T(H.list[[x]]); return(up.T.prob(T.count,T.prior))  } )  
	   # count the 9 transitions for each sample
           trans<-Reduce('+', trans.list)  # add the 9 transitions of 8 samples together
        
           for (j  in 1:(n1+n2) )
	     {                     
		     xx<-gibbs.sample.1.ID(H.list[[j]], Obs=mC.matrix[H.index.list[[j]],1+j], mu, sigma2, trans, initial.pi, lower.TN, upper.TN) # update the H using gibbs.sample
                     sample.H.list[[j]][i, ]<-xx$H
                     H.list[[j]]<-xx$H
                     joint.p.matrix[j,i]<-joint.prob(H.list[[j]], Obs=mC.matrix[H.index.list[[j]],1+j], mu, sigma2, trans, initial.pi, lower.limit=lower.TN, upper.limit=upper.TN)
	     }
        }
        date(); 

        ##################################
        # 2.1.3 summerize the estimated H for each sample
        ##################################
	xx.list<-list()
        for (j  in 1:(n1+n2) )
          {
                     
                    H.summary<-check.post.H(sample.H.list[[j]][(0.5*R+1):R, ]) # summerize H based on last half iterations 
                    xx<-get.H.max.string.ver2(H1=t(H.summary))
                    xx.list<-c(xx.list, list(xx))
	   }

        # plot the joint probability over iterations
	postscript(paste(output.dir, "joint.prob.ps", sep="/"), paper="letter", horizontal=T)
        par(mfrow=c(2,2))
	  for (j in 1:(n1))
	   { plot (1:R, joint.p.matrix[j,], xlab="iterations", ylab="joint prob" , main=paste("Joint Probability control ",j,sep=""))  }
	  for (j   in 1:(n2))
	   { plot (1:R, joint.p.matrix[n1+j,], xlab="iterations", ylab="joint prob" , main=paste("Joint Probability test ",j,sep="")) }	   
        dev.off()

        ##################################
        # 2.1.4 output the HMM results
        ##################################
        H.1<-matrix(NA, nrow=G-1, ncol=0)
        H.2<-matrix(NA, nrow=G-1, ncol=0)
	post.1<-matrix(NA, nrow=G-1, ncol=0)
	post.2<-matrix(NA, nrow=G-1, ncol=0)
	for (j in 1:(n1+n2))
	       {
                    xx.H<-rep(NA, G)
		    xx.H[H.index.list[[j]]]<-xx.list[[j]][,5]                    
		    xx.post<-rep(NA, G)
		    xx.post[H.index.list[[j]]]<-xx.list[[j]][,4]                 
		    H.1<-cbind(H.1,xx.H[1:G-1])
		    H.2<-cbind(H.2,xx.H[2:G])
		    post.1<-cbind(post.1,xx.post[1:G-1])
		    post.2<-cbind(post.2,xx.post[2:G])
	       }
        HMM.results<-cbind(mC.matrix[1:(G-1),1], H.1, H.2, post.1, post.2)

        #################################################################################
        # 2.2 Fisher exact test step
        #################################################################################     
	# check the distance status of CGs
        distance.status<-get.distance.status(positions=mC.matrix[,1], distance.cutoff=dist.combine)
        # perform Fisher's exact test
        HMM.fisher.results<-fisher.exact(HMM.results, mC.matrix, distance.status, n1, n2)
        #################################################################################
        # 2.3 Define DM bases on p-value and mean differences
        #################################################################################
	# calculate the meanDiff for each CG
        mean1<-apply(mC.matrix[,2:(2+n1-1)],1, function(x) { return(mean(x[!is.na(x)]))})
        mean2<-apply(mC.matrix[,(2+n1):(2+n1+n2-1)],1, function(x) { return(mean(x[!is.na(x)]))})
        meanDiff.2<- round((mean1-mean2),4)
        meanDiff<-meanDiff.2[1:dim(HMM.results)[1]]

	# DM status defined based on p-value
        Hyper.index<-(1:dim(HMM.fisher.results)[1])[HMM.fisher.results[,4] <= p.threshold & meanDiff >= 0]
        Hypo.index<-(1:dim(HMM.fisher.results)[1])[HMM.fisher.results[,4] <= p.threshold & meanDiff < 0]
        DM.status<-rep("EM", dim(HMM.fisher.results)[1]) # 0 for EM status
	DM.status[Hyper.index] <- "hyper" # 1 for hyper status (control > test)
	DM.status[Hypo.index] <- "hypo" # 1 for hypo status (control < test)

        # filter results by mean difference
	# DM CGs (hyper anf  hypo) with mean difference < meanDiff.cut are re-classified as EM
	Hyper.2.index<-(1:dim(HMM.fisher.results)[1])[DM.status == "hyper" & abs(meanDiff) < meanDiff.cut]
	Hypo.2.index<-(1:dim(HMM.fisher.results)[1])[DM.status == "hypo" & abs(meanDiff) < meanDiff.cut]
	DM.status.meandiff<-DM.status
	DM.status.meandiff[c(Hyper.2.index,Hypo.2.index )] <- "EM" # re-classified as 0 EM

        # calculate the mean coverage for both groups
        meanCov.test<- round(apply(total.reads[mC.matrix.index,2:(2+n1-1)],1, function(x) { return(mean(x[x>0]))}),2)
        meanCov.control<- round(apply(total.reads[mC.matrix.index,(2+n1):(2+n1+n2-1)],1,function(x) { return(mean(x[x>0]))}),2)

        HMM.fisher.results<-cbind(HMM.fisher.results, meanDiff, DM.status,DM.status.meandiff, c(1:1:dim(HMM.fisher.results)[1]) ,meanCov.test[1:1:dim(HMM.fisher.results)[1]],meanCov.control[1:1:dim(HMM.fisher.results)[1]])
        HMM.fisher.results<-as.data.frame(HMM.fisher.results, stringsAsFactors = FALSE)
	HMM.fisher.results[,c(17,20:22)]<-apply(HMM.fisher.results[,c(17,20:22)], 2, function(x) as.numeric(x))

        # column names
	colnames(HMM.fisher.results)<-c("chr", "pos", "pos2", "p.value", "test.s", "con.s", "test.s2",  "con.s2" , "test.mC",  "con.mC", "test.mC2", "con.mC2", "test.post", "con.post", "test.post2", "con.post2", "meanDiff", "DM.p", "DM.status", "index", "meanCov.test", "meanCov.control")
        
	DM.CGs<-HMM.fisher.results[HMM.fisher.results[,19]!= "EM",]
        write.table(DM.CGs, paste(output.dir, "DM.CG.txt", sep="/"), quote=F, col.names=c("chr", "pos", "pos2", "p.value", "test.s", "con.s", "test.s2",  "con.s2" , "test.mC",  "con.mC", "test.mC2", "con.mC2", "test.post", "con.post", "test.post2", "con.post2", "meanDiff", "DM.p", "DM.status", "index", "meanCov.test", "meanCov.control"), row.names=F, sep="\t")
        write.table(HMM.fisher.results, paste(output.dir, "all.CG.txt", sep="/"), quote=F, col.names=c("chr", "pos", "pos2", "p.value", "test.s", "con.s", "test.s2",  "con.s2" , "test.mC",  "con.mC", "test.mC2", "con.mC2", "test.post", "con.post", "test.post2", "con.post2", "meanDiff", "DM.p", "DM.status", "index", "meanCov.test", "meanCov.control"), row.names=F, sep="\t")

        #################################################################################
        # 3. summarize DM CGs into regions
        #################################################################################
        HMM.fisher.region.1<-chr.DM.region.by.CG.HMM.fisher(HMM.fisher.results, raw.CG=total.reads[,1], chr=chromosome, distance.threshold=max.distance, report.singleCG=singleton, empty.CG=max.empty.CG )

        DMRs<-DM.region.combine.HMM.fisher(HMM.fisher.region.1, HMM.fisher.results,  raw.CG=total.reads[,1],  distance.threshold=max.distance, num.CG.between=max.EM, p.threshold=max.p ,report.singleCG=singleton, empty.CG=max.empty.CG )

        colnames(DMRs)<-c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "meanCov.test", "meanCov.control", "meanDiff")
	write.table(DMRs, paste(output.dir, "DMRs.txt", sep="/"), quote=F, col.names=c("chr", "start", "end", "len", "DM",  "num.CG", "total.CG",  "meanCov.test", "meanCov.control", "meanDiff"), row.names=F, sep="\t")

}

