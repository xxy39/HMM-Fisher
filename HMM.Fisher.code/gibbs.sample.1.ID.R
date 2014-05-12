gibbs.sample.1.ID<-function(H, Obs, mu, sigma2, trans, initial.pi, lower.limit, upper.limit)
{ 
   # This function is used to update the states in the hidden sequence H for each individual: 1, 2, 3, 
   # that is, 0, 1/2, 1, (not methylated, partly methylated, and fullly methylated respectively). 
   
  
   GG<-length(Obs) # This is the number of CG sites. 

   H.condi.prob<-matrix(NA, nrow=3, ncol=GG) 
   # I set up the above matrix to save the conditional probabilities 
   
   for ( g in 1:GG) 
   { 
      # if (g%%20000==0) { cat("-----------------when g is:", g, "\n") } 
      condi.prob.g<-rep(NA,3) 
      for ( h in 1:3) 
      { # h  in 1:3 corresponds to the hidden states {0, 1/2, 1}       
        if ( g==1) 
        { log.trans.g<-log(trans[h, H[g+1]]) 
          # log.emiss.g<-log( dnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h]) ))
          log.emiss.g<-log( dtnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h]), lower=lower.limit[h], upper=upper.limit[h]) )
          condi.prob.g[h]<-log(initial.pi[h])+log.trans.g + log.emiss.g
        }
        if ( g > 1 && g < GG) 
        { log.trans.g<-log(trans[H[g-1], h]) + log(trans[h, H[g+1]])
          # log.emiss.g<-log( dnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h]) ) )
          log.emiss.g<-log( dtnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h]), lower=lower.limit[h], upper=upper.limit[h]) )
          condi.prob.g[h]<-log.trans.g + log.emiss.g 
        }  

        if ( g == GG) 
        { log.trans.g<-log(trans[H[g-1], h]) 
          # log.emiss.g<-log( dnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h])))
          log.emiss.g<-log( dtnorm(Obs[g], mean=mu[h], sd=sqrt(sigma2[h]), lower=lower.limit[h], upper=upper.limit[h]))
          condi.prob.g[h]<-log.trans.g + log.emiss.g
        } 
     }
     
     H.condi.prob[, g]<-exp(condi.prob.g)/sum(exp(condi.prob.g))
     H[g]<-sample(c(1, 2, 3 ), 1, prob=exp(condi.prob.g))
     # cat("when g is", g, "the condional prob is", condi.prob.g, "the sampled H is", H[g], "\n") 
   }

  list(H=H, H.condi.prob=H.condi.prob) 
} 
    
# gibbs.sample.txt 
