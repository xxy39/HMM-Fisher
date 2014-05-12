joint.prob<-function(H, Obs, mu, sigma2, trans, initial.pi, lower.limit, upper.limit)
{ 
   # This function is used to calculate the joint probability P(H,O) 
   # P(H,O)=P(H)*P(O|H)= P(H_{1})*P(O_{1} | H_{1})* [ prod( P(H_{g} |H_{g-1}) * P(O_{g} | H_{g}) ) ]
   

   G<-length(H) # This is the number of CGs

   log.joint.p<-log( initial.pi[H[1]] ) # This is the initial distribution of each hidden state. 
   log.joint.p<-log.joint.p + log(dtnorm(Obs[1], mean=mu[H[1]], sd=sqrt(sigma2[H[1]]), lower=lower.limit[H[1]], upper=upper.limit[H[1]] ))
   # The above two lines give P(H_{1})*P(O_{1} | H_{1})
 
   for ( g in 2:G) 
   { 
         log.trans.g<-log(trans[H[g-1], H[g]]) 
 
         log.emiss.g<-log(dtnorm(Obs[g], mean=mu[H[g]], sd=sqrt(sigma2[H[g]]), lower=lower.limit[H[g]], upper=upper.limit[H[g]] ))
         log.joint.p<-log.joint.p + log.trans.g + log.emiss.g
    }
     
   return(log.joint.p) 
}

# joint.prob.txt 
