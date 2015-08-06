get.distance.status<-function(positions, distance.cutoff)
 {
    # This function is to code each CG site based on the distance between itself and two flanking CG sites
    # coding system:
    # single: the distance between this CG sites and two flanking CG sites are both larger than distance.cutoff
    # left: the left flanking CG site is within distance.cutoff bp, while the right flanking CG site is  >distance.cutoff bp away.
    # right: the left flanking CG site is > distance.cutoff bp away, while the right flanking CG site is within distance.cutoff bp.
    # both: the left and right flanking CG site are both within distance.cutoff bp from this CG sites.

    # parameters:
    # 1) position: vector of physical positions of all CG sites along the chr
    # 2) distance.cutoff: a single number represents the cutoff of distance in coding system

    # output: a vector represents the codes for each CG sites, in the physical order

    GG<-length(positions) # This is the number of CG sites. 

    distance.status<-rep(NA, GG)
    # a vector to store the code for each CG site

    ###########################################
    # the first CG site
    ###########################################
    if (positions[2]-positions[1]<=distance.cutoff) # the fisrt CG site and second CG sites are not too far away
      { distance.status[1]<-"right"} # assign as right 
    else   #  the fisrt CG site and second CG sites are far away, only consider the first CG site itself
      {distance.status[1]<-"single" } # assign as singl
 


    ###########################################
    # 2nd CG site - last second CG site
    ###########################################

    for ( g in 2:(GG-1)) 
      { 
         # if (g%%20000==0) { cat("-----------------when g is:", g, "\n") } 
      
         if (positions[g]-positions[g-1]>=distance.cutoff &&  positions[g+1]-positions[g]<=distance.cutoff)
           { distance.status[g]<-"right"} # assign as right  
         else if ( positions[g]-positions[g-1]<=distance.cutoff &&  positions[g+1]-positions[g]<=distance.cutoff)
           { distance.status[g]<-"both" } # assign as both                   
         else if ( positions[g]-positions[g-1]<=distance.cutoff &&  positions[g+1]-positions[g]>=distance.cutoff)
           { distance.status[g]<-"left" } # assign as left            
         else
          { distance.status[g]<-"single"} # assign as single        
      }


    ###########################################
    # the last CG sites
    ###########################################
    if (positions[GG]-positions[GG-1]<=distance.cutoff) # the last CG site and secondlast CG sites are not too far away
      {  distance.status[GG]<-"left" } # assign as left     
    else   #  the last CG site and secondlast CG sites are far away, only consider the first CG site itself
      { distance.status[GG]<-"single" } # assign as single	 
    
    return (distance.status)
}
