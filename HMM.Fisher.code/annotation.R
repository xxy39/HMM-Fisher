# annotaion.R
#
# R script to provide gene annotation analysis for identified DM CG sites
#
# Usage: 
# R CMD BATCH '--args  input1  input2  distance header1 header2 output'   HMM.Fisher.code/annotation.R
#
###########################################################################################################
# 1) Input1: the DM.CG.txt results from HMM-Fisher program
###########################################################################################################
#
###########################################################################################################
# 2) Input2: The gene reference, from UCSC gene browser
# column names:
# bin  name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
# 585	NR_028269	chr1	-	4224	7502	7502	7502	7	4224,4832,5658,6469,6719,7095,7468,	4692,4901,5810,6631,6918,7231,7502,	0	LOC100288778	unk	unk	-1,-1,-1,-1,-1,-1,-1,  
###########################################################################################################
#
# 3) distance: the distance of the promoter regions. Promoter region for a specific gene is defined 
#              as the distance bp extended from the start and end of the gene.
# 4) header1: Character, whether the inpout1 has a header line. T, true; F, fasle.
#
# 5) header2: Character, whether the inpout2 has a header line. T, true; F, fasle.
#
###########################################################################################################
# output: The annotated CG sites. It contains 7 fields for each CG in DM.CG.txt. 
# column names:
# chr       pos     DM meanDiff.mC   meanCov     genes    promoters      
# chr1    795361  hypo    -0.4652    70.25:69   FAM41C           NA
# chr1    835742  hypo    -0.373  41.25:47.25       NA           NA
# chr1    841778  hypo    -0.4799     30:27.5       NA           NA
# chr1    851081  hypo    -0.6253    8.75:5.5    SAMD11          NA
# chr1    851104  hypo    -0.7374        9:5.5   SAMD11          NA
###########################################################################################################
#
# example command
# R CMD BATCH '--args  DM.CG.txt  refGene.txt  1000  T T   annotation.txt'   HMM.Fisher.code/annotation.R
 
 # function to process each chr
 procChr <- function(In, Ref){
      out<-matrix(NA, nrow=0, ncol=7); 
       len_in<-dim(In)[1];
     
         if (sum(!is.na(Ref))==0){
          out_1to4<-cbind(as.character(In[,1]), as.numeric(In[,2]) ,as.character(In[,19]),as.numeric(In[,17]), paste(as.character(In[,21]), as.character(In[,22]), sep=":"));
            out_full<-cbind(out_1to4, rep("NA", len_in),rep("NA", len_in));
           out<-rbind(out, out_full);
         }
      else{ 
          len_ref<-dim(Ref)[1];
          for ( i in 1:len_in ) {
             # get the position for input 1     
               pos<-as.numeric(In[i, 2]);
              
                # get the index of all genes containing the position i, 
               index<-(1:len_ref)[pos >= Ref[, 5] & pos <= Ref[, 6]];
             # get the index of all promoter regions containing the position i,
               index2<-(1:len_ref)[(pos >= (Ref[, 5]-dis) & pos < Ref[, 5] & as.character(Ref[,4]) == "+")|(pos <= (Ref[, 6]+dis) & pos > Ref[, 6] & as.character(Ref[,4]) == "-")];
             
                geneName<-"NA";
               promGene<-"NA";
               
                if (length(index)>0){
                     geneName<-paste(sort(unique(Ref[index,13])), collapse= ";");
                   }
               
                 if(length(index2)>0) {
                     promGene<-paste(sort(unique(Ref[index2,13])), collapse= ";");
                   }
               
                 out_1to4<-cbind(as.character(In[i,1]), as.numeric(In[i,2]) ,as.character(In[i,19]),as.numeric(In[i,17]), paste(as.character(In[i,21]), as.character(In[i,22]), sep=":") );
               out_full<-cbind(out_1to4, geneName, promGene);
               out<-rbind(out, out_full);
     }
   }
   
  return (out);
 }

 args<-commandArgs(trailingOnly = TRUE);
 
 dis<-as.integer(args[3]);

 header1<-as.logical(args[4]);
 header2<-as.logical(args[5]);

 Input <- read.table(file=args[1], header=header1);
 Refer <- read.table(file=args[2], header=header2); 

 # convert 0-base to 1-base position for txstart;
 Refer[,5]=Refer[,5]+1;
 
 hyper.index<-(1:dim(Input)[1])[Input[,19]=="hyper"]
 hypo.index<-(1:dim(Input)[1])[Input[,19]=="hypo"]
 DM<-rep(NA,dim(Input)[1])
 DM[hyper.index]<-"hyper"
 DM[hypo.index]<-"hypo"

 Input[,19]<-DM
 
 # get the all unqiue chr name
 uniq_chr<-unique(Input[,1]);
 
 out<-matrix(NA, nrow=0, ncol=7); 
 
 # process each chrs
 for (i in 1:length(uniq_chr)){ 
   
   index<-(1:dim(Input)[1])[Input[,1] == uniq_chr[i]];
   index_Ref<-(1:dim(Refer)[1])[as.character(Refer[,3]) == as.character(uniq_chr[i])];
  
  In_chrI<-Input[index,];
   Ref_chrI<-NA;
  
  if (length(index_Ref)>0) {
    Ref_chrI<-Refer[index_Ref,];
  }
 
  out_chrI<-procChr(In_chrI, Ref_chrI);
 
  out<-rbind(out, out_chrI)
  
 }

 write.table(out, file =args[6], quote=F, row.names=F, col.names=F, sep="\t"); rm(list=ls());
