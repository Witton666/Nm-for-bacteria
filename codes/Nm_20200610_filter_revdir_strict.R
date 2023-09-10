rm(list=ls())


Getdropreads<-function(a){
  fwd<-read.delim(paste(a,".fwd_d.txt",sep=""),header = FALSE,stringsAsFactors = FALSE,sep = "\t") #
  rev<-read.delim(paste(a,".rev_d.txt",sep=""),header = FALSE,stringsAsFactors = FALSE,sep = "\t") #
  
  mydata<-as.data.frame(cbind(fwd,rev$V3,myanno$V2),stringsAsFactors = FALSE)
  names(mydata)<-c("ref","position","fwd","rev","refstrand")
  
  mydata1<-mydata$position[which((mydata$fwd>0)&(mydata$refstrand == "-"))] #delete
  mydata2<-mydata$position[which((mydata$rev>0)&(mydata$refstrand == "+"))] #delete
  mydata3<-mydata$position[which((mydata$fwd>0)&(mydata$refstrand == "."))] #add direction
  mydata4<-mydata$position[which((mydata$rev>0)&(mydata$refstrand == "."))] #add direction
  mydata5<-mydata$position[which((mydata$fwd>0)&(mydata$rev>0)&(mydata$refstrand == "."))] #delete
  
  
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # 
  # BiocManager::install("Rsamtools")
  #load library
  library(Rsamtools)
  
  #read in entire BAM file
  mybam<-scanBam(paste("../",a,".sorted.bam",sep = ""))
  
  mybam2<-as.data.frame(mybam)
  
  mybam2$note<-rep(0,dim(mybam2)[1])
  
  mydata5.del<-mybam2$pos[which(mybam2$pos %in% mydata5)]
  mydata1.del<-mybam2$pos[which(mybam2$pos %in% mydata1)]
  mydata2.del<-mybam2$pos[which(mybam2$pos %in% mydata2)]
  
  mybam2$note[which(mybam2$pos %in% mydata5.del)]<-1
  
  mybam2$note[which((mybam2$pos %in% mydata1.del)&(mybam2$strand == "+"))]<-1
  mybam2$note[which((mybam2$pos %in% mydata2.del)&(mybam2$strand == "-"))]<-1
  
  
  dropreads<-mybam2[which(mybam2$note == 1),1]
  
  write.table(dropreads,paste(a,".drop.reads.txt",sep=""),row.names = FALSE,quote = FALSE,sep = "\t",col.names = FALSE)
  
}


#PA
setwd("Nm/Nm/20190821/data/Group_1/PA/20200607_test/analysis/03aligned/023dropReverseWithoutRmMultReads/")
myanno<-read.delim("Workplaces/ref_all/PAO1/PAO1_base_strand.txt",header = FALSE,stringsAsFactors = FALSE,sep = " ")
mylist<-read.table("../../list.txt",header = F,stringsAsFactors = FALSE)$V1

for(k in mylist){
  Getdropreads(k)
}

