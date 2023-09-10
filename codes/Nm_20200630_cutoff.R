#cutoff
#20200630

setwd("/Workplaces/project/Nm/20190821/data/Group_1/PA/20200607_test/analysis/03aligned/023dropReverseWithoutRmMultReads/results/03Calculating/rRNA/R1/parameres/") 
rm(list = ls())

mydata1<-read.table("G1_rRNA_value.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)
mydata2<-read.table("G2_rRNA_value.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)
mydata3<-read.table("G3_rRNA_value.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)


mydata<-as.data.frame(cbind(mydata1[c(2,11,12)],mydata2[12],mydata3[12]),stringsAsFactors = FALSE)

names(mydata)<-c("pos","note","G1","G2","G3")

head(mydata)

mydata$T1<-rep(0,dim(mydata)[1])
mydata$T2<-rep(0,dim(mydata)[1])
mydata$T3<-rep(0,dim(mydata)[1])
mydata$sum<-rep(0,dim(mydata)[1])


mydata$T1[which(mydata$G1>3)]<-1
mydata$T2[which(mydata$G2>3)]<-1
mydata$T3[which(mydata$G3>3)]<-1

mydata$sum<-mydata$T1+mydata$T2+mydata$T3

mydataX<-mydata[which(mydata$sum>=2),]


head(mydataX)
sle<-mydataX$pos
mydata1X<-mydata1[which(mydata1$pos %in% sle),]
mydata2X<-mydata2[which(mydata2$pos %in% sle),]
mydata3X<-mydata3[which(mydata3$pos %in% sle),]

mydata1XFilter<-mydata1X[which((mydata1X$depth.OED>10)&(mydata1X$chi1P<0.001)&(mydata1X$enrich2>=0.9)&(mydata1X$chi2P<0.001)),]
mydata2XFilter<-mydata2X[which((mydata2X$depth.OED>10)&(mydata2X$chi1P<0.001)&(mydata2X$enrich2>=0.9)&(mydata2X$chi2P<0.001)),]
mydata3XFilter<-mydata3X[which((mydata3X$depth.OED>10)&(mydata3X$chi1P<0.001)&(mydata3X$enrich2>=0.9)&(mydata3X$chi2P<0.001)),]

y<-intersect(intersect(mydata1XFilter$pos,mydata2XFilter$pos),mydata3XFilter$pos)

length(y)

length(mydata1XFilter$pos[which(mydata1XFilter$note>0)])
length(mydata2XFilter$pos[which(mydata2XFilter$note>0)])
length(mydata3XFilter$pos[which(mydata3XFilter$note>0)])

mydata1XFilter$pos[which(mydata1XFilter$note>0)]


k<-mydata1$note[which(mydata$note>0)]
length(k)
setdiff(y,k)

As16<-c(225109,2729819,3943146,4036869,4167997,4209485)
As23<-c(228273,2726670,3424269,3946218,4040034,4212557)

intersect(As16,y)

intersect(As23,y)

myresults<-mydataX[which(mydataX$pos %in% y),]

write.table(myresults,"20200620_erich1_3.txt",sep = "\t",row.names = FALSE,quote = FALSE)

write.table(mydata1XFilter,"PA_mydata1XFilter.txt",quote = FALSE,sep = "\t",row.names = FALSE)

