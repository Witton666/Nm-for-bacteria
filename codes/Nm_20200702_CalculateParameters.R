#20200622
#Nm mRNA

#add coverage
#OEDvsInput
#20200622

setwd("Workplaces/project/Nm/20190821/data/Group_1/PA/20200607_test/analysis/03aligned/03dropReverse/results/03Calculating/G_test/")
rm(list=ls())

a<-"PA-OED-3" ##change
b<-"PA-Input-3" ##change

OED<-read.delim(paste("../",a,".1.basic.txt",sep = ""),header=TRUE,stringsAsFactors = FALSE,sep = "\t") #R1
Input<-read.delim(paste("../",b,".1.basic.txt",sep = ""),header=TRUE,stringsAsFactors = FALSE,sep = "\t") #R1

#-----------add coverage------------
OEDc<-read.delim(paste("../../02coverage/",a,".1.ufilter.sorted_d.txt",sep = ""),header=FALSE,stringsAsFactors = FALSE,sep = "\t")
Inputc<-read.delim(paste("../../02coverage/",b,".1.ufilter.sorted_d.txt",sep = ""),header=FALSE,stringsAsFactors = FALSE,sep = "\t")

mydata<-as.data.frame(cbind(OED[c(1,2)],OEDc$V3,OED[c(3,5,6)],Inputc$V3,Input[c(3,5,6)]),stringsAsFactors = FALSE)
names(mydata)<-c("ref","pos","coverage.OED","depth.OED","Dsum.OED","Dave.OED","coverage.Input","depth.Input","Dsum.Input","Dave.Input")


#PA
#23S
G2237 <-c(726340, 4789488, 5265016, 6040500)
C2484 <-c(726587, 4789241, 5264769, 6040253)
U2538 <-c(726641, 4789187, 5264715, 6040199)

#16S
C1395 <-c(723491, 4792337, 5267865, 6043349)

mydata$note<-rep(0,dim(mydata)[1])
mydata$note[which(mydata$pos %in% G2237)]<-"G2237"
mydata$note[which(mydata$pos %in% C2484)]<-"C2484"
mydata$note[which(mydata$pos %in% U2538)]<-"U2538"

mydata$note[which(mydata$pos %in% C1395)]<-"C1395"


write.table(mydata,"G3.txt",quote = FALSE, row.names = FALSE,sep = "\t")  ##change

#------mRNA---------------------

setwd("/home/xindeng/Public/Workplaces/project/Nm/20190821/data/Group_1/PA/20200607_test/analysis/03aligned/03dropReverse/results/03Calculating/G_test/")
rm(list=ls())

G1<-read.delim("G1.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t")
G2<-read.delim("G2.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t")
G3<-read.delim("G3.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t")

S1<-seq(721775,727262)
S2<-seq(4788574,4793732)
S3<-seq(5276152,5281666)
S4<-seq(6051627,6057275)

x<-c(S1,S2,S3,S4)


"%ni%" <- Negate("%in%")

G1_mRNA<-G1[which(G1$pos %ni% x),]
G2_mRNA<-G2[which(G1$pos %ni% x),]
G3_mRNA<-G3[which(G1$pos %ni% x),]


write.table(G1_mRNA,"mRNA/G1_mRNA.txt",quote = FALSE, row.names = FALSE,sep = "\t")  ##change
write.table(G2_mRNA,"mRNA/G2_mRNA.txt",quote = FALSE, row.names = FALSE,sep = "\t")  ##change
write.table(G3_mRNA,"mRNA/G3_mRNA.txt",quote = FALSE, row.names = FALSE,sep = "\t")  ##change

#-------------------calculate parameters----------------
setwd("Workplaces/project/Nm/20190821/data/Group_1/PA/20200607_test/analysis/03aligned/03dropReverse/results/03Calculating/G_test/mRNA/")
rm(list=ls())

mydata<-read.delim("G3_mRNA.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t") ###change

#RATIO
ratio1<-rep(0,dim(mydata)[1])

mydata$enrich1<-log2(ratio1+1)






#chi_test_1
#Chi(Dbase_Nm,DAve_Nm,Dbase_input,DAve_input)
p<-rep(1,dim(mydata)[1])  #initial P
for(i in 1:dim(mydata)[1]){
  if(mydata$depth.OED[i]>0){
    y<-chisq.test(matrix(c(mydata$depth.OED[i],mydata$Dave.OED[i],mydata$depth.Input[i],mydata$Dave.Input[i]),nrow=2), simulate.p.value = TRUE, B = 2000)
    p[i]<-y$p.value 
  }
}
mydata$chi1P<-p


#RATIO2
ratio2<-rep(0,dim(mydata)[1])

#log2(Dbase_Nm/DAve_Nm)>=2
ratio2<-mydata$depth.OED/mydata$Dave.OED


#Log2(Ratio+offset)  offset=1

mydata$enrich2<-log2(ratio2)


#chi_test_2
#Chi(Dbase_Nm,DSum_Nm,DAve_Nm,DSum_Nm)
t1<-proc.time()
q<-rep(1,dim(mydata)[1])  #initial P
for(i in 1:dim(mydata)[1]){
  if(mydata$depth.OED[i]>0){
    g<-chisq.test(matrix(c(mydata$depth.OED[i],mydata$Dsum.OED[i],mydata$Dave.OED[i],mydata$Dsum.OED[i]),nrow=2),simulate.p.value = TRUE, B = 2000)
    q[i]<-g$p.value
  }
}
proc.time()-t1
mydata$chi2P<-q

write.table(mydata,"G3_mRNA_value.txt",quote = FALSE,row.names = FALSE,sep = "\t")  ###change







