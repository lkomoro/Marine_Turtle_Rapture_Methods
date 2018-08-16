setwd("~/Box/Research Projects/NGS Projects and resources/Pacific leatherback rapture project/Dcor_Rapturedata/NEW_RaptureDataAnalysis/All_genotyping_QC_NEWDATA")
library(ggplot2)
library(reshape2)

#-------------
#Read in files:
SNPinfo<-read.csv("../kSNPkey_relposinfo.csv")
QCdf<-read.csv("../No.seq_QC_NEWDATAwithmetadata_added.csv")
QCdf$platewell<-QCdf$X

#Clean up/merge 0.8 PP: ##
All.0.8.geno <- read.table(gzfile("./QCPA_0.8_kSNP.geno.gz"), header=F)
All.0.8.geno[All.0.8.geno=="NN"]=NA
All.0.8.geno$no.samples.NA<-rowSums(is.na(All.0.8.geno))
All.0.8.geno$no.samples.GT<-(length(All.0.8.geno)-3)-All.0.8.geno$no.samples.NA
All.0.8.sumGT<-All.0.8.geno[,c(1,2,1101)]#summary
summary(All.0.8.sumGT$no.samples.GT)

All.0.8.geno$locus<-paste0(All.0.8.geno$V1,"_",All.0.8.geno$V2)
row.names(All.0.8.geno)<-All.0.8.geno$locus
All.0.8.geno<-All.0.8.geno[,c(3:1099)]
All.0.8.genot<-as.data.frame(t(All.0.8.geno))
samplelist<-read.delim("AllQCpassedsamples_bamlist_NEWDATA.txt", header=F)
samplelist$ID<-gsub("*_sortfltr.bam","",samplelist$V1)
All.0.8.genot.1<-cbind(samplelist,All.0.8.genot)
All.0.8.genot.1<-All.0.8.genot.1[,2:2009]
All.0.8.genot.1$no.lociNA<-rowSums(is.na(All.0.8.genot.1))
All.0.8.genot.1$no.loci.GT<-(length(All.0.8.genot.1)-2)-All.0.8.genot.1$no.lociNA
All.0.8.genot.1sum<-All.0.8.genot.1[,c(1,2010)]
vars<-colsplit(All.0.8.genot.1sum$ID, "_", c("species","LABID","DID","plate","well"))
All.0.8.genot.1sum$platewell<-paste0(vars$plate,vars$well)
All.0.8.genot.1.QC<-merge(QCdf,All.0.8.genot.1sum, by="platewell")
All.0.8.genokSNP <- read.table(gzfile("./QCPA_0.8_kSNP.geno.gz"), header=F)
All.0.8.genokSNP$Locus<-paste0(All.0.8.genokSNP$V1,"_",All.0.8.genokSNP$V2)
#row.names(All.0.8.genokSNP)<-All.0.8.genokSNP$locus
All.0.8.genokSNP1<-All.0.8.genokSNP[,c(1100,3:1099)]
All.0.8.genokSNP1[All.0.8.genokSNP1=="NN"]=NA 
All.0.8.genokSNP1$no.samplesNA<-rowSums(is.na(All.0.8.genokSNP1))
All.0.8.genokSNP1$no.samples.GT<-(length(All.0.8.genokSNP1)-2)-All.0.8.genokSNP1$no.samplesNA
All.0.8.genokSNP1.1<-All.0.8.genokSNP1[,c(1,1100)]
All.0.8.genokSNP1.2<-merge(All.0.8.genokSNP1.1,SNPinfo)

p<-ggplot(All.0.8.genokSNP1.2, aes(x=Relative_SNP_position, y=no.samples.GT)) +theme_bw()+
  geom_point(alpha=0.5)  +  
  ylab("No. samples Genotyped (0.8 PP)") +xlab("SNP position in Rapture Locus") +geom_vline(xintercept = 84,colour="red", linetype = "longdash")+
  geom_vline(xintercept = 134,colour="blue", linetype = "longdash")
s13.5<-p + annotate("text", x = 45, y = 50, label = "Only a priori SNP positions (2007 total)")
pdf("All.0.8 No. Samples GT vs. kSNP relative position.pdf", width=6, height=6)
s13.5 #Fig 2B
dev.off()

#Clean up/merge 0.9 PP: ##
All.0.9.geno <- read.table(gzfile("./QCPA_0.9_kSNP.geno.gz"), header=F)
All.0.9.geno[All.0.9.geno=="NN"]=NA
All.0.9.geno$no.samples.NA<-rowSums(is.na(All.0.9.geno))
All.0.9.geno$no.samples.GT<-(length(All.0.9.geno)-3)-All.0.9.geno$no.samples.NA
All.0.9.sumGT<-All.0.9.geno[,c(1,2,1101)]#summary
summary(All.0.9.sumGT$no.samples.GT)
All.0.9.geno$locus<-paste0(All.0.9.geno$V1,"_",All.0.9.geno$V2)
row.names(All.0.9.geno)<-All.0.9.geno$locus
All.0.9.geno<-All.0.9.geno[,c(3:1099)]
All.0.9.genot<-as.data.frame(t(All.0.9.geno))
samplelist<-read.delim("AllQCpassedsamples_bamlist_NEWDATA.txt", header=F)
samplelist$ID<-gsub("*_sortfltr.bam","",samplelist$V1)
All.0.9.genot.1<-cbind(samplelist,All.0.9.genot)
All.0.9.genot.1<-All.0.9.genot.1[,2:2009]
All.0.9.genot.1$no.lociNA<-rowSums(is.na(All.0.9.genot.1))
All.0.9.genot.1$no.loci.GT<-(length(All.0.9.genot.1)-2)-All.0.9.genot.1$no.lociNA
All.0.9.genot.1sum<-All.0.9.genot.1[,c(1,2010)]
vars<-colsplit(All.0.9.genot.1sum$ID, "_", c("species","LABID","DID","plate","well"))
All.0.9.genot.1sum$platewell<-paste0(vars$plate,vars$well)
All.0.9.genot.1.QC<-merge(QCdf,All.0.9.genot.1sum, by="platewell")

#Clean up/merge 0.95 PP: ##
All.0.95.geno <- read.table(gzfile("./QCPA_0.95_kSNP.geno.gz"), header=F)
All.0.95.geno[All.0.95.geno=="NN"]=NA
All.0.95.geno$no.samples.NA<-rowSums(is.na(All.0.95.geno))
All.0.95.geno$no.samples.GT<-(length(All.0.95.geno)-3)-All.0.95.geno$no.samples.NA
All.0.95.sumGT<-All.0.95.geno[,c(1,2,1101)]#summary
summary(All.0.95.sumGT$no.samples.GT)
All.0.95.geno$locus<-paste0(All.0.95.geno$V1,"_",All.0.95.geno$V2)
row.names(All.0.95.geno)<-All.0.95.geno$locus
All.0.95.geno<-All.0.95.geno[,c(3:1099)]
All.0.95.genot<-as.data.frame(t(All.0.95.geno))
samplelist<-read.delim("AllQCpassedsamples_bamlist_NEWDATA.txt", header=F)
samplelist$ID<-gsub("*_sortfltr.bam","",samplelist$V1)
All.0.95.genot.1<-cbind(samplelist,All.0.95.genot)
All.0.95.genot.1<-All.0.95.genot.1[,2:2009]
All.0.95.genot.1$no.lociNA<-rowSums(is.na(All.0.95.genot.1))
All.0.95.genot.1$no.loci.GT<-(length(All.0.95.genot.1)-2)-All.0.95.genot.1$no.lociNA
All.0.95.genot.1sum<-All.0.95.genot.1[,c(1,2010)]
vars<-colsplit(All.0.95.genot.1sum$ID, "_", c("species","LABID","DID","plate","well"))
All.0.95.genot.1sum$platewell<-paste0(vars$plate,vars$well)
All.0.95.genot.1.QC<-merge(QCdf,All.0.95.genot.1sum, by="platewell")

#Combined GT by sequence fragments for Comparisons between different posterior probability thresholds####
All.0.9.genot.1.QC$PPT<-0.9
All.0.8.genot.1.QC$PPT<-0.8
All.0.95.genot.1.QC$PPT<-0.95
All.comb<-rbind(All.0.9.genot.1.QC,All.0.8.genot.1.QC)
All.comb<-rbind(All.comb, All.0.95.genot.1.QC)
All.comb$PPT<-factor(All.comb$PPT)
All.comb<-All.comb[,c(1,2,21:23)]
p<-ggplot(All.comb, aes(x=No_seq, y=no.loci.GT,fill=PPT)) +theme_bw()+
  geom_point(alpha=0.75,shape=21,size=3)  +   
  #geom_jitter()+
  scale_fill_hue(l=40)+
  ylab("No. loci Genotyped ") +xlab("No. Sequenced Fragments/Individual")+geom_vline(aes(xintercept=150000),linetype="dashed",color="red")+geom_vline(aes(xintercept=500000),linetype="dashed",color="red")
s17<-p + annotate("text", x = 1000000, y = 1650, label = "kSNP Genotyped")
s17
pdf("AllQCpassed samples PP combined No. loci Genotyped at kSNP by No. seqs_wline NEWDATA.pdf", width=6, height=6)
s17 #Fig 2A
dev.off()

#examine asmyptote more closely:
ggplot(All.comb, aes(x=No_seq, y=no.loci.GT,fill=PPT)) +theme_bw()+
  geom_point(alpha=0.75,shape=21,size=3)  +   
  scale_fill_hue(l=40)+xlim(0,200000)+
  ylab("No. loci Genotyped ") +xlab("No. Sequenced Fragments/Individual")