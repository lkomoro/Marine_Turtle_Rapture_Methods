#Setup####
setwd("~/...")
list.files()

library(ggplot2)
library(tidyr)
library(data.table)
#Multiplot Functions####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##### I. Initial QC checks of raw read data ########
###Data read in, basic clean up:####
datarenamed<-read.table("Linecount_bysample_renamedNEWDATA_LK.txt", header=T, sep="")
datarenamed$log_No_seqs<-log(datarenamed$No_seq)

samplesonly<-subset(datarenamed, species !="emptywell")
levels(samplesonly$species<-factor(samplesonly$species))#double check for typos in categories

emptywells<-subset(datarenamed, species =="emptywell")
levels(emptywells$species<-factor(emptywells$species))

#Summaries by plate for emptywells:
mean<-tapply(emptywells$No_seq, emptywells$species,mean)
min<-tapply(emptywells$No_seq, emptywells$species,min)
max<-tapply(emptywells$No_seq, emptywells$species,max)
median<-tapply(emptywells$No_seq, emptywells$species,median)
summary_emptywell<-data.frame(mean,median,min,max)

#Summaries by plate:
mean<-tapply(samplesonly$No_seq, samplesonly$plate,mean)
min<-tapply(samplesonly$No_seq, samplesonly$plate,min)
max<-tapply(samplesonly$No_seq, samplesonly$plate,max)
median<-tapply(samplesonly$No_seq, samplesonly$plate,median)
summary_by_plate<-data.frame(mean,median,min,max)

#####Histograms of raw reads by plate, with binning for thresholds####
sample.reads.bin<-samplesonly #from above
sample.reads.binRA<-subset(sample.reads.bin,sequence_direction=="RA.fastq")
sample.reads.binRA$Bin_seq <- ifelse(sample.reads.binRA$No_seq <= 10000, 10000, ifelse((sample.reads.binRA$No_seq >10000) & (sample.reads.binRA$No_seq < 100000), sample.reads.binRA$No_seq, 100000))
pdf("Turtle Rapture-No. Reads Binned.pdf", width=10, height=6)
ggplot(sample.reads.binRA, aes(x=Bin_seq)) + 
  geom_histogram(colour="black", fill="blue",binwidth = 10000) +theme_bw()+
  facet_wrap(~plate,scales="free")
dev.off()

sample.reads.binRA$QC_cat <- ifelse(sample.reads.binRA$No_seq <= 10000, "failed", ifelse((sample.reads.binRA$No_seq >10000) & (sample.reads.binRA$No_seq <= 100000), "good", "high"))
sample.reads.binRAshort<-sample.reads.binRA[,c(2:7,10,11)]
sample.reads.binRAshort$Bin_seqlow<-ifelse(sample.reads.binRAshort$No_seq <= 20000, sample.reads.binRAshort$No_seq, 20000)#bin everything high to be able to look at low distribution and confirm we think 10K is reasonable cutoff for pass/fail
pdf("Turtle Rapture-No. Reads Binned LOW.pdf", width=10, height=6)
ggplot(sample.reads.binRAshort, aes(x=Bin_seqlow)) + 
  geom_histogram(colour="black", fill="blue",binwidth = 1000) +theme_bw()+
  facet_wrap(~plate,scales="free")
dev.off()

failed<-subset(sample.reads.binRAshort, QC_cat=="failed")
table(failed$species)

#Patterns with metadata for failed samples?####
QCdf<-read.csv("No.seq_QC_NEWDATAwithmetadata_added.csv")
plot(QCdf$Year_Collected,QCdf$No_seq,las=2)#many samples are missing year info so undet., but samples known to be collected prior to 1994 didn't work; then variable (likely locations and other covariates in there)
plot(QCdf$Country,QCdf$No_seq,las=2)
plot(QCdf$Loc_combined,QCdf$No_seq,las=2)
plot(QCdf$Sex,QCdf$No_seq,las=2)
plot(QCdf$Capture_type,QCdf$No_seq,las=2)
plot(QCdf$Collection_Method,QCdf$No_seq,las=2)
plot(QCdf$Stage_Class,QCdf$No_seq,las=2)

#####II. Mapping Stats########
#Comparing filtered (rmdup) and unfiltered####
#Read in data and combine####
library(reshape2)
list.files(pattern="All_head")
mapped.raw<-read.delim("All_head_flagstat_reformat.txt",header=T)
mapped.fltr<-read.delim("All_head_sortfltr_flagstat_reformat.txt",header=T)
mapped.raw$ID<-gsub("_sort_flagstat.txt", "", mapped.raw$sampleID)
mapped.fltr$ID<-gsub("_sortfltr_flagstat.txt", "", mapped.fltr$sampleID)
mapped.comb<-merge(mapped.raw,mapped.fltr, by="ID")
mapped.comb<-mapped.comb[,c(1,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24,13,25)]
vars<-colsplit(mapped.comb$ID, "_", c("species","LABID","DID","plate","well"))
mapped.comb<-cbind(mapped.comb,vars)
mapped.comb<-mapped.comb[,c(1,24:28,2:23)]

mapped.comb<-subset(mapped.comb, species!="emptywell")
mapped.comb$species<-factor(mapped.comb$species)#reset
mapped.comb$prop.NDreads<-mapped.comb$mapped.reads.QC.passed.y/mapped.comb$total.reads.QC.passed.x
summary(mapped.comb$prop.NDreads)#22-27% of raw reads are left after mapping/filtering
sd(mapped.comb$prop.NDreads)
hist(mapped.comb$prop.NDreads)
mapped.comb$prop.dupreads<-(1-mapped.comb$prop.NDreads)
summary(mapped.comb$prop.dupreads)
sd(mapped.comb$prop.dupreads)
hist(mapped.comb$prop.dupreads)

mapped.comb.short<-mapped.comb[,c(1:12,29)]
str(mapped.comb.short)
names(mapped.comb.short) <- c("sampleID", "species", "LABID",	"D_ID", "plate"	,"well", "UF_total_reads_QC_passed"	,"F_total_reads_QC_passed","UF_mapped_reads_QC_passed","F_mapped_reads_QC_passed",	"UF_percent_reads_mapped","F_percent_reads_mapped","prop")

#add in QC category info from section I:
SRA1<-sample.reads.binRAshort[,c(1,5:9)]
SRA1$platewell<-paste(SRA1$plate,SRA1$well)
mapped.comb.short$platewell<-paste(mapped.comb.short$plate,mapped.comb.short$well)
CFS<-merge(mapped.comb.short,SRA1,by="platewell")
CFS.ok<-subset(CFS,QC_cat!="failed")

#Input concentration matters:
DNAconc<-read.csv("Rapture_inputDNAconckey.csv")
CFS.ok1<-merge(DNAconc,CFS.ok)
CFS.ok1$concentration<-factor(CFS.ok1$concentration, levels=c("<=5","5","10","20","40","variable"))
DNAbox<-ggplot(CFS.ok1, aes(x=concentration, y=prop)) +theme_bw()+
  geom_boxplot(colour="black", fill="red", alpha=0.7)  +  
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  ylab("Proportion Filtered Mapped/Sequenced Fragments") +xlab("DNA Input Original Concentration (ng/ul)") 

pdf("Turtle Rapture-prop filtered mapped-seq frag by DNA input.pdf", width=6, height=6)
DNAbox #Figure 1b
dev.off()

#### Section III. Bedtools coverage summaries and graphics####
# Using relative position 20 and bedtools coverage estimates (see upstream markdown and bash scripts):####
data<-read.table("./All.combined.pos20.2.bed", header=FALSE)
names(data)<-c("scaffold", "position","coverage","sampleID")
data1 <- data[grep("emptywell", data$sampleID, invert=T), ]

#Calculate coverage by locus, all QC Passed Samples  ####
#(based on 10K read threshold)
QC<-read.csv("No.seq_QC_NEWDATA.csv")
QC$platewell<-paste(QC$plate,QC$well,sep="")
qc.1<-QC[,c(8,10)]
z<-colsplit(data1$sampleID,"_",c("species","LABID","DID","plate","well"))
data1$platewell<-paste(z$plate,z$well,sep="")
data1.qc<-merge(data1,qc.1,by="platewell")  
data1.ok<-subset(data1.qc, QC_cat!="failed")
data1.ok1<-data1.ok[,2:5]

datawide_pos20<-spread(data1.ok1,sampleID,coverage,convert=FALSE) 
no.samp<-ncol(datawide_pos20)-2
colnames(datawide_pos20)[(no.samp+2)] #check, should be last sample
datawide_pos20$gequal4x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 4)
datawide_pos20$gequal8x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 8)
datawide_pos20$gequal16x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 16)
datawide_pos20$gequal24x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 24)

#Calculate coverage thresholds by locus, QC ok samples only ####
totcov_pos20<-datawide_pos20[,c(1,2,(no.samp+3),(no.samp+4),(no.samp+5),(no.samp+6))]
totcov_pos20$prop4x<-totcov_pos20$gequal4x/no.samp
totcov_pos20$prop8x<-totcov_pos20$gequal8x/no.samp
totcov_pos20$prop16x<-totcov_pos20$gequal16x/no.samp
totcov_pos20$prop24x<-totcov_pos20$gequal24x/no.samp

a20<-ggplot(totcov_pos20, aes(x=prop4x)) + 
  geom_histogram(binwidth=0.05,colour="black", fill="skyblue") +
  theme_bw()+xlab("Proportion of QCok samples")+ylab("Loci covered >=4x")+coord_cartesian(xlim=c(-0.05,1),ylim=c(0,800))
a20
b20<-ggplot(totcov_pos20, aes(x=prop8x)) + 
  geom_histogram(binwidth=0.05,colour="black", fill="skyblue") +
  theme_bw()+xlab("Proportion of QCok samples")+ylab("Loci covered >=8x")+coord_cartesian(xlim=c(-0.05,1),ylim=c(0,800))
b20
c20<-ggplot(totcov_pos20, aes(x=prop16x)) + 
  geom_histogram(binwidth=0.05,colour="black", fill="skyblue") +
  theme_bw()+xlab("Proportion of QCok samples")+ylab("Loci covered >=16x")+coord_cartesian(xlim=c(-0.05,1),ylim=c(0,815))
c20
d20<-ggplot(totcov_pos20, aes(x=prop24x)) + 
  geom_histogram(binwidth=0.05,colour="black", fill="skyblue") +
  theme_bw()+xlab("Proportion of QCok samples")+ylab("Loci covered >=24x")+coord_cartesian(xlim=c(-0.05,1),ylim=c(0,815))
d20

pdf("Coverage by locus GENREF hist.QCok.pdf", width=6, height=6)
multiplot(a20,c20,b20,d20,  cols=2) #Figure S2
dev.off()
############################################################

#Coverage pos 20 file by mapped alignment #s:####
samp<-read.csv("Relpos20_GENREF_bysample.allsamp.NEWDATA.csv")
flst<-read.csv("Flagstat_combined_filt_unfilt_NEWDATA.csv")
flst$platewell<-paste(flst$plate,flst$well,sep="")
flst$sampleID<-flst$ID
flst<-flst[,c(33,32,3:31)]
samp_flst<-merge(flst,samp)
QC<-read.csv("No.seq_QC_NEWDATA.csv")
QC$platewell<-paste(QC$plate,QC$well,sep="")
qc.1<-QC[,c(10,1,7:9)]
FLG<-merge(qc.1,samp_flst,by="platewell")  
colnames(FLG)[(ncol(flst))+(ncol(qc.1)-1)] #check-this should be the last column before individual scaffold coverage columns
FLG_short<-FLG[,c(1:((ncol(flst))+(ncol(qc.1)-1)),(ncol(FLG)-3),(ncol(FLG)-2),(ncol(FLG)-1),(ncol(FLG)))]
FLG_short.fail<-subset(FLG_short, QC_cat=="failed")
FLG_short.ok<-subset(FLG_short, QC_cat!="failed")

#Graph the # of Mapped Alignments vs. the No. of Loci that were covered at >=4X####
#The goals of this are to:####
#1. general patterns of what coverage we are getting across samples
#2. see if there's a clear threshold for No. mapped filtered alignments you need to get good enough coverage across targets
#3. is our 10K read threshold generally discerning of what we should call a 'failed' sample?

#take out outlier to be able to visualize rest of patterns better:####
FLG_short1<-subset(FLG_short, DID !=15048)
FLG_short1.fail<-subset(FLG_short1, QC_cat=="failed")
FLG_short1.ok<-subset(FLG_short1, QC_cat!="failed")

all.lined<-ggplot(FLG_short1, aes(x=mapped.reads.QC.passed.y,y=gequal4x)) +  geom_point(alpha=0.4,size=3,fill="blue",shape=21)+theme_bw()+geom_vline(aes(xintercept=25000),linetype="dashed")+geom_vline(aes(xintercept=75000),linetype="dashed")+ylab("loci (relative pos 20) covered at >=4x")+xlab("No. Filtered Alignments")
all.lined
pdf("Coverage pos 20-allsamplesforMS_OTRM_NEWDATA.pdf", width=8, height=6)
all.lined #Figure S3A
dev.off()

#Look at samples that had >10K reads, but <=1/4 loci covered at 4X?####
FLG_short.meh<-subset(FLG_short,gequal4x<=500)
#the relationship between these suggests that our filter of 'failed' at 10K is a good one because although there are some samples with higher #s of reads that aren't getting good recovery, there are some that have only ~12K reads that have 400+ loci covered.
#one option is to add a secondary filter-i.e., 'failed' = <10K reads OR <5K filtered alignments.

#graph to demonstrate this:
meh4<-ggplot(FLG_short.meh, aes(x=mapped.reads.QC.passed.y,y=No_seq,fill=gequal4x)) +  geom_point(alpha=0.8,size=2,shape=21) +theme_bw()+ylab("No. Raw Reads")+xlab("Filtered Alignments")+annotate("text", x = 12000, y = 97000, label = "Samples <=500 loci at 4X ")+ geom_hline(aes(yintercept=10000),linetype="dashed")+ geom_vline(aes(xintercept=5000),linetype="dashed")
pdf("Coverage pos 20.samp<500loci4x_thresholdlines_NEWDATA.pdf", width=8, height=5)
meh4  #Figure S3B
dev.off()

#Extract new QC passed samples:
NF<-subset(FLG_short, QC_cat=="failed" | mapped.reads.QC.passed.y <=5000)
QCdf_short<-QCdf[,c(10:19)]
colnames(QCdf_short)[1]<-"platewell"
NF1<-merge(NF,QCdf_short)
Passed<-subset(FLG_short, QC_cat!="failed" & mapped.reads.QC.passed.y >=5000)
P1<-merge(Passed,QCdf_short)

#Graph with species colored:####
library(RColorBrewer)
FLG_short1$species<-gsub("Ei","EI", FLG_short1$species)
FLG_short1$species<-factor(FLG_short1$species,levels=c("Ccar","CM","EI","hybrid","LK","LV","DC" ))
f<-brewer.pal(8, "Set1")

all<-ggplot(FLG_short1, aes(x=mapped.reads.QC.passed.y,y=gequal4x)) + geom_point(alpha=0.8,size=3,shape=21,aes(fill = species)) +xlab("No. Filtered Alignments")+theme_bw()+ylab("relative position 20 covered >=4x")+scale_fill_brewer(palette="Set1")
#+annotate("text", x = 250000, y = 600, label = "All Samples except outlier ")

FLG_short_noDC<-subset(FLG_short1,species!="DC")

all_noDC<-ggplot(FLG_short_noDC, aes(x=mapped.reads.QC.passed.y,y=gequal4x)) + geom_point(alpha=0.8,size=3,shape=21,aes(fill = species)) +xlab("No. Filtered Alignments")+theme_bw()+ylab("")+scale_fill_brewer(palette="Set1")
#+annotate("text", x = 100000, y = 400, label = "All other species ")

pdf("Coverage pos 20.byspecies_NEWDATA.pdf", width=6, height=8)
multiplot(all,all_noDC,  cols=1) #Figure 3
#this graph shows nicely the fallout of different species at the lower loci #;note the green turtles up high and two hybrids-would be interesting to go back and see if those hybrids are green mixes.
dev.off()
##Samples maxing out at ~1750 despite higher coverage

##### IV. Assessing proportion of reads on target ########
bait_ref1<-read.csv("RAD_Rapture_combined_flagstat.csv")#this is from flagstat combined of RAD and Rapture data that were aligned to the target baits; note: not filtered, just read 1's
bait_ref2<-read.delim("All_Rapture_head_sort_baitref_combined.flagstat_reformatNEWDATA.txt",header=T)
bait_ref2$Rapture_or_RAD<-"Rapture_NEWDATA"
bait_ref3<-rbind(bait_ref1,bait_ref2)
bait_ref3$QC_cat <- ifelse(bait_ref3$total.reads.QC.passed <= 10000, "failed", ifelse((bait_ref3$total.reads.QC.passed >10000) & (bait_ref3$total.reads.QC.passed <= 100000), "good", "high"))
bait_ref3.1<-subset(bait_ref3, sampleID  !="EI_15048_pG09_wA08_baitref_sort_flagstat.txt")
bait_ref3.1<-subset(bait_ref3.1,sampleID !="EI_0_15048_pG09_wA08_sort_flagstat.txt")

p1.1<-ggplot(bait_ref3.1, aes(x=total.reads.QC.passed, y=percent.reads.mapped)) +theme_bw()+
  geom_point(aes(fill=Rapture_or_RAD),alpha=0.7,size=2.5,shape=21)  +  
  ylab("Prop. Mapped of Total Reads") +xlab("Total Reads") +scale_fill_manual(values=c("red","blue","yellow"))+guides(fill=FALSE)
pdf("On-off Target Capture_no_outliertop.pdf", width=6, height=6)
p1.1 #Figure 1A
dev.off()

#V. Effects of Bait GC content ####
data_pos20<-read.table("./All.combined.pos20.baitref.2.bed", header=FALSE)
names(data_pos20)<-c("scaffold", "position","coverage","sampleID")
no.loci<-2007
QC<-read.csv("No.seq_QC_NEWDATA.csv")#Exclude 'failed' samples based on 10K read threshold
QC$platewell<-paste(QC$plate,QC$well,sep="")
qc.1<-QC[,c(8,10)]
z<-colsplit(data_pos20$sampleID,"_",c("species","LABID","DID","plate","well"))
data_pos20$platewell<-paste(z$plate,z$well,sep="")
data_pos20.qc<-merge(data_pos20,qc.1,by="platewell")  
data_pos20.ok<-subset(data_pos20.qc, QC_cat!="failed")
data_pos20.ok1<-data_pos20.ok[,2:5]

datawide_pos20<-spread(data_pos20.ok1,sampleID,coverage,convert=FALSE) 
no.samp<-ncol(datawide_pos20)-2
datawide_pos20$gequal4x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 4)
datawide_pos20$gequal8x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 8)
datawide_pos20$gequal16x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 16)
datawide_pos20$gequal24x<-rowSums(datawide_pos20 [3:(no.samp+2)] >= 24)
datawide_pos20$scaf_pos<-paste(datawide_pos20$scaffold,"_",datawide_pos20$position)

baitAT<-read.csv("../OLD_RaptureDataAnalysis//Bait_ATcontent.csv")
dw.short_pos20.2<-datawide_pos20[,c((no.samp+7),(no.samp+3),(no.samp+4),(no.samp+5),(no.samp+6))]
dw.short_pos20.2$Bait<-gsub("_Priority*. _ 20", "", dw.short_pos20.2$scaf_pos)
cov<-dw.short_pos20.2[,c(6,2:5)]
df2<-merge(baitAT,cov,by="Bait")

pdf("GC content by coverage_newdata.pdf", width=6, height=6)
p1<-ggplot(df2, aes(x=GC_content,y=gequal8x)) +  geom_point(fill="blue",alpha=0.6,size=2,shape=21) +ylab(">=8x coverage")+theme_bw()+xlab("Bait GC content")
p1
dev.off()
#repeated with old data to combine for Figure S1 a&b
