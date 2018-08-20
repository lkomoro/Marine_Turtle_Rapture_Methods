setwd("/Users/lkomoroske/Box/Research Projects/NGS Projects and resources/Pacific leatherback rapture project/Dcor_Rapturedata/NEW_RaptureDataAnalysis/NGS_admix_files")
library(ggplot2)
library(tidyr)
library(reshape2)
annot<-read.csv("../NEWDATAmetadata_groupings.csv")

#Figure 4C:#### 
IDs<-read.table("./hardshell_QCP.bamlist",header=FALSE)
IDs$bamlistorder<-as.numeric(rownames(IDs))
IDs$V1<-gsub("_sortfltr.bam","",IDs$V1)
z<-colsplit(IDs$V1,"_",c("species","LABID","DID","plate","well"))
IDs$platewell<-paste0(z$plate,z$well)
m<-merge(IDs,annot,by="platewell") 
str(m)
m[m=="EI"]="Ei"
pop<-m[,c(1,3,5,6,7,22,14:16)]
pop1<-pop[order(pop$bamlistorder),]
#K=6 #N.B. Ran across range of K and compared across to determine best fit for visualization
d<-read.table("HS_K6.qopt")
d$platewell<-(pop1[,1])
d$order<-(pop1[,2])
d1<-gather(d,"K","proportion",1:6)
d2<-merge(d1,pop1)
d2$species<-factor(d2$species)
d2$ID<-factor(d2$platewell, levels=unique(d2$platewell[order(d2$species)]), ordered=TRUE)
d2$species<-factor(d2$species, levels=c("CM" ,"Ei","hybrid", "Ccar","LK","LV"))
#d2<-subset(d2, species !="LV" & species !="LK")
d2$IDc<-paste0(d2$Loc_combined,d2$order)
d2$facetkey<-as.factor(paste0(d2$species,"_",d2$Ocean_Region))
d2$facetkey<-factor(d2$facetkey, levels=c("CM_Med", "CM_Eatl",  "CM_Watl", "CM_Wpac", "CM_Cpac", "CM_Epac", "Ei_Watl","Ei_Wind" ,"Ei_Wpac","Ei_Cpac",  "Ei_Epac", "Ccar_Udet","Ccar_Watl", "Ccar_Eatl","Ccar_Wpac","Ccar_Upac","Ccar_Cpac", "hybrid_Udet","LK_Udet","LV_Udet" ))

p2 <- ggplot() + geom_bar(aes(y = proportion, x = IDc, fill = K), data = d2,stat="identity",width=1)+  scale_fill_manual(values=c("blue","gray","yellow","red","forestgreen","skyblue"))
p3<-p2+facet_grid(.~facetkey, scales="free_x",space="free")+
  theme(panel.spacing.x = unit(0, "lines"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6))
pdf("HS_alldata_NGSadmixK6_NEWDATA.pdf", width=12, height=6)
p3
dev.off()

#Extract Kemp's/Hybrids for new SNP disco bamlist:####
d3<-subset(pop1,bamlistorder>=1&bamlistorder<=20)
d3<-subset(d3,bamlistorder!=11&bamlistorder!=12&bamlistorder!=13&bamlistorder!=14)
d3.1<-merge(d3,IDs)
d3.2<-paste0(d3.1$V1,"_sortfltr.bam")
#write.table(d3.2,"ridley-hybridcheck.bamlist",row.names = F,col.names = F,quote = F)
d3.2<-paste0(d3.1$V1,"_sortfltr_20000.bam")#note a couple have dropped out-just removed manually to do quickly for now
#write.table(d3.2,"ridley-hybridcheck20k.bamlist",row.names = F,col.names = F,quote = F)

#########################################
#Figure 4D to examine Cm/Ei influence/weights:
IDs<-read.table("./HS_noEiorCM_QCP.bamlist",header=FALSE)
IDs$bamlistorder<-as.numeric(rownames(IDs))
IDs$V1<-gsub("_sortfltr.bam","",IDs$V1)
z<-colsplit(IDs$V1,"_",c("species","LABID","DID","plate","well"))
IDs$platewell<-paste0(z$plate,z$well)
m<-merge(IDs,annot,by="platewell") 
pop<-m[,c(1,3,5,6,7,22,14:16)]
pop1<-pop[order(pop$bamlistorder),]

#K=4
d<-read.table("Hybcheck_KQCP_K4.qopt")
d$platewell<-(pop1[,1])
d$order<-(pop1[,2])
d1<-gather(d,"K","proportion",1:4)
d2<-merge(d1,pop1)
d2$bamlistorder<-factor(d2$bamlistorder)
d2<-subset(d2, bamlistorder!=31 &bamlistorder!=32 &bamlistorder!=33& bamlistorder!=30)
p2 <- ggplot() + geom_bar(aes(y = proportion, x = bamlistorder, fill = K), data = d2,stat="identity",width=1)+  scale_fill_manual(values=c("blue","yellow","red","green"))
p3<-p2+facet_grid(.~species, scales="free_x",space="free")+
  theme(panel.spacing.x = unit(0, "lines"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("hybcheck_NGSadmixK4_NEWDATA.pdf", width=8, height=6)
p3
dev.off()

