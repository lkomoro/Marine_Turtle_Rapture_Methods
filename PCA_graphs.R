#PCA custom plots using PC scores generated in ANGSD/NGSCovar
#-------------------#####
#READ IN 
setwd("~/.../")
library(ggplot2)
library(reshape2)
#----
#Read in annotation file####
annot<-read.csv("../NEWDATAmetadata_groupings.csv")
str(annot)
#list.files()
#--------------
#20K subsampled ####
#Global hardshell PCAs####
PC_scores<-read.table("HS_20-nocall_KQCP_PC_scores.txt",header=TRUE)
IDs<-read.table("HS_20KQCP.clst",header=TRUE)
IDs$FID<-gsub("_sortfltr_20000","",IDs$FID)
PC_scoresIDs<-cbind(IDs,PC_scores)
PC_scoresIDs<-PC_scoresIDs[,c(1,4:6)]
z<-colsplit(PC_scoresIDs$FID,"_",c("species","LABID","DID","plate","well"))
PC_scoresIDs$platewell<-paste0(z$plate,z$well)
m<-merge(PC_scoresIDs,annot,by="platewell") 
g1<-ggplot(m, aes(x=PC1,y=PC2)) +  geom_point(aes(fill=species),alpha=0.7,size=4,shape=21) +xlab("PC1")+theme_bw()+ylab("PC2")
g1.1<-g1+scale_fill_hue(l=50) #Figure 4A
g2<-ggplot(m, aes(x=PC1,y=PC3)) +  geom_point(aes(fill=species),alpha=0.8,size=3,shape=21) +xlab("PC1")+theme_bw()+ylab("PC3")
g2.1<-g2+scale_fill_hue(l=50)
g3<-ggplot(m, aes(x=PC3,y=PC2)) +  geom_point(aes(fill=species),alpha=0.8,size=3,shape=21) +xlab("PC3")+theme_bw()+ylab("PC2")
g3.1<-g3+scale_fill_hue(l=50)
pdf("HS_20k-NOCALL_global_PCA_NEWDATA.pdf", width=6, height=8)
multiplot(g1.1,g2.1,g3.1,  cols=1)
dev.off()

#Lo/Lv and Hybrid Seperation####
PC_scores<-read.table("Hybcheck_nocall_20KQCP_PC_scores.txt",header=TRUE)
IDs<-read.table("Hybcheck_nocall_20KQCP.clst",header=TRUE)
IDs$FID<-gsub("_sortfltr_20000","",IDs$FID)
PC_scoresIDs<-cbind(IDs,PC_scores)
PC_scoresIDs<-PC_scoresIDs[,c(1,4:6)]
z<-colsplit(PC_scoresIDs$FID,"_",c("species","LABID","DID","plate","well"))
PC_scoresIDs$platewell<-paste0(z$plate,z$well)
m<-merge(PC_scoresIDs,annot,by="platewell") 
m<-subset(m,LABID!=37285 & LABID!=104313 & LABID!=49124  & LABID!=49122)

g1<-ggplot(m, aes(x=PC1,y=PC2)) +  geom_point(aes(fill=species),alpha=0.7,size=4,shape=21) +xlab("PC1")+theme_bw()+ylab("PC2")
g1.1<-g1+scale_fill_hue(l=50) #Figure 4B
g2<-ggplot(m, aes(x=PC1,y=PC3)) +  geom_point(aes(fill=species),alpha=0.8,size=3,shape=21) +xlab("PC1")+theme_bw()+ylab("PC3")
g2.1<-g2+scale_fill_hue(l=50)
g3<-ggplot(m, aes(x=PC3,y=PC2)) +  geom_point(aes(fill=species),alpha=0.8,size=3,shape=21) +xlab("PC3")+theme_bw()+ylab("PC2")
g3.1<-g3+scale_fill_hue(l=50)
pdf("Hybcheck_20K_undetHybs_nocall_PCA_NEWDATA.pdf", width=6, height=8)
multiplot(g1.1,g2.1,g3.1,  cols=1)
dev.off()

########################################
#40K subsampled####
#Global green turtle PCAs
PC_scores<-read.table("CM_nocall_40kQCP_PC_scores.txt",header=TRUE)
IDs<-read.table("CM_nocall_40kQCP.clst",header=TRUE)
IDs$FID<-gsub("_sortfltr_40000","",IDs$FID)
PC_scoresIDs<-cbind(IDs,PC_scores)
PC_scoresIDs<-PC_scoresIDs[,c(1,4:6)]
z<-colsplit(PC_scoresIDs$FID,"_",c("species","LABID","DID","plate","well"))
PC_scoresIDs$platewell<-paste0(z$plate,z$well)
m<-merge(PC_scoresIDs,annot,by="platewell") 
str(m)
g1<-ggplot(m, aes(x=PC1,y=PC2)) +  geom_point(aes(fill=Ocean_Region),alpha=0.8,size=2,shape=21) +xlab("PC1")+theme_bw()+ylab("PC2")
g1 #Figure S5A
g2<-ggplot(m, aes(x=PC1,y=PC3)) +  geom_point(aes(fill=Ocean_Region),alpha=0.8,size=2,shape=21) +xlab("PC1")+theme_bw()+ylab("PC3")
g2
g3<-ggplot(m, aes(x=PC3,y=PC2)) +  geom_point(aes(fill=Ocean_Region),alpha=0.8,size=2,shape=21) +xlab("PC3")+theme_bw()+ylab("PC2")
g3
pdf("CM_40k_global_nocall_PCA_NEWDATA.pdf", width=6, height=8)
multiplot(g1,g2,g3,  cols=1)
dev.off()
#--------------
#Western Pacific C. mydas PCA####
PC_scores<-read.table("wpac_CM_nocall_40kQCP_PC_scores.txt",header=TRUE)
IDs<-read.table("wpac_CM_nocall_40kQCP.clst",header=TRUE)
IDs$FID<-gsub("_sortfltr_40000","",IDs$FID)
PC_scoresIDs<-cbind(IDs,PC_scores)
PC_scoresIDs<-PC_scoresIDs[,c(1,4:6)]
z<-colsplit(PC_scoresIDs$FID,"_",c("species","LABID","DID","plate","well"))
PC_scoresIDs$platewell<-paste0(z$plate,z$well)
m<-merge(PC_scoresIDs,annot,by="platewell") 
g1<-ggplot(m, aes(x=PC1,y=PC2)) +  geom_point(aes(fill=Loc_combined),alpha=0.8,size=2,shape=21) +xlab("PC1")+theme_bw()+ylab("PC2")
g1 #Figure S5B
g2<-ggplot(m, aes(x=PC1,y=PC3)) +  geom_point(aes(fill=Loc_combined),alpha=0.8,size=2,shape=21) +xlab("PC1")+theme_bw()+ylab("PC3")
g2
g3<-ggplot(m, aes(x=PC3,y=PC2)) +  geom_point(aes(fill=Loc_combined),alpha=0.8,size=2,shape=21) +xlab("PC3")+theme_bw()+ylab("PC2")
g3
pdf("CM_40k_Wpac_nocall_PCA_NEWDATA.pdf", width=6, height=8)
multiplot(g1,g2,g3,  cols=1)
dev.off()
