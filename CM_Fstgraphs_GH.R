setwd("~/...")
list.files()
library(tidyr)
library(gplots)
library(RColorBrewer)
library(plyr)
library(reshape2)

#old datasheet
#data<-read.delim("CM_pwFst_list1.txt", header=T)#self groups removed
#data_weight<-data[,c(1,2,4)]
#data_w_wide<-data_weight %>% spread(Population2,Fst.Weight)

d<-read.csv("CM_pwFst_list_MArch2018.csv",header=T)
levels(d$Loc1)
d$Loc1<-factor(d$Loc1, levels=c( "Cyprus","Gabon","STX","CostaRica" , "Mexico" ,"Galapagos","FFS" ,"Japan","AmericanSamoa" , "Malaysia"  ,"FrenchPolynesia","FSM","RMI","Australia"))
e<-as.data.frame(acast(d, Loc1~Loc2, value.var="weighted_Fst"))
f<-e[,c(4,8,14,3,12,9,5,10,1,11,6,7,13,2)]
#write.csv(f,"CM_weighted_Fst_March2018.csv")
data_w_wide1<-d %>% spread(Loc2,weighted_Fst)
data_w_wide1[data_w_wide1< 0] <- NA
data_w_wide1<-data_w_wide1[,c(1,5,9,15,4,13,10,6,11,2,12,7,8,14,3)]
data_w_wide1[is.na(data_w_wide1)] <- 100 #set self crosses to number so not NA issue
#heatmap setup:####
real.data2<-data_w_wide1
rnames <- real.data2[,1] 
mat_real.data2 <- data.matrix(real.data2[,2:ncol(real.data2)])
rownames(mat_real.data2) <- rnames
breaks = c(seq(0,0.1,by=0.01), 
           seq(0.101,0.5,by=0.01)
           )
my_palette <- colorRampPalette(c("aliceblue","deepskyblue3","blue","blue4"))(n=length(breaks)-1)
pdf("CM_weightedFst_March2018.pdf", width=10, height=8)
heatmap.2(mat_real.data2, 
          main = "Weighted Fst", 
          notecol="black",
          density.info="histogram",
          trace="none",
          margins =c(8,8),
          col=my_palette,
          breaks=breaks,
          key=TRUE,
          symkey=FALSE,
          scale="none",
          cexRow=1, cexCol=1,
          dendrogram="none",
          Colv="NA",
          Rowv='NA',
         na.color="blue"
)
dev.off()



