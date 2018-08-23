#### Set up ####
setwd("~/...")
library(tidyr)
library(plyr)
library(reshape2)
library(ggplot2)

#multiplot function ####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
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

#Across Populations:####
#RF Pr2 only ####
setwd("./theta_RFPR2_acrosspops/")
pestPG<-list.files(path=".")
for (i in 1:length(pestPG)){
  file<-pestPG[i]
  sp<-paste0(strsplit(file, "_")[[1]][1])
  assign(sp, read.table(pestPG[i], header=F))
}

colnames(Ccar) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
Ccar$species<-"Ccar"
colnames(Cm) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
Cm$species<-"CM"
colnames(EI) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
EI$species<-"EI"
colnames(LK) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
LK$species<-"LK"
colnames(LV) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
LV$species<-"LV"
colnames(DC) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DC$species<-"DC"

comb<-rbind(Ccar,Cm,EI,LK,LV,DC)

#Watterson's Estimator:
mean_RFPr2<-tapply(comb$Watt_Est,comb$species,mean)
n_RFPr2<-tapply(comb$Watt_Est,comb$species,length)
median_RFPr2<-tapply(comb$Watt_Est,comb$species,median)

mean_all<-as.data.frame(tapply(comb$Watt_Est,comb$species,mean))
names(mean_all)[1] <- "mean_Watt_Est"

#Pi:
mean_PiRFPr2<-tapply(comb$Pi,comb$species,mean)
n_PiRFPr2<-tapply(comb$Pi,comb$species,length)
median_PiRFPr2<-tapply(comb$Pi,comb$species,median)

#boxplots
w<-ggplot(comb, aes(x=species, y=Watt_Est)) + 
  geom_violin()
w + geom_boxplot(width=0.05)

p<-ggplot(comb, aes(x=species, y=Pi)) + 
  geom_violin()
p + geom_boxplot(width=0.05)


#Within Representative Populations:####
#RF Pr2 only ####
setwd("../theta_RFPR2_withinpops/")
list.files()
pestPG<-list.files(path=".")
for (i in 1:length(pestPG)){
  file<-pestPG[i]
  sp<-paste0(strsplit(file, "_")[[1]][1])
  assign(sp, read.table(pestPG[i], header=F))
}

colnames(CcarJapan) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
CcarJapan$species_pop<-"Ccar_Japan"
colnames(CmAmSam) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
CmAmSam$species_pop<-"Cm_AmSam"
colnames(CmCR) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
CmCR$species_pop<-"Cm_CR"
colnames(CmRMI) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
CmRMI$species_pop<-"Cm_RMI"
colnames(DCBraz) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DCBraz$species_pop<-"DC_Braz"
colnames(DCCR) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DCCR$species_pop<-"DC_CR"
colnames(DCGab) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DCGab$species_pop<-"DC_Gab"
colnames(DCIndo) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DCIndo$species_pop<-"DC_Indo"
colnames(DCMex) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
DCMex$species_pop<-"DC_Mex"
colnames(EICR) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
EICR$species_pop<-"EI_CR"
colnames(EIHI) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
EIHI$species_pop<-"EI_HI"
colnames(EINic) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
EINic$species_pop<-"EI_Nic"
colnames(LK) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
LK$species_pop<-"LK"
colnames(LV) <- c("Pos_stop_starts",	"Chr",	"WinCenter",	"Watt_Est","Pi",	"tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")
LV$species_pop<-"LV"

comb<-rbind(CcarJapan,CmAmSam,CmCR,CmRMI,DCBraz,DCCR,DCGab,
            DCIndo,DCMex,EICR,EIHI,EINic,LK,LV)
comb$species_pop<-factor(comb$species_pop,c("Ccar_Japan","Cm_AmSam","Cm_RMI","Cm_CR",
                                            "EI_CR","EI_Nic","EI_HI","LK","LV", "DC_Braz","DC_Gab",
                                            "DC_CR","DC_Indo","DC_Mex"))   

comb$species<-comb$species_pop
comb<-comb %>% separate(species, c("species", "pop"))
comb<-comb[,1:16]

#Watterson's Estimator:
mean_RFPr2<-tapply(comb$Watt_Est,comb$species_pop,mean)
n_RFPr2<-tapply(comb$Watt_Est,comb$species_pop,length)
median_RFPr2<-tapply(comb$Watt_Est,comb$species_pop,median)

mean_all<-as.data.frame(tapply(comb$Watt_Est,comb$species.pop,mean))
names(mean_all)[1] <- "mean_Watt_Est"

#Pi:
mean_PiRFPr2<-tapply(comb$Pi,comb$species_pop,mean)
n_PiRFPr2<-tapply(comb$Pi,comb$species_pop,length)
median_PiRFPr2<-tapply(comb$Pi,comb$species_pop,median)

#boxplots
col<-c("darkorange2","forestgreen","firebrick","gold1","blue2","darkviolet")
w<-ggplot(comb, aes(x=species_pop, y=Watt_Est, fill=species)) + 
  geom_violin(alpha=0.5)
S<-w + geom_boxplot(width=0.05)+theme_bw()+ylab("S")+theme(axis.title.x=element_blank(),
                                                          axis.text.x=element_blank(),
                                                          axis.ticks.x=element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(values=col)

p<-ggplot(comb, aes(x=species_pop, y=Pi, fill=species)) + 
  geom_violin(alpha=0.5)
Pi<-p + geom_boxplot(width=0.05)+theme_bw()+ylab("Pi")+
  xlab("Species & Population Group")+theme(legend.position="none")+
  scale_fill_manual(values=col)
multiplot(S,Pi,  cols=1)

pdf("Thetas_Representative_Pops.pdf", width=12, height=6)
multiplot(S,Pi,  cols=1) #Figure 6
dev.off()