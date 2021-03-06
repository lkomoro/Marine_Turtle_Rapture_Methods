setwd("~/...")
getwd()

ls()
list.files()
Pac<-read.table("TURTLE-Priority1-probes-120-filtration.txt", header=TRUE)
str(Pac)
Pac.stringent<-subset(Pac,Stringent=="pass")
Pac.s.no.RM<-subset(Pac.stringent, X.RM==0)
Pac.s.no.RM_GC<-subset(Pac.s.no.RM,X.GC<=65)
Pac.s.no.RM_GC<-subset(Pac.s.no.RM_GC,X.GC>=31)
write.csv(Pac.s.no.RM_GC,file="Priority_1_final.csv")

All<-read.table("TURTLE-Priority2-probes-120-filtration.txt", header=TRUE)
str(All)
All.stringent<-subset(All,Stringent=="pass")
All.s.no.RM<-subset(All.stringent, X.RM==0)
All.s.no.RM_GC<-subset(All.s.no.RM,X.GC<=65)
All.s.no.RM_GC<-subset(All.s.no.RM_GC,X.GC>=35)
write.csv(All.s.no.RM_GC,file="Priority_2_final.csv")
