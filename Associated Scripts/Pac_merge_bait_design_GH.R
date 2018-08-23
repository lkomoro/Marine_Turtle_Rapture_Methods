#playing around with Rapture design lists
setwd("~/...")
getwd()
pac_site84<-read.table("Pac_bait_sites_84.txt", header=TRUE)
pac_site84_92adj<-read.table("Pac_bait_sites_84-92adj.txt", header=TRUE)
pac_site100<-read.table("Pac_bait_sites_100.txt", header=TRUE)
pac_site150<-read.table("Pac_bait_sites_150.txt", header=TRUE)

#concantenate first
pac_site84$SCAFFOLD_SNP_POSITION<-paste(pac_site84$SCAFFOLD,"_",pac_site84$SNP_POSITION)
pac_site84_92adj$SCAFFOLD_SNP_POSITION<-paste(pac_site84_92adj$SCAFFOLD,"_",pac_site84_92adj$SNP_POSITION)
pac_site100$SCAFFOLD_SNP_POSITION<-paste(pac_site100$SCAFFOLD,"_",pac_site100$SNP_POSITION)
pac_site150$SCAFFOLD_SNP_POSITION<-paste(pac_site150$SCAFFOLD,"_",pac_site150$SNP_POSITION)

#quick syntax check
str(pac_site84_92adj)
#merge to figure out common/unique ones:
pac_site84and84_92adj<-merge(pac_site84,pac_site84_92adj,"SCAFFOLD_SNP_POSITION",all=TRUE)
pac_site84and84_92adj_100<-merge(pac_site84and84_92adj,pac_site100,"SCAFFOLD_SNP_POSITION",all=TRUE)
pac_site84and84_92adj_100_150<-merge(pac_site84and84_92adj_100,pac_site150,"SCAFFOLD_SNP_POSITION",all=TRUE)

description="Pac.bait_candidates_merged"
write.csv(pac_site84and84_92adj_100_150, file = paste(description, ".csv", sep = ""))

#All group:
setwd("~/Desktop/leatherbackRAD work Davis March 2016/All_bait_candidate_files")
All_site84<-read.table("All_sites_84.txt", header=TRUE)
All_site84_92adj<-read.table("All_sites_84-92.txt", header=TRUE)
All_site100<-read.table("All_sites_100.txt", header=TRUE)
All_site150<-read.table("All_sites_150.txt", header=TRUE)

#concantenate first
All_site84$SCAFFOLD_SNP_POSITION<-paste(All_site84$SCAFFOLD_84,"_",All_site84$SNP_POSITION_84)
All_site84_92adj$SCAFFOLD_SNP_POSITION<-paste(All_site84_92adj$SCAFFOLD_84.92,"_",All_site84_92adj$SNP_POSITION_84.92)
All_site100$SCAFFOLD_SNP_POSITION<-paste(All_site100$SCAFFOLD_100,"_",All_site100$SNP_POSITION_100)
All_site150$SCAFFOLD_SNP_POSITION<-paste(All_site150$SCAFFOLD_150,"_",All_site150$SNP_POSITION_150)

#quick syntax check
str(All_site84_92adj)
#merge to figure out common/unique ones:
All_site84and84_92adj<-merge(All_site84,All_site84_92adj,"SCAFFOLD_SNP_POSITION",all=TRUE)
All_site84and84_92adj_100<-merge(All_site84and84_92adj,All_site100,"SCAFFOLD_SNP_POSITION",all=TRUE)
All_site84and84_92adj_100_150<-merge(All_site84and84_92adj_100,All_site150,"SCAFFOLD_SNP_POSITION",all=TRUE)

description="All.bait_candidates_merged"
write.csv(All_site84and84_92adj_100_150, file = paste(description, ".csv", sep = ""))

