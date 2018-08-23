#fin <- commandArgs(T)

setwd("~/...")
#list.files()
library(ape)
library(phangorn)

#regular trees:
fin <-"./spcompV4.tree"
fout <- paste(fin,".pdf",sep="",collaps="")

pdf(file=fout)

plot(read.tree(fin), cex=0.5)

invisible(dev.off())

#plotTreeBoots.R:
 # fin <- commandArgs(T)
fin <-"./spcompV4.boot.tree"
fout <- paste(fin,".pdf",sep="",collaps="")

#library(ape)
#library(phangorn)

trees <- read.tree(fin, skip=2)
tree <- read.tree(fin)[[1]]

pdf(file=fout)

plotBS(tree, trees, type="phylo", cex=0.5, p=10) #change p argument here if want more support values plotted-default is supposed to be only <80 but 55 shows up

invisible(dev.off())
