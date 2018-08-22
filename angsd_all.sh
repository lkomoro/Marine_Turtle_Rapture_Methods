#!/bin/bash -l

nInd=96

minInd=1

angsd -bam bamlist_all -out all -minQ 20 -minMapQ 10 -minInd $minInd -GL 1 -doMajorMinor 1 -doMaf 2


