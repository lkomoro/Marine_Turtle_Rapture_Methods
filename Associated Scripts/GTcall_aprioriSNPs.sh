#!/bin/bash -l

ulimit -n 1500
angsd -bam AllQCpassedsamples_bamlist_NEWDATA.txt -out QCPA_0.9_kSNP -minQ 20 -minMapQ 10  -GL 1 -doMajorMinor 1 -doMaf 2 -doGeno 4 -doPost 1 -postCutoff 0.9 -rf ./Rapture_sitekey_2017-14-02_forANGSD_SNPpos_targets.txt

angsd -bam AllQCpassedsamples_bamlist_NEWDATA.txt -out QCPA_0.8_kSNP -minQ 20 -minMapQ 10  -GL 1 -doMajorMinor 1 -doMaf 2 -doGeno 4 -doPost 1 -postCutoff 0.8 -rf ./Rapture_sitekey_2017-14-02_forANGSD_SNPpos_targets.txt

angsd -bam AllQCpassedsamples_bamlist_NEWDATA.txt -out QCPA_0.95_kSNP -minQ 20 -minMapQ 10  -GL 1 -doMajorMinor 1 -doMaf 2 -doGeno 4 -doPost 1 -postCutoff 0.95 -rf ./Rapture_sitekey_2017-14-02_forANGSD_SNPpos_targets.txt
