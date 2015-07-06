#RSVSim.R
# Simulate 1129 deletions, 490 insertions, 202 inversions and 145 tandem duplications
# between 500bp and 10kb on hg19.fa

library(RSVSim)

SVminSize=50
SVmaxSize=10000

delNum=2943

dupNum=532
max_dups=3

insNum=503

inv_minSize=500
invNum=24

delSizes = estimateSVSizes(n=delNum, minSize=SVminSize, maxSize=SVmaxSize, default="deletions", hist=FALSE)
dupSizes = estimateSVSizes(n=dupNum, minSize=SVminSize, maxSize=SVmaxSize, default="tandemDuplications", hist=FALSE)
insSizes = estimateSVSizes(n=insNum, minSize=SVminSize, maxSize=SVmaxSize, default="insertions", hist=FALSE)
invSizes = estimateSVSizes(n=invNum, minSize=inv_minSize, maxSize=SVmaxSize, default="inversions", hist=FALSE)

setwd("/home/tjiang/")
sim = simulateSV(output="./Rtest/", genome="init_hg19.fa", dels=delNum, sizeDels=delSizes, ins=insNum, sizeIns=insSizes, 
                 invs=invNum, sizeInvs=invSizes, dups=dupNum, sizeDups=dupSizes, maxDups = max_dups,
                 percCopiedIns=1, repeatBias=TRUE, random=c(TRUE, TRUE,TRUE,TRUE, TRUE), verbose=TRUE)
