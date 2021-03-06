library(openxlsx)
library(GenomicRanges)
library(IRanges)

source("cnvAnalysis_Group_plotSummary.R")

# Before running this script, save files containing genomic segments and called CN status ("balanced", "gain" or "loss" for each segment. 
# This can be done using R package conumee.
# For this in github.com/dstichel/conumee 3 different versions of baseline correction are avaiable:
# MAD: Centering the baseline according to median average deviation of intensities of all samples
# MAXDENS: Centering the baseline to the intensitiy level of the most probes in the sample
# BAF: Finding diploid level based on evaluation of SNP probes included in the chip

anno.genome.size0<-read.xlsx("anno.genome.size.xlsx",colNames=FALSE)
anno.genome.pq0<-read.xlsx("anno.genome.pq.xlsx",colNames=FALSE)
anno.genome.chr0<-read.xlsx("anno.genome.chr.xlsx"),colNames=FALSE)
chromosome_positions<-read.xlsx("chromosome_positions.xlsx"))
chromosome_positions$chromosome<-paste0("chr",chromosome_positions$chromosome)

path.out <- ...  # path for output

allCases<-read.xlsx(...) #data.frame containing columns:
                           #"ID" (=sample identifier)
                           #"group" (=group annotation)
                           #"filename.for.segment.data" 

#Plot für alle samples
allCases$group.for.plot<-"all"
groups <-split(allCases,allCases$group.for.plot)
plotCNASummary(groups,path.out,TRUE,TRUE)
      
#Plot für einzelne Gruppen
allCases$group.for.plot<-allCases$group
groups <-split(allCases,allCases$group)
plotCNASummary(groups,path.out,TRUE,TRUE)


