library(openxlsx)
library(GenomicRanges)
library(IRanges)

source("cnvAnalysis_Group_plotSummary.R")
#source("removeArtifacts.R")

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


