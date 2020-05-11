plotCNASummary <-function(groups,outpath,remove.artifacts,plot.title=TRUE){

lapply(groups,function(x) {
  groupname<-x$group.for.plot[1]
  no.of.cases<-length(x$ID)
  message(paste0("Now evaluating group ", groupname, " with ", no.of.cases, " cases."))

  files<-x$filename.for.segment.data
  file.colnames <- c(colnames(read.csv(files[1],sep="",header=TRUE)),"origin")
  
  seg <- do.call(rbind,lapply(files, function(file) {
    if(file.exists(file)){
      act.file <- read.csv(file,sep="",header=TRUE)
      act.file$origin<-file
      colnames(act.file) <- file.colnames
      act.file
    }else {print("File doesn't exist:")                  
      print(file)}
  }))
  
  message(paste0("Chromothripsis found: ", length(which(seg$alteration=="chromothripsis"))))
  if (length(which(seg$alteration=="chromothripsis"))>0){
    seg.chromothripsis<-seg[which(seg$alteration=="chromothripsis"),]
    seg<-seg[-which(seg$alteration=="chromothripsis"),]
    write.xlsx(seg.chromothripsis, paste(outpath,"CNA_chromothripsis_detected_",groupname,"_",no.of.cases,".xlsx")) 
  }

 colnames(seg)
    seg<-seg[,which(colnames(seg)%in%c("chrom","loc.start","loc.end","alteration"))]
    colnames(seg)
    newrow1p<-c("chr1",121392717,144500871,"balanced")
    newrow9p<-c("chr9",44817471,66534425,"balanced")
    newrow13p<-c("chr13",100000,17500000,"balanced")
    newrow14p<-c("chr14",100000,17500000,"balanced")
    newrow15p<-c("chr15",100000,18500000,"balanced")
    newrow16p<-c("chr16",35017911,46517891,"balanced")
    newrow21p1<-c("chr21",100000,10848939,"balanced")  
    newrow21p2<-c("chr21",11109342,14407229,"balanced")   
    newrow22p<-c("chr22",100000,14500000,"balanced")
    seg<-rbind(seg,newrow1p,newrow9p,newrow13p,newrow14p,newrow15p,newrow16p,newrow21p1,newrow21p2,newrow22p)                          
    
    GR<- GRanges(seqnames=seg$chrom,ranges=IRanges(as.numeric(seg$loc.start),as.numeric(seg$loc.end)))
    disj<-as.data.frame(GenomicRanges::disjoin(GR))
    
  count <- lapply(1:nrow(disj),function(i){
    x <- disj[i,]
    act.seq <- seg[seg$chrom == x$seqnames & as.numeric(seg$loc.start) <= as.numeric(x$start) & as.numeric(seg$loc.end) >= as.numeric(x$end),]
    balanced<-sum(act.seq$alteration =="balanced")
    gain <- sum(act.seq$alteration =="gain")+sum(act.seq$alteration =="amp")
    loss <- sum(act.seq$alteration =="loss")+sum(act.seq$alteration =="homodel")
    list(gain,loss,balanced)
  })
  
  disj$gains=unlist(lapply(count,function(x)x [1]))*100/length(files)    # / Anzahl Fälle in Gruppe
  disj$loss=unlist(lapply(count,function(x)x [2]))*100/length(files)
  disj$balanced=unlist(lapply(count,function(x)x [3]))*100/length(files)  
  
if (remove.artifacts){
  #!!! Remove peaks, which turned out to be artifacts
  disj<-removeBadData(disj,"chr5",140149999,140825000)    # Protocadherins?
  disj<-removeBadData(disj,"chr8",6924999,8100001)        # Defensines (immune...)
  disj<-removeBadData(disj,"chr11",50049999,56675000)     # Olfactory receptors
  disj<-removeBadData(disj,"chr12",34374999,38178347)     # 12p11.1, da ist gar nix
  disj<-removeBadData(disj,"chr14",17499999,20450001)     # Olfactory receptors
  disj<-removeBadData(disj,"chr16",32074999,33986576)     # ???
  }

  disj2<-disj[rep(1:nrow(disj),each=2),]
  odd_indexes<-seq(1,nrow(disj2),2)
  even_indexes<-seq(2,nrow(disj2),2)                        
  disj2$xpos<-NA
  disj2$xpos[odd_indexes]<-disj2$start[odd_indexes]
  disj2$xpos[even_indexes]<-disj2$end[even_indexes]
  
  disj2$start<-as.numeric(disj2$start)
  disj2$end<-as.numeric(disj2$end)
  disj2$width<-as.numeric(disj2$width)
  disj2$xpos<-as.numeric(disj2$xpos)
  disj2$loss<-as.numeric(disj2$loss)
  disj2$gains<-as.numeric(disj2$gains)
  disj2$balanced<-as.numeric(disj2$balanced)
  
  disj2$chromnum<-gsub("chr","",disj2$seqnames)
  disj2$chromnum<-gsub("X","23",disj2$chromnum)
  disj2$chromnum<-gsub("Y","24",disj2$chromnum)
  disj2$chromnum<-as.numeric(disj2$chromnum)
  disj2<-disj2[order(disj2$chromnum,as.numeric(disj2$xpos)),]
  
  newrow1<-c("chr1",10,20,10,"*",0,0,0,10,1)
  newrowY<-c(as.character(disj2$seqnames[nrow(disj2)]),as.numeric(disj2$end[nrow(disj2)])+10,as.numeric(disj2$end[nrow(disj2)])+20,10,"*",0,0,0,as.numeric(disj2$end[nrow(disj2)])+10,ifelse(as.character(disj2$seqnames[nrow(disj2)])=="chr22",22,ifelse(as.character(disj2$seqnames[nrow(disj2)])=="chrX","chrX","chrY")))
  disj2<-rbind(disj2,newrow1,newrowY)
  disj2$start<-as.numeric(disj2$start)
  disj2$end<-as.numeric(disj2$end)
  disj2$xpos<-as.numeric(disj2$xpos)
  disj2$loss<-as.numeric(disj2$loss)
  disj2$gains<-as.numeric(disj2$gains)
  disj2$balanced<-as.numeric(disj2$balanced)
  disj2$chromnum<-gsub("chr","",disj2$seqnames)
  disj2$chromnum<-gsub("X","23",disj2$chromnum)
  disj2$chromnum<-gsub("Y","24",disj2$chromnum)
  disj2$chromnum<-as.numeric(disj2$chromnum)
  disj2<-disj2[order(disj2$chromnum,as.numeric(disj2$xpos)),]
  
  if ("chrX" %in% disj2$seqnames){disj2<-disj2[-which(disj2$seqnames=="chrX"),]}
  if ("chrY" %in% disj2$seqnames){disj2<-disj2[-which(disj2$seqnames=="chrY"),]} 
  disj2<-rbind(disj2,disj2[nrow(disj2),] )
  disj2[nrow(disj2),c(6:9)]<-c(0,0,0,  disj2[nrow(disj2),9]+1)
  anno.genome.size<-anno.genome.size0[1:22,,drop=FALSE]
  anno.genome.pq<-anno.genome.pq0[1:22,,drop=FALSE]
  anno.genome.chr<-anno.genome.chr0[1:22,,drop=FALSE]

    
  
message("plot.title=", plot.title)
  if (plot.title){
   main.title<-paste(groupname, " (n=", no.of.cases, ")",sep="")
  }else{main.title<-""}
  
  pdf(paste(outpath,groupname,"_",no.of.cases,".pdf", sep=""), width = 18, height = 6)
  PLOTa<-dev.cur()
  png(paste(outpath,groupname,"_",no.of.cases,".png", sep=""), width = 1800, height = 600)
  dev.control("enable")
  par(mfrow = c(1, 1), mgp=c(4,1,0),mar = c(4, 8, 4, 4), oma = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, sum(as.numeric(anno.genome.size$X1))), ylim = c(-94, 94), xaxs = "i", xaxt = "n", yaxt = "n",
       xlab = NA, ylab = "rate of alterations [%]", main = main.title,cex=3,cex.lab=3,cex.axis=3,cex.main=3)
  abline(v = cumsum(as.numeric(head(anno.genome.size$X1, -1))),  col = "black")
  abline(v = cumsum(c(0, head(anno.genome.size$X1, -1))) + anno.genome.pq$X1, col = "black", lty = 2)
  abline(h = 0, col="black")

  axis(1, at = (cumsum(c(0, head(anno.genome.size$X1, -1))) + anno.genome.size$X1/2), labels = gsub("chr","",anno.genome.chr$X1), las = 1, col.axis="black",cex.axis=1)
  axis(2, las = 2, lty=3, at = seq(0, 100, 20), labels=abs(seq(0, 100, 20)),col.axis="green",cex.axis=2)
  axis(2, las = 2, lty=3, at = seq(-100, 0, 20), labels=abs(seq(-100, 0, 20)),col.axis="red",cex.axis=2)
  axis(2, las = 2, lty=3, at = seq(0, 0, 20), labels=abs(seq(0, 0, 20)),col.axis="black",cex.axis=2)
  polygon(cumsum(c(0, head(anno.genome.size$X1, -1)))[match(disj2$seqnames, anno.genome.chr$X1)] + disj2$xpos, disj2$gains,col="grey",lwd=3)
  polygon(cumsum(c(0, head(anno.genome.size$X1, -1)))[match(disj2$seqnames, anno.genome.chr$X1)] + disj2$xpos, -disj2$loss,col="grey",lwd=3)
  dev.copy(which=PLOTa)
  dev.off() 
  dev.off()
  
})  #end of lapply on groups

} #end of function



removeBadData<-function(x,chrom,start,end){
  bad.data.index<-which(x$seqnames==chrom & as.numeric(x$start)>start&as.numeric(x$end)<end)
  
  if(length(bad.data.index)>0){
    if (x$seqnames[min(bad.data.index)-1]==chrom){
      x$end[min(bad.data.index)-1]<-x$end[max(bad.data.index)]
      x<-x[-bad.data.index,]
  }else{
    x$start[max(bad.data.index)+1]<-x$start[min(bad.data.index)]
    x<-x[-bad.data.index,]
  }
  }
  return(x)}
        
