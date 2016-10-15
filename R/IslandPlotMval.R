#' Mvalue Ploting
#'
#' Plot multiple CpG M-values for an island,shore,and shelf.
#' @param dat contains a dataframe with illumina 450k methylation data
#' @param IslName name of CpG island you want to graph
#' @param Clockset dataframe including info on CpG clocks on Islands
#' @keywords M_val islands
#' @export
#' @examples
#' IslandPlot_Mval(beM.64, "chr7:77648451-77648668",  "CpGs.hypo.drift10.isl")
# plots M-values for different people at different CpGs
IslandPlot_Mval = function(dat, IslName, ClockSet, bottomleft = "") {

  gofish = as.character(IslName)
  cpgn.GENE = manifestData[manifestData$Islands_Name==gofish,"Name"]
  gene.GENE = manifestData[manifestData$Islands_Name==gofish,"UCSC_RefGene_Name"]
  pos.GENE = manifestData[manifestData$Islands_Name==gofish,"pos"]
  isl.GENE = manifestData[manifestData$Islands_Name==gofish,"Relation_to_Island"]
  refgrp.GENE = manifestData[manifestData$Islands_Name==gofish,"UCSC_RefGene_Group"]

  np = ncol(dat)
  ido = order(pos.GENE)
  cpgn = cpgn.GENE[ido]
  labels = refgrp.GENE[ido]

  ncpg = length(cpgn)

  # colors
  cols = rep(1,ncpg)
  cols[isl.GENE[ido]=="Island"]=2
  cols[grepl("Shore",isl.GENE[ido])]=6
  cols[grepl("Shelf",isl.GENE[ido])]=4

  ids = na.omit(match(cpgn,rownames(dat)))
  inc=(1:ncpg)[!is.na(match(cpgn,rownames(dat)))]
  plot(1:ncpg,rep(10,ncpg),pch="",ylim=c(-8,8),xaxt='n',xlab="CpG",ylab="M-value")
  axis(1,at=c(1:ncpg),labels=FALSE)
  text(x=c(1:ncpg),y=rep(par("usr")[3],9), labels = labels, srt = 45, adj = c(1.2,2), xpd = TRUE, cex=.5)
  for(i in 1:length(inc)) {points(rep(inc[i],np),dat[ids[i],],col=cols[inc[i]],pch=19,cex=.3)}

  ## gene memberships; indicate range by CpGs
  genes = unique(unlist(strsplit(gene.GENE,";")))
  for(i in 1:length(genes)) {
    i1 = i-1
    idum1 = (1:ncpg)[grepl(genes[i],gene.GENE[ido])]
    idum2 = (1:ncpg)[grepl("Island", isl.GENE[ido])]
    lines(sort(idum1),rep(7.5-0.8*i1,length(idum1)),lwd=10,col='darkgrey')
    lines(sort(idum2),rep(7.5,length(idum2)),lwd=5,col='red')
    text(mean(idum1),6.8-0.8*i1,genes[i],cex=1)
  }

  ## indicate clock CpGs
  idum = na.omit(match(ClockSet,cpgn))
  ids = na.omit(match(cpgn[idum],rownames(dat)))
  points(idum,rep(7.5,length(idum)),pch='|',cex=0.7)

  # for(i in 1:length(idum)) {points(rep(idum[i],np),dat[ids[i],],col="orange",pch=19,cex=.2)}
  legend('bottomleft',bottomleft,bty='n',y.intersp=1.25,cex=1.2)

  return(list(cpgn=cpgn,genes=genes))
}
