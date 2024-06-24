library(circlize)
#setwd('E:/VMShared/Justine/')
setwd('/data1/chris/Justine/')
options("scipen"=6)
mycex <- 0.3
mycexlegend <- 0.5
axisFontSize <- 0.375
glabelpos <- 1.0
labelXpos <- 0.5
labelYpos <- 5.8
t.min <- 0.3
t.max <- 0.8
pad <- 2000
firstTrackMax <- t.max

x <- read.table('all_GANTC_fused.txt',header=T)
c <- read.table('all_coverage_fused_GC.txt',header=T)
a.all <- read.csv('fused.csv',header=TRUE,row.names=1)
genes <- read.table('Atu_degs_Mut-WT_with_GANTC.51.fused.txt',sep="\t",header=T)
genes$Chromosome <- gsub('.[0-9]$', '', genes$Chromosome)
ori <- read.table('ori_fused.txt',sep="\t",header=T)

motifs <- c( "GAATC +", "GATTC -", "GATTC +", "GAATC -", "GAGTC +", "GACTC -", "GACTC +", "GAGTC -")
h <- rep(seq(0,0.95,by = 1.0/length(motifs)))
s <- rep(c(1.0,0.8),length(motifs))
v <- rep(c(1.0,0.8),length(motifs))
colors <- hsv(h,s,v)

iAp <- which(x$motif=='GAATC' & x$strand=='+')
iAm <- which(x$motif=='GAATC' & x$strand=='-')
iCp <- which(x$motif=='GACTC' & x$strand=='+')
iCm <- which(x$motif=='GACTC' & x$strand=='-')
iGp <- which(x$motif=='GAGTC' & x$strand=='+')
iGm <- which(x$motif=='GAGTC' & x$strand=='-')
iTp <- which(x$motif=='GATTC' & x$strand=='+')
iTm <- which(x$motif=='GATTC' & x$strand=='-')

pdf('GANTC_circos_fused.pdf',paper="a4r",width=0,height=0)
op <- par(mar=c(1,0,1,0))

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = c(1,1,3))
circos.initialize(xlim = a.all)
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iAp], x=x$pos[iAp], y=x$JC2307_neg_ipd[iAp]/x$JC2140_neg_ipd[iAp],pch=1,cex=mycex,col=colors[1])
circos.trackPoints(x$chr[iTm], x=x$pos[iTm], y=x$JC2307_neg_ipd[iTm]/x$JC2140_neg_ipd[iTm],pch=4,cex=mycex,col=colors[2])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iTp], x=x$pos[iTp], y=x$JC2307_neg_ipd[iTp]/x$JC2140_neg_ipd[iTp],pch=1,cex=mycex,col=colors[3])
circos.trackPoints(x$chr[iAm], x=x$pos[iAm], y=x$JC2307_neg_ipd[iAm]/x$JC2140_neg_ipd[iAm],pch=4,cex=mycex,col=colors[4])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iGp], x=x$pos[iGp], y=x$JC2307_neg_ipd[iGp]/x$JC2140_neg_ipd[iGp],pch=1,cex=mycex,col=colors[5])
circos.trackPoints(x$chr[iCm], x=x$pos[iCm], y=x$JC2307_neg_ipd[iCm]/x$JC2140_neg_ipd[iCm],pch=4,cex=mycex,col=colors[6])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iCp], x=x$pos[iCp], y=x$JC2307_neg_ipd[iCp]/x$JC2140_neg_ipd[iCp],pch=1,cex=mycex,col=colors[7])
circos.trackPoints(x$chr[iGm], x=x$pos[iGm], y=x$JC2307_neg_ipd[iGm]/x$JC2140_neg_ipd[iGm],pch=4,cex=mycex,col=colors[8])

legend(-1,.93,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
       title='IPD of JC2307_neg vs JC2140_neg on fused C58',
       fill=colors[seq(1,length(motifs))],
       legend=motifs)
circos.yaxis(side='left',track.index=1,sector.index=rownames(a.all)[1],labels.cex=axisFontSize)

for(i in 1:nrow(a.all))
{
	n <- rownames(a.all)[i]
	circos.axis(sector.index = n, track.index=1, major.at = seq(0,a.all[i,2],by=100000), labels.cex=axisFontSize, labels.facing='clockwise', direction='outside')
	circos.text(a.all[i,2]/2, firstTrackMax, n, track.index = 1, sector.index = n, facing = "inside", cex = axisFontSize + 0.25,adj = c(labelXpos, -labelYpos))
	gi <- which(genes$Chromosome == n & genes$Strand == "+")
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 1, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 2, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 3, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 4, sector.index = n, col='pink', border='pink')
	#circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),paste(genes$Gene_ID,as.character(genes$Fold_Change)),track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize,adj=c(0.5,-0.5),col='grey30')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
	p95 <- mean(x$JC2307_neg_ipd/x$JC2140_neg_ipd,na.rm=T)+1.959964*sd(x$JC2307_neg_ipd/x$JC2140_neg_ipd,na.rm=T)
	p99 <- mean(x$JC2307_neg_ipd/x$JC2140_neg_ipd,na.rm=T)+2.575829*sd(x$JC2307_neg_ipd/x$JC2140_neg_ipd,na.rm=T)
	circos.lines(c(1,a.all[i,2]), c(p95,p95), track.index = 1, sector.index = n, col='grey')
	circos.lines(c(1,a.all[i,2]), c(p99,p99), track.index = 1, sector.index = n, col='pink')
	circos.lines(c(1,a.all[i,2]), c(p95,p95), track.index = 2, sector.index = n, col='grey')
	circos.lines(c(1,a.all[i,2]), c(p99,p99), track.index = 2, sector.index = n, col='pink')
	circos.lines(c(1,a.all[i,2]), c(p95,p95), track.index = 3, sector.index = n, col='grey')
	circos.lines(c(1,a.all[i,2]), c(p99,p99), track.index = 3, sector.index = n, col='pink')
	circos.lines(c(1,a.all[i,2]), c(p95,p95), track.index = 4, sector.index = n, col='grey')
	circos.lines(c(1,a.all[i,2]), c(p99,p99), track.index = 4, sector.index = n, col='pink')
	gi <- which(genes$Chromosome == n & genes$Strand == "-")
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 1, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 2, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 3, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max,genes$Gene_End[gi]+pad,t.max+0.05, track.index = 4, sector.index = n, col='cyan', border='cyan')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
	go <- which(ori$Chromosome == n)
	circos.rect(ori$Gene_Start[go]-pad,t.max,ori$Gene_End[go]+pad,t.max+0.05, track.index = 1, sector.index = n, col='green', border='green')
	circos.rect(ori$Gene_Start[go]-pad,t.max,ori$Gene_End[go]+pad,t.max+0.05, track.index = 2, sector.index = n, col='green', border='green')
	circos.rect(ori$Gene_Start[go]-pad,t.max,ori$Gene_End[go]+pad,t.max+0.05, track.index = 3, sector.index = n, col='green', border='green')
	circos.rect(ori$Gene_Start[go]-pad,t.max,ori$Gene_End[go]+pad,t.max+0.05, track.index = 4, sector.index = n, col='green', border='green')
	circos.text(ori$Gene_Start[go],rep(glabelpos,length(go)),ori$Gene_ID,
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
}

circos.clear()

t.min <- 1.0
t.max <- 7.0
firstTrackMax <- t.max
glabelpos <- 9.0

x <- read.table('all_GANTC_fused_WG.txt',header=T)

iAp <- which(x$motif=='GAATC' & x$strand=='+')
iAm <- which(x$motif=='GAATC' & x$strand=='-')
iCp <- which(x$motif=='GACTC' & x$strand=='+')
iCm <- which(x$motif=='GACTC' & x$strand=='-')
iGp <- which(x$motif=='GAGTC' & x$strand=='+')
iGm <- which(x$motif=='GAGTC' & x$strand=='-')
iTp <- which(x$motif=='GATTC' & x$strand=='+')
iTm <- which(x$motif=='GATTC' & x$strand=='-')

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = c(1,1,3))
circos.initialize(xlim = a.all)
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iAp], x=x$pos[iAp], y=x$JC2140_neg_ipd[iAp],pch=1,cex=mycex,col=colors[1])
circos.trackPoints(x$chr[iTm], x=x$pos[iTm], y=x$JC2140_neg_ipd[iTm],pch=4,cex=mycex,col=colors[2])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iTp], x=x$pos[iTp], y=x$JC2140_neg_ipd[iTp],pch=1,cex=mycex,col=colors[3])
circos.trackPoints(x$chr[iAm], x=x$pos[iAm], y=x$JC2140_neg_ipd[iAm],pch=4,cex=mycex,col=colors[4])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iGp], x=x$pos[iGp], y=x$JC2140_neg_ipd[iGp],pch=1,cex=mycex,col=colors[5])
circos.trackPoints(x$chr[iCm], x=x$pos[iCm], y=x$JC2140_neg_ipd[iCm],pch=4,cex=mycex,col=colors[6])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iCp], x=x$pos[iCp], y=x$JC2140_neg_ipd[iCp],pch=1,cex=mycex,col=colors[7])
circos.trackPoints(x$chr[iGm], x=x$pos[iGm], y=x$JC2140_neg_ipd[iGm],pch=4,cex=mycex,col=colors[8])

legend(-1,.93,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
       title='IPD of JC2140_neg on fused C58',
       fill=colors[seq(1,length(motifs))],
       legend=motifs)
circos.yaxis(side='left',track.index=1,sector.index=rownames(a.all)[1],labels.cex=axisFontSize)

for(i in 1:nrow(a.all))
{
	n <- rownames(a.all)[i]
	circos.axis(sector.index = n, track.index=1, major.at = seq(0,a.all[i,2],by=100000), labels.cex=axisFontSize, labels.facing='clockwise', direction='outside')
	circos.text(a.all[i,2]/2, firstTrackMax, n, track.index = 1, sector.index = n, facing = "inside", cex = axisFontSize + 0.25,adj = c(labelXpos, -labelYpos))
	gi <- which(genes$Chromosome == n & genes$Strand == "+")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='pink', border='pink')
	#circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),paste(genes$Gene_ID,as.character(genes$Fold_Change)),track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize,adj=c(0.5,-0.5),col='grey30')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
	gi <- which(genes$Chromosome == n & genes$Strand == "-")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='cyan', border='cyan')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
}

circos.clear()

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = c(1,1,3))
circos.initialize(xlim = a.all)
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iAp], x=x$pos[iAp], y=x$JC2307_neg_ipd[iAp],pch=1,cex=mycex,col=colors[1])
circos.trackPoints(x$chr[iTm], x=x$pos[iTm], y=x$JC2307_neg_ipd[iTm],pch=4,cex=mycex,col=colors[2])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iTp], x=x$pos[iTp], y=x$JC2307_neg_ipd[iTp],pch=1,cex=mycex,col=colors[3])
circos.trackPoints(x$chr[iAm], x=x$pos[iAm], y=x$JC2307_neg_ipd[iAm],pch=4,cex=mycex,col=colors[4])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iGp], x=x$pos[iGp], y=x$JC2307_neg_ipd[iGp],pch=1,cex=mycex,col=colors[5])
circos.trackPoints(x$chr[iCm], x=x$pos[iCm], y=x$JC2307_neg_ipd[iCm],pch=4,cex=mycex,col=colors[6])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iCp], x=x$pos[iCp], y=x$JC2307_neg_ipd[iCp],pch=1,cex=mycex,col=colors[7])
circos.trackPoints(x$chr[iGm], x=x$pos[iGm], y=x$JC2307_neg_ipd[iGm],pch=4,cex=mycex,col=colors[8])

legend(-1,.93,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
       title='IPD of JC2307_neg on fused C58',
       fill=colors[seq(1,length(motifs))],
       legend=motifs)
circos.yaxis(side='left',track.index=1,sector.index=rownames(a.all)[1],labels.cex=axisFontSize)

for(i in 1:nrow(a.all))
{
	n <- rownames(a.all)[i]
	circos.axis(sector.index = n, track.index=1, major.at = seq(0,a.all[i,2],by=100000), labels.cex=axisFontSize, labels.facing='clockwise', direction='outside')
	circos.text(a.all[i,2]/2, firstTrackMax, n, track.index = 1, sector.index = n, facing = "inside", cex = axisFontSize + 0.25,adj = c(labelXpos, -labelYpos))
	gi <- which(genes$Chromosome == n & genes$Strand == "+")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='pink', border='pink')
	#circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),paste(genes$Gene_ID,as.character(genes$Fold_Change)),track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize,adj=c(0.5,-0.5),col='grey30')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
	gi <- which(genes$Chromosome == n & genes$Strand == "-")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='cyan', border='cyan')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
}

circos.clear()

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = c(1,1,3))
circos.initialize(xlim = a.all)
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iAp], x=x$pos[iAp], y=x$JC2307_pos_ipd[iAp],pch=1,cex=mycex,col=colors[1])
circos.trackPoints(x$chr[iTm], x=x$pos[iTm], y=x$JC2307_pos_ipd[iTm],pch=4,cex=mycex,col=colors[2])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iTp], x=x$pos[iTp], y=x$JC2307_pos_ipd[iTp],pch=1,cex=mycex,col=colors[3])
circos.trackPoints(x$chr[iAm], x=x$pos[iAm], y=x$JC2307_pos_ipd[iAm],pch=4,cex=mycex,col=colors[4])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iGp], x=x$pos[iGp], y=x$JC2307_pos_ipd[iGp],pch=1,cex=mycex,col=colors[5])
circos.trackPoints(x$chr[iCm], x=x$pos[iCm], y=x$JC2307_pos_ipd[iCm],pch=4,cex=mycex,col=colors[6])
circos.track(ylim = c(t.min,t.max), track.height = 1.0/6)
circos.trackPoints(x$chr[iCp], x=x$pos[iCp], y=x$JC2307_pos_ipd[iCp],pch=1,cex=mycex,col=colors[7])
circos.trackPoints(x$chr[iGm], x=x$pos[iGm], y=x$JC2307_pos_ipd[iGm],pch=4,cex=mycex,col=colors[8])

legend(-1,.93,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
       title='IPD of JC2307_pos on fused C58',
       fill=colors[seq(1,length(motifs))],
       legend=motifs)
circos.yaxis(side='left',track.index=1,sector.index=rownames(a.all)[1],labels.cex=axisFontSize)

for(i in 1:nrow(a.all))
{
	n <- rownames(a.all)[i]
	circos.axis(sector.index = n, track.index=1, major.at = seq(0,a.all[i,2],by=100000), labels.cex=axisFontSize, labels.facing='clockwise', direction='outside')
	circos.text(a.all[i,2]/2, firstTrackMax, n, track.index = 1, sector.index = n, facing = "inside", cex = axisFontSize + 0.25,adj = c(labelXpos, -labelYpos))
	gi <- which(genes$Chromosome == n & genes$Strand == "+")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='pink', border='pink')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='pink', border='pink')
	#circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),paste(genes$Gene_ID,as.character(genes$Fold_Change)),track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize,adj=c(0.5,-0.5),col='grey30')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
	gi <- which(genes$Chromosome == n & genes$Strand == "-")
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 1, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 2, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 3, sector.index = n, col='cyan', border='cyan')
	circos.rect(genes$Gene_Start[gi]-pad,t.max+0.15,genes$Gene_End[gi]+pad,t.max+0.65, track.index = 4, sector.index = n, col='cyan', border='cyan')
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
	            paste(genes$Gene_ID,format(round(genes$Fold_Change,2),nsmall=2)),
	            track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
}

circos.clear()

par(op)
op <- par(mar=c(4,4,4,4)+0.1)
for(i in 1:nrow(a.all))
{
  n <- rownames(a.all)[i]
  ci <- which(c$chr == n)
  plot(c$start[ci],c$JC2140_neg[ci]/median(c$JC2140_neg[ci])+0.25,type="l",ylim=c(0.0,1.8),main=paste("Coverage/median(Coverage)",n))
  points(c$start[ci],c$JC2307_pos[ci]/median(c$JC2307_pos[ci])+0.5,type="l",col="green")
  points(c$start[ci],c$JC2307_neg[ci]/median(c$JC2307_neg[ci]),type="l",col="red")
  points(c$start[ci],c$GC[ci],type="l",col="magenta")
  abline(h=1.25)
  abline(h=1.5,col="green")
  abline(h=1.0,col="red")
  if(i == 1)
  {
    # could grab coordinates from ori table
    abline(v=2197137)
    abline(v=3867344)
  }
  legend(a.all[i,2]/5,.1,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
         title='Coverage & GC',fill=c("black","red","green","magenta"),legend=c("JC2140_neg + 0.25","JC2307_neg","JC2307_pos + 0.5","%GC"))
}

for(i in 1:nrow(a.all))
{
  n <- rownames(a.all)[i]
  xi <- which(x$chr == n)
  ci <- which(c$chr == n)
  plot(x$pos[xi],sqrt(x$JC2140_neg_ipd[xi]/mean(x$JC2140_neg_ipd[xi]))+0.3,type="l",ylim=c(0.4,2.0),main=paste("sqrt(IPD/mean(IPD))",n))
  points(x$pos[xi],sqrt(x$JC2307_neg_ipd[xi]/mean(x$JC2307_neg_ipd[xi])),type="l",col="red")
  points(x$pos[xi],sqrt(x$JC2307_pos_ipd[xi]/mean(x$JC2307_pos_ipd[xi]))+0.6,type="l",col="green")
  points(c$start[ci],c$GC[ci],type="l",col="magenta")
  abline(h=1.3)
  abline(h=1.0,col="red")
  abline(h=1.6,col="green")
  if(i == 1)
  {
    # could grab coordinates from ori table
    abline(v=2197137)
    abline(v=3867344)
  }
  legend(a.all[i,2]/5,.425,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
         title='IPD & GC',fill=c("black","red","green","magenta"),legend=c("JC2140_neg + 0.3","JC2307_neg","JC2307_pos + 0.6","%GC"))
}

dev.off()

#----------------------------------------------------------------------------------
