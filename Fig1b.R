#  This script generates the final version of Fig 1b
#  BICC - Nicolas Guex, 2023.09.07

    setwd('/Users/nguex/Documents/_Collaborations/JustineCollier/Revisions')
    x <- read.table('all_GANTC_fused_WG.txt',header=T)
    c <- read.table('all_coverage_fused_GC.txt',header=T)
    a.all <- read.csv('fused.csv',header=TRUE,row.names=1)
    genes <- read.table('Atu_degs_Mut-WT_with_GANTC.51.fused.txt',sep="\t",header=T)
    genes$Chromosome <- gsub('.[0-9]$', '', genes$Chromosome)
    ori <- read.table('ori_fused.txt',sep="\t",header=T)
    ori <- cbind(ori,Gene_Center = ori$Gene_Start+0.5*(ori$Gene_End-ori$Gene_Start))
    ori$Gene_ID <- c("ter2L", "ter2R", "ter1", "ori1", "gcrA", "ccrM", "ori2", "RepC", "RepC")

    plotConfig  <- data.frame(colname=c('JC2140_neg_ipd','JC2307_pos_ipd','JC2307_neg_ipd'),color=c('#5DBFC2','#F5A460','#CE5C5C'),ofsetY=c(-1,+1,+1))

    # loess computed out of the plotting loop to do it once and allow playing with plotting parameters without recomputing.
    whichChr <- 1
    data <- x[which(x$chr==rownames(a.all)[whichChr]),]

    loess.predict <- list()
    for (i in 1:3)
    {
        ipd_colidx <- which(colnames(data) == plotConfig$colname[i])
        x1 <- data[which(x$chr==rownames(a.all)[whichChr]),'pos']
        x1Pred <- sort(c(x1,ori$Gene_Center[1:7]))
        y1 <- data[,ipd_colidx]
        loess.fit = loess(y1 ~ x1)
        loess.predict[[i]] = predict(loess.fit, x1Pred, se = TRUE)
    }


    mm2inch <- 0.03937008
    pdf('Fig1b.pdf',width=mm2inch*180,height=mm2inch*90,pointsize=7)
    # plot frame
    plot(0,type='n',main='',xlim=c(a.all[whichChr,'start'],a.all[whichChr,'end']),ylim=c(0,5.5),las=1,ylab='GANTC IPD',xlab='dicentric chromosome position [Mbp]',frame=FALSE,xaxt='n')
    axis(1, at = seq(0,5000000,by=1000000), round(seq(0,5,by=1),digits=1), tick = TRUE, line = NA,pos = NA, outer = FALSE, font = NA, lty = "solid",lwd = 1, lwd.ticks = 1, col = NULL, col.ticks = NULL,hadj = NA, padj = NA, gap.axis = NA)
    abline(h=1:5,lty=2,col='grey')
    for (i in 1:7)
    {
        rect(ori$Gene_Start[i], 0, ori$Gene_End[i], 5, border = "grey",col="grey")
        text(ori$Gene_Start[i]+0.5*(ori$Gene_End[i]-ori$Gene_Start[i]),5.25,ori$Gene_ID[i],cex=1,font=3)
    }

    # plot lines
    plotConfig  <- data.frame(colname=c('JC2140: WT -IPTG','JC22307: DccrM Ptac-ccrM +IPTG','JC22307: DccrM Ptac-ccrM -IPTG'),color=c('#5DBFC2','#F5A460','#CE5C5C'),ofsetY=c(-1,+1,+1))
    for (i in 1:3)
    {
        p.fit <- loess.predict[[i]]$fit
        p.se  <- loess.predict[[i]]$se.fit
        p.df  <- loess.predict[[i]]$df
        print(qt(0.975,p.df))
        d1 = c(x1Pred, rev(x1Pred))
        d2 = c(    p.fit + qt(0.975,p.df) * p.se,
               rev(p.fit - qt(0.975,p.df) * p.se))
        polygon(d1, d2, col = plotConfig$color[i], border = NA)
        lines(x1Pred, p.fit)

        for (g in 1:7)
        {
            xpos <- ori$Gene_Start[g]+0.5*(ori$Gene_End[g]-ori$Gene_Start[g])
            ypos <- p.fit[which(x1Pred == xpos)]
            if (is.na(ypos)) { ypos <- p.fit[2] } # use leftmost data.
            rect(ori$Gene_Start[g], 0, ori$Gene_End[g], 5, border = "grey",col="grey")
            text(xpos+80000,ypos+plotConfig$ofsetY[i]*0.25,round(ypos,digits=1),cex=0.75,col=plotConfig$color[i],font=2)
        }
    }
    
    l1 <- 'JC2140: WT -IPTG'
    l2 <- expression('JC2307: '*Delta*italic('ccrM',1,2)~'P'*italic('tac-ccrM')~' +IPTG')
    l3 <- expression('JC2307: '*Delta*italic('ccrM',1,2)~'P'*italic('tac-ccrM')~' -IPTG')
    legend(50000,0.75,fill=plotConfig$color,legend=c(l1,l2,l3),cex=0.66)
    dev.off()
    

# ---------------------------------------------------------------------------------------------------------------------

