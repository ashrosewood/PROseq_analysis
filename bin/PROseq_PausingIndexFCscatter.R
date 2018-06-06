args <- commandArgs()

help <- function(){
    cat("PROseq_PausingIndexFCscatter.R :
- Compare pausing index or body promoter occupancy ratio fold changes for 2 conditions relative to a control.
- Gives a scatter plot for all possible comparisons of sample fold changes and a table with the fold changes. \n")
    cat("Usage: \n")
    cat("--Table   : path to table of averge body and tss coverages                           [required]
                     columns sample names followed by _pi or _prr, last column is gene_id.            \n")
    cat("--Control : Name of Control sample to calculate FC over                               [optional]\n")
    cat("--PRR     : Plot the promoter release ratio (body/promoter) instead of pi (0/1)      [default = 0 plot pausing index]\n")
    cat("--FC      : Fold change to highlight genes in each quadrant                          [default=2]\n")

    cat("--outName : Prefix for plot name and table. Each comparison appended to the end.     [required]
                     (ex. /home/ashley/Pol2_293T_rep1)                                                 \n")
    cat("--q1Col   : Quadrant 1 color (defaults are RdBu colors )                             [optional]\n")
    cat("--q2Col   : Quadrant 2 color                                                         [optional]\n")
    cat("--q3Col   : Quadrant 3 color                                                         [optional]\n")
    cat("--q4Col   : Quadrant 4 color                                                         [optional]\n")
    cat("--plotPcc : Flag to add pearson correlation (0/1)                                    [default=1]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args))){
    help()
} else {
    Table   <- sub( '--Table=', '', args[grep('--Table=', args)])
    Control <- sub( '--Control=', '', args[grep('--Control=', args)])
    outName <- sub( '--outName=', '',args[grep('--outName=',args)])
    PRR     <- sub( '--PRR=', '',args[grep('--PRR=',args)])
    FC      <- sub( '--FC=', '',args[grep('--FC=',args)])
    q1Col   <- sub( '--q1Col=', '',args[grep('--q1Col=',args)])
    q2Col   <- sub( '--q2Col=', '',args[grep('--q2Col=',args)])
    q3Col   <- sub( '--q3Col=', '',args[grep('--q3Col=',args)])
    q4Col   <- sub( '--q4Col=', '',args[grep('--q4Col=',args)])
    plotPcc <- sub( '--plotPcc=', '',args[grep('--plotPcc=',args)])
}

if (identical(PRR,character(0))){
    PRR     <- 0
}else{
    PRR     <- as.numeric(PRR)
}

if (identical(plotPcc,character(0))){
    plotPcc <- 1
}else{
    plotPcc <- as.numeric(plotPcc)
}

if (identical(FC,character(0))){
    FC      <- 2
}else{
    FC      <- as.numeric(FC)
}

if (identical(q1Col,character(0))){
    q1Col   <- "#67001F"
}
if (identical(q2Col,character(0))){
    q2Col   <- "#B2182B"
}
if (identical(q3Col,character(0))){
    q3Col   <- "#053061"
}
if (identical(q4Col,character(0))){
    q4Col   <- "#2166AC"
}

cat("Table:" ,    Table ,    sep="\n")
cat("outName:" ,  outName ,  sep="\n")
cat("Control:" ,  Control ,  sep="\n")
cat("FC:" ,       FC ,       sep="\n")
cat("PRR:" ,      PRR ,      sep="\n")

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

df           <- read.table(Table, sep="\t", header=TRUE)
rownames(df) <- paste(df$gene_id)
dim(df)

Tss          <- df[ , grep("tss$", names(df)) ]
Body         <- df[ , grep("body$", names(df)) ]

## be sure they are in the same order
Order        <- sub( "_tss", "", names(Tss) )
Tss.or       <- Tss[ , paste(Order,"tss", sep="_") ]
Body.or      <- Body[ , paste(Order,"body", sep="_") ]

stopifnot( sub("_tss", "", colnames(Tss.or)) == sub("_body", "", colnames(Body.or)))
stopifnot( rownames(Tss.or) == rownames(Body.or) )

## put the control first and set experiments
exps <- Order[grep(Control, Order, invert=TRUE)]
exps

con  <- Order[grep(Control, Order, invert=FALSE)]
con

if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))(500)[seq(1,500,length.out=length(exps))]

if(PRR==0){
    ## pausing index
    print("Calculating ratio of promoter over body")
    Ratio <- do.call(cbind, lapply(Order, function(x){
        pseudo       <- min( min(Tss[Tss!=0]) , min(Body[Body!=0]) )
        pi           <- (Tss.or[, paste0(x, "_tss")]+pseudo) / (Body.or[, paste0(x, "_body")]+pseudo) 
        pi           <- as.data.frame(pi)
        rownames(pi) <- rownames(Tss.or)
        colnames(pi) <- x
        pi
    }
    ))    
    fc <- do.call(cbind, lapply(1:length(exps), function(i){
        print(paste("FC sample", i, exps[i]))
        FC           <- Ratio[,exps[i]] / Ratio[,con]
        FC           <- as.data.frame(FC)
        rownames(FC) <- rownames(Ratio)
        colnames(FC) <- sub("$", "_fc",exps[i])
        FC
    }))
    for (i in 1:length(exps)) {
        for (j in 1:length(exps)){
            q1 <- log2(fc[,i]) > log2(FC) & log2(fc[,j]) > log2(FC)
            q2 <- log2(fc[,j]) > log2(FC) & log2(fc[,i]) < -log2(FC)
            q3 <- log2(fc[,i]) < -log2(FC) & log2(fc[,j]) < -log2(FC)
            q4 <- log2(fc[,i]) > log2(FC) & log2(fc[,j]) < -log2(FC)
            png(paste(outName, names(fc)[i], "vs", names(fc)[j], "PausingIndexFC.png", sep="."), 500, 500)
            par(mar=c(5,5,5,5))
            plot(log2(fc[,i]), log2(fc[,j]),
                 xlab=paste("log2(",names(fc[i]),")",sep=""),
                 ylab=paste("log2(",names(fc[j]),")",sep=""),
                 main="Pausing Index (Tss/Body)\nFold Change",
                 cex.axis=1.1,
                 cex.lab=1.1,
                 col="#C0C0C070",
                 pch=19)
            if(plotPcc==1){
                legend("topleft",legend= round(cor(fc[,i],fc[,j],method="pearson"), 4), cex=1.1)
            }
            points(log2(fc[ q1, i]), log2(fc[q1, j]), col=q1Col, pch=19,cex=1)
            points(log2(fc[ q3, i]), log2(fc[q3, j]), col=q3Col, pch=19,cex=1)
            points(log2(fc[ q2, i]), log2(fc[q2, j]), col=q2Col, pch=19,cex=1)
            points(log2(fc[ q4, i]), log2(fc[q4, j]), col=q4Col, pch=19,cex=1)
            legend("bottomright", title=paste("FC:", FC)
                  ,legend=c(paste("q1",sum(q1))
                           ,paste("q2",sum(q2))
                           ,paste("q3",sum(q3))
                           ,paste("q4",sum(q4))
                           )
                 ,pch=c(19,19,19,19), col=c(q1Col,q2Col,q3Col,q4Col),cex=0.9,)
            abline(h=0, lty="dashed", col="black")
            abline(v=0, lty="dashed", col="black")
            dev.off()
        }
    }
    Ratio$gene_id <- rownames(Ratio)
    fc$gene_id    <- rownames(fc)
    write.table(fc, file=paste(outName, "PausingIndexFC.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(Ratio, file=paste(outName, "PausingIndex.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)    
    ## violin plot
    df.long <- melt(fc[,grep("gene_id", names(fc), invert=TRUE)])
    df.long$variable <- sub("_fc$", "", df.long$variable)
    pdf(paste(outName, "violinPlotPausingIndexFC.pdf", sep="."), 6, 5)
    print({
    p <-
        ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
        geom_violin() +
        geom_boxplot(width=0.1, colour = "black", outlier.colour=NA, notch=TRUE) +
        theme_classic() +
        scale_fill_manual(values=cols) +
        ylab("log2(Pausing Index FC)")+
        xlab("")+
        ggtitle(basename(outName))+
        theme(panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank(), 
              panel.background = element_blank()
             ,axis.line = element_line(colour = "black"))+
        theme(text = element_text(size=12),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size=11),
              plot.title=element_text(size=11))
    })
    dev.off()
    ## boxPlot
    pdf(paste(outName, "boxPlotPausingIndexFC.pdf", sep="."), 6, 5)
    print({
        p <- ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
            stat_boxplot(geom ='errorbar')+
            geom_boxplot(colour = "black", outlier.colour=NA, notch=TRUE) +
            theme_classic() +
            scale_fill_manual(values=cols) +
            ylab("log2(PausingIndex FC)")+
            xlab("")+
            ggtitle(basename(outName))+
            theme(panel.grid.major = element_blank()
                 ,panel.grid.minor = element_blank(), 
                  panel.background = element_blank()
                 ,axis.line = element_line(colour = "black"))+
            theme(text = element_text(size=12),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size=12),
                  plot.title=element_text(size=11))
    })
    dev.off()    
}else{
    print("Calculating ratio of body over promoter")

    Ratio <- do.call(cbind, lapply(Order, function(x){
        pseudo       <- min( min(Tss[Tss!=0]) , min(Body[Body!=0]) )
        pi           <-  (Body.or[, paste0(x, "_body")]+pseudo) / (Tss.or[, paste0(x, "_tss")]+pseudo)
        pi           <- as.data.frame(pi)
        rownames(pi) <- rownames(Tss.or)
        colnames(pi) <- x
        pi
    }
    ))
    fc <- do.call(cbind, lapply(1:length(exps), function(i){
        print(paste("FC sample", i, exps[i]))
        FC           <- Ratio[,exps[i]] / Ratio[,con]
        FC           <- as.data.frame(FC)
        rownames(FC) <- rownames(Ratio)
        colnames(FC) <- sub("$", "_fc",exps[i])
        FC
    }))
    for (i in 1:length(exps)) {
        for (j in 1:length(exps)){
            q1 <- log2(fc[,i]) > log2(FC) & log2(fc[,j]) > log2(FC)
            q2 <- log2(fc[,j]) > log2(FC) & log2(fc[,i]) < -log2(FC)
            q3 <- log2(fc[,i]) < -log2(FC) & log2(fc[,j]) < -log2(FC)
            q4 <- log2(fc[,i]) > log2(FC) & log2(fc[,j]) < -log2(FC)
            png(paste(outName, names(fc)[i], "vs", names(fc)[j], "PromterReleaseFC.png", sep="."), 500, 500)
            par(mar=c(5,5,5,5))
            plot(log2(fc[,i]), log2(fc[,j]),
                 xlab=paste("log2(",names(fc[i]),")",sep=""),
                 ylab=paste("log2(",names(fc[j]),")",sep=""),
                 main="Promoter Release Ratio (Body/Promoter) \nFold Change",
                 cex.axis=1.1,
                 cex.lab=1.1,
                 col="#C0C0C070",
                 pch=19)
            if(plotPcc==1){
                legend("topleft",legend= round(cor(fc[,i],fc[,j],method="pearson"), 4), cex=1.1)
            }
            points(log2(fc[ q1, i]), log2(fc[q1, j]), col=q1Col, pch=19,cex=1)
            points(log2(fc[ q3, i]), log2(fc[q3, j]), col=q3Col, pch=19,cex=1)
            points(log2(fc[ q2, i]), log2(fc[q2, j]), col=q2Col, pch=19,cex=1)
            points(log2(fc[ q4, i]), log2(fc[q4, j]), col=q4Col, pch=19,cex=1)
            legend("bottomright", title=paste("FC:", FC)
                  ,legend=c(paste("q1",sum(q1))
                           ,paste("q2",sum(q2))
                           ,paste("q3",sum(q3))
                           ,paste("q4",sum(q4))
                           )
                 ,pch=c(19,19,19,19), col=c(q1Col,q2Col,q3Col,q4Col),cex=0.9,)
            abline(h=0, lty="dashed", col="black")
            abline(v=0, lty="dashed", col="black")
            dev.off()
        }
    }
    Ratio$gene_id <- rownames(Ratio)
    fc$gene_id    <- rownames(fc)
    write.table(fc, file=paste(outName, "PromterReleaseFC.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(Ratio, file=paste(outName, "PromterRelease.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)    
    ## violin plot
    df.long <- melt(fc[,grep("gene_id", names(fc), invert=TRUE)])
    df.long$variable <- sub("_fc$", "", df.long$variable)
    pdf(paste(outName, "violinPlotPromterReleaseFC.pdf", sep="."), 6, 5)
    print({
        p <-
            ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
            geom_violin() +
            geom_boxplot(width=0.1, colour = "black", outlier.colour=NA, notch=TRUE) +
            theme_classic() +
            scale_fill_manual(values=cols) +
            ylab("log2(PRR FC)")+
            xlab("")+
            ggtitle(basename(outName))+
            theme(panel.grid.major = element_blank()
                 ,panel.grid.minor = element_blank(), 
                  panel.background = element_blank()
                 ,axis.line = element_line(colour = "black"))+
            theme(text = element_text(size=12),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size=11),
                  plot.title=element_text(size=11))
    })
    dev.off()
    ## boxPlot
    pdf(paste(outName, "boxPlotPromterReleaseFC.pdf", sep="."), 6, 5)
    print({
        p <- ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
            stat_boxplot(geom ='errorbar')+
            geom_boxplot(colour = "black", outlier.colour=NA, notch=TRUE) +
            theme_classic() +
            scale_fill_manual(values=cols) +
            ylab("log2(PRR FC)")+
            xlab("")+
            ggtitle(basename(outName))+
            theme(panel.grid.major = element_blank()
                 ,panel.grid.minor = element_blank(), 
                  panel.background = element_blank()
                 ,axis.line = element_line(colour = "black"))+
            theme(text = element_text(size=12),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size=12),
                  plot.title=element_text(size=11))
    })
    dev.off()
}

print("done")
