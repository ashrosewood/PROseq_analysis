args <- commandArgs()

help <- function(){
    cat("PROseq_TssFcPlots.R :
- Compare average Tss coverage fold changes for 2 conditions relative to a control.
- Gives a scatter plot for all possible comparisons of sample fold changes, violin plot, box plot and a table with the fold changes. \n")
    cat("Usage: \n")
    cat("--Table   : path to table of average tss coverages (columns must have _tss)          [required]\n")
    cat("--Control : Name of Control sample to calculate FC over                               [optional]\n")
    cat("--FC      : Fold change to highlight genes in each quadrant                          [default=2]\n")
    cat("--outName : Prefix for plot name and table. Each comparison appended to the end.     [required]
                     (ex. /home/ashley/Pol2_293T_rep1)                                                  \n")
    cat("--q1Col   : Quadrant 1 color (defaults are RdBu colors )                             [optional]\n")
    cat("--q2Col   : Quadrant 2 color                                                         [optional]\n")
    cat("--q3Col   : Quadrant 3 color                                                         [optional]\n")
    cat("--q4Col   : Quadrant 4 color                                                         [optional]\n")
    cat("--plotPcc : Flag to add Pearson correlation (0/1)                                    [default=1]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args))){
    help()
} else {
    Table      <- sub( '--Table=', '', args[grep('--Table=', args)])
    Control    <- sub( '--Control=', '', args[grep('--Control=', args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    FoldChange <- sub( '--FC=', '',args[grep('--FC=',args)])
    q1Col      <- sub( '--q1Col=', '',args[grep('--q1Col=',args)])
    q2Col      <- sub( '--q2Col=', '',args[grep('--q2Col=',args)])
    q3Col      <- sub( '--q3Col=', '',args[grep('--q3Col=',args)])
    q4Col      <- sub( '--q4Col=', '',args[grep('--q4Col=',args)])
    plotPcc    <- sub( '--plotPcc=', '',args[grep('--plotPcc=',args)])
}

if (identical(plotPcc,character(0))){
    plotPcc <- 1
}else{
    plotPcc <- as.numeric(plotPcc)
}

if (identical(FoldChange,character(0))){
    FoldChange      <- 2
}else{
    FoldChange      <- as.numeric(FoldChange)
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

cat("Table:" ,    Table ,      sep="\n")
cat("outName:" ,  outName ,    sep="\n")
cat("Control:" ,  Control ,    sep="\n")
cat("FC:" ,       FoldChange , sep="\n")

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

df            <- read.table(Table, sep="\t", header=TRUE)
rownames(df)  <- paste(df$gene_id)
dim(df)

Tss           <- df[ , grep("tss$", names(df)) ]
print(head(Tss))
colnames(Tss) <- sub("_tss$", "", colnames(Tss))

## get sample names
Order        <- sub( "_tss", "", names(Tss) )

## put the control first and set experiments
exps <- Order[grep(Control, Order, invert=TRUE)]
exps

con  <- Order[grep(Control, Order, invert=FALSE)]
con

if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

fc <- do.call(cbind, lapply(1:length(exps), function(i){
    print(paste("FC sample", i, exps[i]))
    pseudo       <- min(Tss[Tss!=0])
    print(paste("pseudo count:", pseudo))
    FC           <- (Tss[,exps[i]]+pseudo) / (Tss[,con]+pseudo)
    FC           <- as.data.frame(FC)
    rownames(FC) <- rownames(Tss)
    colnames(FC) <- sub("$", "_fc",exps[i])
    FC
}))

for (i in 1:length(exps)) {
    for (j in 1:length(exps)){        
        q1 <- log2(fc[,i]) > log2(FoldChange) & log2(fc[,j]) > log2(FoldChange)
        q2 <- log2(fc[,j]) > log2(FoldChange) & log2(fc[,i]) < -log2(FoldChange)
        q3 <- log2(fc[,i]) < -log2(FoldChange) & log2(fc[,j]) < -log2(FoldChange)
        q4 <- log2(fc[,i]) > log2(FoldChange) & log2(fc[,j]) < -log2(FoldChange)
        png(paste(outName, names(fc)[i], "vs", names(fc)[j], "TssCovFC.png", sep="."), 500, 500)
        plot(log2(fc[,i]), log2(fc[,j]),
             xlab=paste("log2(",names(fc[i]),")",sep=""),
             ylab=paste("log2(",names(fc[j]),")",sep=""),
             main="Tss average coverage\nFold Change",
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
        legend("bottomright", title=paste("FC:", FoldChange)
              ,legend=c(paste("q1",sum(q1))
                       ,paste("q2",sum(q2))
                       ,paste("q3",sum(q3))
                       ,paste("q4",sum(q4))
                        )
              ,pch=c(19,19,19,19), col=c(q1Col,q2Col,q3Col,q4Col),cex=0.9,)
        abline(h=0, lty="dashed", col="black")
        abline(v=0, lty="dashed", col="black")
        dev.off()
        png(paste(outName, names(fc)[i], "vs", names(fc)[j], "densityTssCovFC.png", sep="."), 500, 500)
        smoothScatter(log2(fc[,i]), log2(fc[,j]),
                     ,colramp=colorRampPalette(c("white", brewer.pal(9, "Blues")))
                     ,col=NA
                     ,xlab=paste("log2(",names(fc[i]),")",sep="")
                     ,ylab=paste("log2(",names(fc[j]),")",sep="")
                     ,main="Tss average coverage\nFold Change"
                     ,nbin=300
                     ,cex=2
                      )
        abline(h=0, lty=2, col="black")
        abline(v=0, lty=2, col="black")
        points(log2(fc[ q1, i]), log2(fc[q1, j]), col=q1Col, pch=19,cex=1)
        points(log2(fc[ q3, i]), log2(fc[q3, j]), col=q3Col, pch=19,cex=1)
        points(log2(fc[ q2, i]), log2(fc[q2, j]), col=q2Col, pch=19,cex=1)
        points(log2(fc[ q4, i]), log2(fc[q4, j]), col=q4Col, pch=19,cex=1)
        legend("bottomright", title=paste("FC:", FoldChange)
              ,legend=c(paste("q1",sum(q1))
                       ,paste("q2",sum(q2))
                       ,paste("q3",sum(q3))
                       ,paste("q4",sum(q4))
                        )
              ,pch=c(19,19,19,19), col=c(q1Col,q2Col,q3Col,q4Col),cex=0.9,)
        dev.off()
    }
}

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))(500)[seq(1,500,length.out=length(exps))]

df.long <- melt(fc)
df.long$variable <- sub("_fc$", "", df.long$variable)

pdf(paste(outName, "violinPlotTssCovFC.pdf", sep="."), 6, 5)
print({
    p <- ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
        geom_violin() +
        geom_boxplot(width=0.1, colour = "black", outlier.colour=NA, notch=TRUE) +
        theme_classic() +
        scale_fill_manual(values=cols) +
        ylab("log2(Tss avg FC)")+
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


#minMax <- log2(apply(fc, 2, function(x)boxplot.stats(x)$stats))#[c(1, 5)]))

#stats <- boxplot.stats(fc)$stats
#df <- data.frame(x="label1", ymin=stats[1], lower=stats[2], middle=stats[3], 
#                 upper=stats[4], ymax=stats[5])


pdf(paste(outName, "boxPlotTssCovFC.pdf", sep="."), 6, 5)
print({
    p <- ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
    stat_boxplot(geom ='errorbar')+
    geom_boxplot(colour = "black", outlier.colour=NA, notch=TRUE) +
        theme_classic() +
        scale_fill_manual(values=cols) +
        ylab("log2(Tss avg FC)")+
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




fc$gene_id    <- rownames(fc)
write.table(fc, file=paste(outName, "TssCovFC.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("done")
