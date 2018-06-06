args <- commandArgs()

help <- function(){
    cat("PROseq_TssFC_vs_BodyFC_scatter.R :
- Compare body fold changes vs tss fold changes and return a scatter plot and table with 
- with all the average coverages and fold changes. \n")
    cat("Usage: \n")
    cat("--Table   : path to table of averge body and tss coverages                           [required]
                     columns sample names followed by _pi or _prr, last column is gene_id.            \n")
    cat("--Control : Name of Control sample to calculate FC over                              [required]\n")
    cat("--FC      : Fold change to highlight genes in each quadrant                          [default=2]\n")

    cat("--outName : Prefix for plot name and table. Each comparison appended to the end.     [required]
                     (ex. /home/ashley/Pol2_rep1)                                                 \n")
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
    FC      <- sub( '--FC=', '',args[grep('--FC=',args)])
    q1Col   <- sub( '--q1Col=', '',args[grep('--q1Col=',args)])
    q2Col   <- sub( '--q2Col=', '',args[grep('--q2Col=',args)])
    q3Col   <- sub( '--q3Col=', '',args[grep('--q3Col=',args)])
    q4Col   <- sub( '--q4Col=', '',args[grep('--q4Col=',args)])
    plotPcc <- sub( '--plotPcc=', '',args[grep('--plotPcc=',args)])
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

print("Calculating ratios for body and promoter")

BodyFC <- do.call(cbind, lapply(exps, function(x){
    pseudo       <- min(Body[Body!=0]) 
    pi           <-  (Body.or[, paste0(x, "_body")]+pseudo) / (Body.or[, paste0(con, "_body")]+pseudo)
    pi           <- as.data.frame(pi)
    rownames(pi) <- rownames(Tss.or)
    colnames(pi) <- x
    pi
}
))
names(BodyFC) <- sub("$", "_bodyFC", names(BodyFC))

TssFC <- do.call(cbind, lapply(exps, function(x){
    pseudo       <- min(Tss[Tss!=0]) 
    pi           <-  (Tss.or[, paste0(x, "_tss")]+pseudo) / (Tss.or[, paste0(con, "_tss")]+pseudo)
    pi           <- as.data.frame(pi)
    rownames(pi) <- rownames(Tss.or)
    colnames(pi) <- x
    pi
}
))
names(TssFC) <- sub("$", "_tssFC", names(TssFC))

for (i in 1:length(exps)) {
    q1 <- log2(BodyFC[,i]) > log2(FC) & log2(TssFC[,i]) > log2(FC)
    q2 <- log2(BodyFC[,i]) < -log2(FC) & log2(TssFC[,i]) > log2(FC)
    q3 <- log2(BodyFC[,i]) < -log2(FC) & log2(TssFC[,i]) < -log2(FC)
    q4 <- log2(BodyFC[,i]) > log2(FC) & log2(TssFC[,i]) < -log2(FC)
    png(paste(outName, names(BodyFC)[i], "vs", names(TssFC)[i], "scatter.png", sep="."), 500, 500)
    par(mar=c(5,5,5,5))
    plot(x=log2(BodyFC[,i]), y=log2(TssFC[,i]),
         xlab=paste("log2(",names(BodyFC[i]),")",sep=""),
         ylab=paste("log2(",names(TssFC[i]),")",sep=""),
         main="Body vs Promoter \nFold Change",
         cex.axis=1.1,
         cex.lab=1.1,
         col="#C0C0C070",
         pch=19)
    if(plotPcc==1){
        legend("topleft",legend= round(cor(TssFC[,i],BodyFC[,i],method="pearson"), 4), cex=1.1)
    }
    points(log2(BodyFC[ q1, i]), log2(TssFC[q1, i]), col=q1Col, pch=19,cex=1)
    points(log2(BodyFC[ q3, i]), log2(TssFC[q3, i]), col=q3Col, pch=19,cex=1)
    points(log2(BodyFC[ q2, i]), log2(TssFC[q2, i]), col=q2Col, pch=19,cex=1)
    points(log2(BodyFC[ q4, i]), log2(TssFC[q4, i]), col=q4Col, pch=19,cex=1)
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

stopifnot(rownames(BodyFC) == rownames(TssFC))
out.df <- cbind(df, TssFC, BodyFC) 

write.table(out.df, file=paste(outName, "FoldChange.txt", sep="."), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


print("done")
