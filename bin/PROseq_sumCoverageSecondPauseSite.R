args <- commandArgs()

help <- function(){
    cat("PROseq_sumCoverageSecondPauseSite.R :
- From a provided GRnages sum the coverage downstream of the max for a set number of bases while excluding the max position.
- GRanges object must have a tssMaxStart column to reposition ranges to.
- We are excluding the maximum position, so it will not be included in the sum
- You will be returned a table with the ranges, sums and fold changes and a boxplot.\n")
 
    cat("Usage: \n")
    cat("--regions    : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--downStream : distance downstream of the tssMax site to sum the coverage     [default = 50]\n")
    cat("--assembly   : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")
    cat("--bwFiles    : path to bigWig Files                                           [required]\n")
    cat("--bwPattern  : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [default = all .bw in path]\n")
    cat("--Control    : bw prefix of the control sample used to calculate fold changes [default = all .bw in path]
                          relative to (No .minus.bw or .plus.bw extension)                                         \n")
    cat("--numCores   : number of cores to use                                         [default = 10 ]\n")
    cat("--outName    : prefix to your out file names (No .extension)                  [required]\n")
    cat("--Cols       : table of colors for box plot (need a sample and color column)  [default = rainbow(n)]\n")
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args)) ){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFiles    <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern  <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Control    <- sub( '--Control=', '', args[grep('--Control=', args)])
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
}

bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
if (! identical(bwPattern,character(0))){
    bws <- bws[grep(bwPattern,bws,invert=FALSE)]
}
bws

if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(downStream,character(0))){
    print("downStream not provided use 50bp")
    downStream <- 50
}else{
    downStream <- as.numeric(downStream)
    print(paste("downStream provided as", downStream))
}

cat("regions:", regions, sep="\n")
cat("outName:", outName, sep="\n")
cat("bwPattern:", bwPattern, sep="\n")
cat("bwFiles:", bwFiles, sep="\n")
cat("Control:", Control, sep="\n")
cat("downStream:", downStream, sep="\n")
cat("assembly:", assembly, sep="\n")

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens"
    species <- "Homo sapiens"}
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus"
    species <- "Mus musculus"}
if (assembly == "sacCer3") { organismStr <- "Scerevisiae"
    species <- "Saccharomyces cerevisiae"}
if (assembly == "dm3") { organismStr <- "Dmelanogaster"
    species <- "Drosophila melanogaster"}
if (assembly == "rn6") { organismStr <- "Rnorvegicus"
    species <- "Rattus norvegicus"}

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)
library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(gtools)

if (assembly == "hg19") {
    organism <- Hsapiens
    }
if (assembly == "mm9") {
    organism <- Mmusculus
}
if (assembly == "mm10") {
    organism <- Mmusculus
}
if (assembly == "sacCer3") {
    organism <- Scerevisiae
}
if (assembly == "dm3") {
    organism <- Dmelanogaster
}

## read in transcripts
if( length(grep(".bed$", regions)) > 0 ){
    Model          <- import.bed(regions)
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}else{
    Model          <- get(load(regions))
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}

## define TSS

print("use max as Tss")
maxTss                                                      <- Model
print("fix start and remove first position")
start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart + 1 
print("fix end and remove first position")
end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart - 1
print("define tss")
Model.win                                                   <- promoters(maxTss ,upstream=0 ,downstream=downStream)
fname                                                       <- sub("$", paste0("_sumCovMaxPlus1_down", downStream), outName)
ranges(Model.win)

## make sure the output directory exists 
Dir <- dirname(outName)
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

###############################
## calculate average coverage
###############################
calcCov <- function(bw,model){
    print(basename(bw))
    plus      <- import.bw(paste0(bw, ".plus.bw"),  RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(paste0(bw, ".minus.bw"), RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    print("sum coverage")
    winCov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                signif(sum(minus.cov[[seqname]][start:end]))
            }else{
                signif(sum(plus.cov[[seqname]][start:end]))
            }
        }   
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })
    winCov           <- data.frame(winCov)
    rownames(winCov) <- model$gene_id
    colnames(winCov) <- sub(".bw", "", basename(bw))
    winCov
}

bwBase            <- unique(sub(".minus.bw|.plus.bw", "", bws))
avgCov            <- do.call(cbind,mclapply(mixedsort(bwBase), model=Model.win, calcCov, mc.cores=1))

#############################
## calculate the fold change
#############################
SAMPLES <- colnames(avgCov)
SAMPLES

exps <- SAMPLES[SAMPLES!=Control]
exps

pseudo <- min(avgCov[!avgCov==0])

fc <- do.call(cbind, lapply(1:length(exps), function(i){
    print(paste("FC sample", i, exps[i]))
    FC               <- log2( (avgCov[,exps[i]] + pseudo) / (avgCov[,Control] + pseudo) )
    df.return        <- data.frame(FC)
    names(df.return) <- exps[i]
    df.return
}))

names(fc) <- sub("$", "_log2FC", names(fc))
names(avgCov) <- sub("$", "_sum", names(avgCov))

## combine resized regions and average coverage
df <- cbind(as.data.frame(Model.win), signif(avgCov, 4), signif(fc, 4))

print("make directory for output file if it does not exist.")
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}


write.table(df
           ,file=paste0(fname, ".txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )

###############################
## make boxplot
###############################
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(ggsignif)

print("make fold change boxplot")

df.long          <- melt(fc)
df.long$variable <- sub("_log2FC", "", df.long$variable)

df.long$variable <-factor(df.long$variable, levels=exps, ordered = TRUE)

if (identical(cols,character(0))){
    Cols <- rainbow(ncol(fc))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}

#minMax <- do.call(rbind, lapply(1:ncol(avgCov), function(x){
#     iQr  <- IQR(avgCov[,x])
#     Ymax <- round(as.numeric(quantile(avgCov[, x])[4] + 1.5 * iQr) *1.05)
#     Ymin <- floor(as.numeric(quantile(avgCov[, x])[2] - 1.5 * iQr) *1.05)
#     cbind(Ymin, Ymax)
#}
#))

minMax <- apply(fc, 2, function(x)boxplot.stats(x)$stats[c(1, 5)])

mround <- function(x,base){ 
    base*round(x/base) 
} 

Ymin <- mround(min(minMax)-0.25, 0.5)
Ymax <- mround(max(minMax)+0.25, 0.5)



pdf(file=sub("$", ".boxPlot.pdf", fname),width=6,height=5)
print({
    p <-
    ggplot(df.long, aes(x=variable, y=value, fill=variable)) + 
        stat_boxplot(geom ='errorbar')+
        geom_boxplot(notch=TRUE,outlier.colour=NA) +
        theme_classic() +
        scale_fill_manual("sample", values=Cols)+
        #coord_cartesian(ylim = c(Ymin, Ymax)) +
        ylab( "log2( Coverage Sum )" ) + 
    xlab("") +
    scale_y_continuous(limits=c(Ymin, Ymax), breaks=scales::pretty_breaks(n = Ymax-Ymin))+#, breaks=seq(Ymin,Ymax, by=0.5))+
        ggtitle( paste( downStream, "downstream Pause Site" ) ) +
        theme(text = element_text(size=8),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=8))
})
dev.off()

print("done")

