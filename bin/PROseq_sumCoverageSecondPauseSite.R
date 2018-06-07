args <- commandArgs()

help <- function(){
    cat("PROseq_sumCoverageSecondPauseSite.R :
- From a provided GRnages sum the coverage downstream of the max for a set number of bases while excluding a set number of positions
\t after the max position. The total range used will downStream - tssMaxDown. tssMaxDown must be < downStream.
- GRanges object must have a tssMaxStart column to reposition ranges to.
- You can exclude maximum position and a number of bases after the tssMax, so it will not be included in the sum
- You will be returned a table with the ranges, sums and fold changes and a boxplot.
- Both tssMaxDown and downStream will be appended to the end of outName.\n")
 
    cat("Usage: \n")
    cat("--regions    : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--tssMaxDown : distance downstream of the tssMax to start coverage sum        [default = 10]
                           use 0 if you wish to include the tssMax                                    \n")
    cat("--downStream : distance downstream of the tssMax to end coverage sum          [default = 100]\n")
    cat("--assembly   : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")
    cat("--bwFiles    : path to bigWig Files                                           [required]\n")
    cat("--bwPattern  : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [default = all .bw in path]\n")
    cat("--Control    : bw prefix of the control sample used to calculate fold changes [default = all .bw in path]
                          relative to (No .minus.bw or .plus.bw extension)                                         \n")
    cat("--numCores   : number of cores to use                                         [default = 10 ]\n")
    cat("--outName    : prefix to your out file names (No .extension)                  [required]\n")
    cat("--cols       : table of colors for box plot (need a sample and color column)  [default = rainbow(n)]\n")
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args)) ){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    tssMaxDown <- sub( '--tssMaxDown=', '', args[grep('--tssMaxDown=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFiles    <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern  <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Control    <- sub( '--Control=', '', args[grep('--Control=', args)])
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols       <- sub( '--cols=', '',args[grep('--cols=',args)])
}

if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(tssMaxDown,character(0))){
    print("tssMaxDown not provided use 50bp")
    tssMaxDown <- 10
}else{
    print(paste("tssMaxDown provided as", tssMaxDown))
    tssMaxDown <- as.numeric(tssMaxDown)
}

if (identical(downStream,character(0))){
    print("downStream not provided use 50bp")
    downStream <- 100
}else{
    downStream <- as.numeric(downStream)
    print(paste("downStream provided as", downStream))
}

cat("regions:", regions, sep="\n")
cat("tssMaxDown:", tssMaxDown, sep="\n")
cat("downStream:", downStream, sep="\n")
cat("outName:", outName, sep="\n")
cat("bwPattern:", bwPattern, sep="\n")
cat("bwFiles:", bwFiles, sep="\n")
cat("Control:", Control, sep="\n")
cat("assembly:", assembly, sep="\n")

bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
if (! identical(bwPattern,character(0))){
    bws <- bws[grep(bwPattern,bws,invert=FALSE)]
}
bws

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
print(paste("fix start to tssMaxStart"))
start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart 
print("fix end and remove first position")
end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart 

##
print("take window 0,100 at tssMax")
Model.win                                                   <- promoters(maxTss ,upstream=0 ,downstream=downStream)
print(ranges(Model.win))

##
print(paste("trim off the first", tssMaxDown, "bases"))
Model.win                                                   <-  resize(Model.win, fix='end',width=width(Model.win) - tssMaxDown )
print(ranges(Model.win))

##
fname                                                       <- sub("$", paste0("_sumCovMaxPlus",tssMaxDown, "_down", downStream), outName)
ranges(Model.win)

## make sure the output directory exists
print("make directory for output file if it does not exist.")
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

write.table(df
           ,file=paste0(fname, ".txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )

###############################
## make boxplots
###############################
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(ggsignif)

## log2FC
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

minMax <- apply(fc, 2, function(x)boxplot.stats(x)$stats[c(1, 5)])

mround <- function(x,base){ 
    base*round(x/base) 
} 

Ymin <- mround(min(minMax)-0.25, 0.5)
Ymax <- mround(max(minMax)+0.25, 0.5)

pdf(file=sub("$", ".boxPlotFC.pdf", fname),width=6,height=5)
print({
    p <-
    ggplot(df.long, aes(x=variable, y=value, fill=variable)) + 
        stat_boxplot(geom ='errorbar')+
        geom_boxplot(notch=TRUE,outlier.colour=NA) +
        theme_classic() +
        scale_fill_manual("sample", values=Cols)+
        #coord_cartesian(ylim = c(Ymin, Ymax)) +
        ylab( "log2( Coverage Sum FC )" ) + 
    xlab("") +
    scale_y_continuous(limits=c(Ymin, Ymax), breaks=scales::pretty_breaks(n = Ymax-Ymin))+#, breaks=seq(Ymin,Ymax, by=0.5))+
        ggtitle( paste( tssMaxDown, "downstream Pause Site to plus", downStream ) ) +
        theme(text = element_text(size=8),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=8))
})
dev.off()

## avgCov
print("make average coverage boxplot")

df.long          <- melt(avgCov + pseudo)
df.long$variable <- sub("_sum", "", df.long$variable)

df.long$variable <-factor(df.long$variable, levels=SAMPLES, ordered = TRUE)

minMax <- apply(log10(avgCov+pseudo), 2, function(x)boxplot.stats(x)$stats[c(1, 5)])

Ymin <- mround(min(minMax)-0.25, 0.5)
Ymax <- mround(max(minMax)+0.25, 0.5)

pdf(file=sub("$", ".boxPlotSum.pdf", fname),width=6,height=5)
print({
    p <-
    ggplot(df.long, aes(x=variable, y=log10(value), fill=variable)) + 
    stat_boxplot(geom ='errorbar')+
    geom_boxplot(notch=TRUE,outlier.colour=NA) +
    theme_classic() +
    scale_fill_manual("sample", values=c("white", Cols))+
    #coord_cartesian(ylim = c(Ymin, Ymax)) +
    ylab( "log10( Coverage Sum )" ) + 
    xlab("") +
    scale_y_continuous(limits=c(Ymin, Ymax), breaks=scales::pretty_breaks(n = Ymax-Ymin))+#, breaks=seq(Ymin,Ymax, by=0.5))+
        ggtitle( paste( tssMaxDown, "downstream Pause Site to plus", downStream ) ) +
        theme(text = element_text(size=8),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=8))
})
dev.off()

print("done")

