args <- commandArgs()

help <- function(){
    cat("PROseq_tssCoverageBoxplots.R :
- From a provided GRnages object of transcripts calcualte the average or total coverage in a defined promoter region
\t either from the defined TSS or the max coverage of PROseq in the tss region (column named tssMaxStart).\n")
 
    cat("Usage: \n")
    cat("--regions     : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--type        : Peaks, Tss, or Tes                                             [required]\n")
    cat("--window      : total size of region around types above (use for type = Peaks) [default = 0]\n")
    cat("--upStream    : distance upstream of the region to take (use for Tss and Tes)  [default = 50]\n")
    cat("--downStream  : distance downstream of the region to take (use for Tss and Tes)[default = 50]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")
    cat("--bwFiles     : path to bigWig Files                                           [required]\n")
    cat("--bwPattern   : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [default = all .bw in path]\n")
    cat("--tssMax      : use the max postion of coverage in the annotated tss (0/1)     [default = 0; use annotated ]
                           column must be named tssMaxStart\n")
    cat("--numCores    : number of cores to use                                         [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                  [required]\n")
    cat("--cols        : need the same number as samples separated by comma             [default = rainbow]\n")
    cat("--calcSum     : calculate the coverage sum rather than the average (0/1)       [default = 0]\n")    
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args)) ){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    Type       <- sub( '--type=', '', args[grep('--type=', args)] )
    Window     <- sub( '--window=', '', args[grep('--window=', args)] )
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFiles    <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern  <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols       <- sub( '--cols=', '',args[grep('--cols=',args)])
    tssMax     <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])
    calcSum    <- sub( '--calcSum=', '',args[grep('--calcSum=',args)])
}

outName
bwPattern
bwFiles

bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
if (! identical(bwPattern,character(0))){
    bws <- bws[grep(bwPattern,bws,invert=FALSE)]
}
bws

## test files
#setwd("/projects/b1025/arw/analysis/yuki/degrons/DLD1")
#regions <- "tables/PRO_0h_NELFcAID_DLD_1057.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#bwFiles <- "data_proseq"
#bwPattern <- "NELFcAID_DLD_1057"
#upStream <- 25
#downStream <- 25
#numCores <- 10
#tssMax <- 1
#outName <- "plots/pausing_index/PRO_0h_NELFcAID_DLD_1057_Tss25bp"
#Cores <- 10
#calcSum <- 1

if (identical(tssMax,character(0))){
   tssMax <- 0
}else{
    tssMax <- as.numeric(tssMax)
}

if (identical(calcSum,character(0))){
   calcSum <- 0
}else{
    calcSum <- as.numeric(calcSum)
}

if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}


if (identical(upStream,character(0))){
    print("upStream not provided use 50bp")
    upStream <- 50
}else{
    upStream <- as.numeric(upStream)
    print(paste("upStream provided as", upStream))
}
if (identical(downStream,character(0))){
    print("downStream not provided use 50bp")
    downStream <- 50
}else{
    downStream <- as.numeric(downStream)
    print(paste("downStream provided as", downStream))
}

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
if ( tssMax== 1 ){
    print("use max as Tss")
    maxTss                                                      <- Model
    print("fix start")
    start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart
    print("fix end")
    end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart
    print("define tss")
    Model.win                                                   <- promoters(maxTss ,upstream=upStream ,downstream=downStream)
    fname <- sub("$", paste0("_MaxTss", "up", upStream, "down", downStream), outName)
}else{
    print("use annotated Tss")
    Model.win   <- promoters(Model ,upstream=upStream ,downstream=downStream)
    fname <- sub("$", paste0("_Tss", "up", upStream, "down", downStream), outName)
}
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
    if (calcSum == 0 ){
        print("average coverage")
        winCov  <- with(as.data.frame(model),{
            mcmapply(function(seqname,start,end,strand){
                if(strand == '-'){
                    signif(mean(minus.cov[[seqname]][start:end]))
                }else{
                    signif(mean(plus.cov[[seqname]][start:end]))
                }
            }   
           ,mc.cores=Cores
           ,as.character(seqnames),start,end,as.character(strand))
        })
    }else{
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
    }
    winCov           <- data.frame(winCov)
    rownames(winCov) <- model$gene_id
    colnames(winCov) <- sub(".bw", "", basename(bw))
    winCov
}

bwBase            <- unique(sub(".minus.bw|.plus.bw", "", bws))
avgCov            <- do.call(cbind,mclapply(bwBase, model=Model.win, calcCov, mc.cores=1))

#colnames(avgCov) <- sub("$", paste0("_", Type), colnames(avgCov))

## combine resized regions and average coverage
df <- cbind(as.data.frame(Model.win), signif(avgCov, 4))

## save the gnModel as a granges object and a tab delimited file

if (calcSum == 0 ){
    fname <- paste0(fname, "_avgCov")
}else{
    fname <- paste0(fname, "_sumCov")
}
fname
    
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

df.long <- melt(avgCov)

if (identical(cols,character(0))){
    Cols <- rainbow(ncol(avgCov))
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

minMax <- apply(avgCov, 2, function(x)boxplot.stats(x)$stats[c(1, 5)])

Ymin <- floor(log2(min(minMax[1,]) +0.001))
Ymax <- ceiling(log2(max(minMax[2,]) +0.001))

if(calcSum > 0){
    yLab <- "log2( Coverage Sum )"
}else{
    yLab <- "log2( Coverage Average )"
}

if(tssMax > 0){
    plotTitle <- paste(upStream, "upstream and", downStream, "downstream Pause Site")
}else{
    plotTitle <- paste(upStream, "upstream and", downStream, "downstream Annotated TSS")
}

pdf(file=sub("$", ".boxPlot.pdf", fname),width=6,height=5)
print({
    p <-
    ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
        stat_boxplot(geom ='errorbar')+
        geom_boxplot(notch=TRUE,outlier.colour=NA) +
        theme_classic() +
        scale_fill_manual("sample", values=Cols)+
        #coord_cartesian(ylim = c(Ymin, Ymax)) +
        ylab( yLab ) + 
        xlab("") +
        ggtitle(plotTitle) +
        theme(text = element_text(size=8),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=8))
})
dev.off()









