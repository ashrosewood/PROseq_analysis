args <- commandArgs()

help <- function(){
    cat("PROseq_findSecondPauseSite.R :
- Some recent studies indicate that there is evidence of a second pausing site slightly downstream.
- From a provided GRnages find the second pause site (max coverage) after the tssMaxsite.
\t after the max position. The total range used will downStream - tssMaxDown. tssMaxDown must be < downStream.
- Table returned containing the resized regions, new max positions and the absolute difference (delta).\n")
 
    cat("Usage: \n")
    cat("--regions    : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--tssMaxDown : distance downstream of the tssMax to start looking for the max [default = 10]\n")
    cat("--downStream : distance downstream of the tssMax to stop looking for the max  [default = 100]\n")
    cat("--assembly   : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")
    cat("--bwFiles    : path to bigWig Files                                           [required]\n")
    cat("--bwPattern  : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [default = all .bw in path]\n")
    cat("--numCores   : number of cores to use                                         [default = 10 ]\n")
    cat("--outName    : prefix to your out file names (No .extension)                  [required]\n")
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
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
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
#cat("Control:", Control, sep="\n")
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
print("take window 0 to downStream of tssMax")
Model.win                                                   <- promoters(maxTss ,upstream=0 ,downstream=downStream)
print(ranges(Model.win))

##
print(paste("trim off the first", tssMaxDown, "bases"))
Model.win                                                   <-  resize(Model.win, fix='end',width=width(Model.win) - tssMaxDown )
print(ranges(Model.win))

##
fname                                                       <- sub("$", paste0("_secondPauseMaxPlus",tssMaxDown, "_down", downStream), outName)

## make sure the output directory exists
print("make directory for output file if it does not exist.")
Dir <- dirname(outName)
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

###############################
## calculate average coverage
###############################
getMax <- function(bw, model){
    print(basename(bw))
    plus      <- import.bw(paste0(bw, ".plus.bw"),  RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(paste0(bw, ".minus.bw"), RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    print("find the second max in coverage")
    max.cov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand,tx_name){
            if(strand == '-'){
                bwMax       <- max(minus.cov[[seqname]][start:end])
                maxPos      <- which.max(minus.cov[[seqname]][start:end])
                inputMax    <- minus.cov[[seqname]][start+maxPos-1]@values
                NewStart    <- as.integer(start+maxPos-1)
                NewStart
            }else{
                bwMax       <- max(plus.cov[[seqname]][start:end])
                maxPos      <- which.max(plus.cov[[seqname]][start:end])
                inputMax    <- plus.cov[[seqname]][start+maxPos-1]@values
                NewStart    <- start+maxPos-1
                NewStart
            }            
        }   
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand),as.character(tx_name))
    })
    secondPause           <- as.data.frame(max.cov, row.names=model$gene_id)
    colnames(secondPause) <- sub("$", "_max2", basename(bw))
    secondPause
}

bwBase            <- unique(sub(".minus.bw|.plus.bw", "", bws))
max2              <- do.call(cbind,mclapply(mixedsort(bwBase), model=Model.win, getMax, mc.cores=1))
print(head(max2))

df <- cbind(as.data.frame(Model.win), max2)

#############################
## calculate delta
#############################
SAMPLES <- basename(bwBase)
SAMPLES
print(SAMPLES)

delta <- do.call(cbind, lapply(1:length(SAMPLES), function(i){
    print(paste("calculate sample 2nd max delta", i, SAMPLES[i]))
    Num              <- abs(df[,paste0(SAMPLES[i], "_max2")] - df$tssMaxStart)
    df.return        <- data.frame(Num, row.names=rownames(df))
    names(df.return) <- SAMPLES[i]
    df.return
}))
names(delta) <- sub("$", "_max2delta", names(delta))

## combine resized regions and average coverage
df.all <- cbind(df, delta)

write.table(df.all
           ,file=paste0(fname, ".txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )

print("done")

