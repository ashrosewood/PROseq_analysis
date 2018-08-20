args <- commandArgs()

help <- function(){
    cat("PROseq_metaGeneMatrix.R :
- From a provided GRnages obnject of genes of interest make metaGene matrix of defined windows of upstream of the TSS
  and down stream the Tes.
- This script is for strand separated bigWigs. Only the coverage on the same strand as the gene is used (sense strand).
- By default it approximates the body length to 400 positions and up and downstream to 100 (table has 600 columns).
- If Bins are provided, it will multiply these approx lengths by these bins, so the resulting table is still 600 columns. \n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest           [required]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                           [default = hg19 ]\n")    
    cat("--plusBw      : bigWigFile for plus strand                                      [required]\n")
    cat("--minusBw     : bigWigFile for minus strand                                     [required]\n")    
    cat("--upStream    : how many bases upstream of the Tss                              [default = 2000 ]\n")
    cat("--downStream  : how many bases downstream of the Tes                            [default = 2000 ]\n")
    cat("--numCores    : number of cores to use                                          [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                   [default = basename(bigWigFile) ]\n")
    cat("--tssMax      : use the max postion of coverage in the annotated tss (0/1)      [default = 0; use annotated ]\n")
    cat("--Bins        : Number of bins to take the average coverage after approximating [default = 1; no binning ]
                            genes to the same length.                                                                 \n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    plusBw     <- sub( '--plusBw=', '', args[grep('--plusBw=', args)])
    minusBw    <- sub( '--minusBw=', '', args[grep('--minusBw=', args)])
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    tssMax     <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])
    Bins       <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])
}

if (identical(upStream,character(0))){
   upStream <- 2000
}else{
   upStream <- as.numeric(upStream)
}

if (identical(downStream,character(0))){
   downStream <- 2000
}else{
   downStream <- as.numeric(downStream)
}

if (identical(Cores,character(0))){
   Cores <- 6
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bw", "", basename(bwFile))
}

if (identical(tssMax,character(0))){
    tssMax <- 0
}else{
    tssMax <- as.numeric(tssMax)
}

if (identical(Bins,character(0))){
    Bins      <- 1
    ApproxLen <- 400
}else{
    Bins      <- as.numeric(Bins)
    ApproxLen <- 400 * Bins
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

## take annotated or maxTss sites to center
if ( tssMax== 1 ){
    print("use max as Tss")
    maxTss                                                      <- get(load(regions))
    seqinfo(maxTss)                                             <- seqinfo(organism)[seqlevels(maxTss)]
    print("fix start")
    start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart
    print("fix end")
    end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart
    maxTss$length                                               <- width(maxTss)
    Model                                                       <- maxTss
    seqlevels(Model,force=TRUE)                                 <- seqlevels(Model)[grep("chrM",seqlevels(Model), invert=TRUE)]
}else{
    print("use annotated Tss")
    Model                       <- get(load(regions))
    seqinfo(Model)              <- seqinfo(organism)[seqlevels(Model)]
    Model$length                <- width(Model)
    seqlevels(Model,force=TRUE) <- seqlevels(Model)[grep("chrM",seqlevels(Model), invert=TRUE)]
}

TSS <- promoters(Model, upstream = upStream, downstream = 0)
TES <- resize(Model, fix = 'end', width = 1)
TES <- promoters(TES,upstream=0,downstream=downStream)
Body <- Model

head(ranges(TSS))
head(ranges(Body))
head(ranges(TES))

fname <- sub("$", ".metaGene.rda", outName)

## only want to use the sense strand to make metagenes with for stranded data
Bin <- function(model,Size){
    ## read in coverage for regions of interest
    plus                  <- import.bw(plusBw,RangedData=FALSE,selection = BigWigSelection(model))
    minus                 <- import.bw(minusBw,RangedData=FALSE,selection = BigWigSelection(model))
    ## take the score column as coverage
    plus.cov              <- coverage(plus, weight='score')
    minus.cov             <- coverage(minus, weight='score') *-1
    ## be sure coverage has seqlengths
    seqlengths(plus.cov)  <- seqlengths(organism)[seqlevels(plus.cov)]
    seqlengths(minus.cov) <- seqlengths(organism)[seqlevels(minus.cov)]
    ## approximate the coverage and flip if minus
    Cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                r <- minus.cov[[seqname]][start:end]
                r <- rev(r)
            }else{
                r <- plus.cov[[seqname]][start:end]
            }
            r <- approx(r,n=Size)$y
            return(r)
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })    
    cov           <- as.data.frame(t(Cov))
    if( Bins > 1){
        window.cov <- function(row){
                window <- as.integer(ncol(cov)/Bins)
                window.coverage <- lapply(0:(window-1), function(jump)
                    rowMeans(row[as.numeric(jump*Bins+1):as.numeric((jump*Bins)+Bins)])
                    )
                t(as.matrix(unlist(window.coverage)))
        }
        win <- mclapply(1:nrow(cov), function(i)
            window.cov(cov[i,]),mc.cores=Cores)    
        bin.mat      <- do.call(rbind, mclapply(win, as.numeric, mc.cores=Cores))        
        df           <- data.frame(bin.mat)
        rownames(df) <- model$gene_id
    }else{
        df           <- cov
        rownames(df) <- model$gene_id        
    }
    df
}

tssMat  <- Bin(model=TSS,  Size=ApproxLen/4)    
bodyMat <- Bin(model=Body, Size=ApproxLen)    
tesMat  <- Bin(model=TES,  Size=ApproxLen/4)    

stopifnot(rownames(tssMat)==rownames(tesMat))
stopifnot(rownames(tssMat)==rownames(bodyMat))

df <- data.frame(cbind(tssMat, bodyMat, tesMat))
rownames(df) <- rownames(tssMat)
save(df,file=fname)
