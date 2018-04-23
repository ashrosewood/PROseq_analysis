args <- commandArgs()

help <- function(){
    cat("PROseq_metaGeneMatrix.R :
- From a provided GRnages obnject of genes of interest make metaGene matrix of defined windows of upstream of the TSS
  and down stream the Tes.
- This script is for strand separated bigWigs. Only the coverage on the same strand as the gene is used (sense strand).\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest         [required]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                         [default = hg19 ]\n")    
    cat("--plusBw      : bigWigFile for plus strand                                    [required]\n")
    cat("--minusBw     : bigWigFile for minus strand                                   [required]\n")    
    cat("--upStream    : how many bases upstream of the Tss                            [default = 2000 ]\n")
    cat("--downStream  : how many bases downstream of the Tes                          [default = 2000 ]\n")
    cat("--numCores    : number of cores to use                                        [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                 [default = basename(bigWigFile) ]\n")
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
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])

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

Model                       <- get(load(regions))
seqinfo(Model)              <- seqinfo(organism)[seqlevels(Model)]
seqlevels(Model,force=TRUE) <- seqlevels(Model)[grep("chrM",seqlevels(Model), invert=TRUE)]

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
    Cov           <- as.data.frame(t(Cov))
    rownames(Cov) <- model$gene_id
    Cov
}

tssMat  <- Bin(model=TSS,  Size=100)    
bodyMat <- Bin(model=Body, Size=400)    
tesMat  <- Bin(model=TES,  Size=100)    

stopifnot(rownames(tssMat)==rownames(tesMat))
stopifnot(rownames(tssMat)==rownames(bodyMat))

df <- data.frame(cbind(tssMat, bodyMat, tesMat))
rownames(df) <- rownames(tssMat)
save(df,file=fname)