args <- commandArgs()

help <- function(){
    cat("PROseq_defineMaxCovForGeneList.R :
- For a provided GRanges object or bed file, define the pause site (max coverage in the TSS) for a bigWig file.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest              [required]\n")    
    cat("--assembly     : genome (hg19, mm9, mm10, dm3, sacCer3)                            [default = hg19]\n")
    cat("--txTssDown    : number of nucleotides downstream tss for tx picking and filtering [default = 500]\n")    
    cat("--plusBw       : bigWigFile for plus strand                                        [required]\n")
    cat("--minusBw      : bigWigFile for minus strand                                       [required]\n")
    cat("--numCores     : number of cores to use                                            [default = 10 ]\n")
    cat("--outName      : prefix to your out file names                                     [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions      <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly     <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    txTssDown    <- sub( '--txTssDown=', '', args[grep('--txTssDown=', args)])
    plusBw       <- sub( '--plusBw=', '', args[grep('--plusBw=', args)])
    minusBw      <- sub( '--minusBw=', '', args[grep('--minusBw=', args)])
    Cores        <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName      <- sub( '--outName=', '',args[grep('--outName=',args)])
}

## for testing purposes
#regions         <- "/projects/b1025/arw/analysis/yuki/degrons/DLD1/tables/PRO_0h_Parental_DLD_1093.filteredProteinCodingTx.rda"
#plusBw           <- "/projects/b1025/arw/analysis/yuki/degrons/DLD1/data_proseq/PRO_0h_NELFcAID_DLD_1057.plus.bw"
#minusBw          <- "/projects/b1025/arw/analysis/yuki/degrons/DLD1/data_proseq/PRO_0h_NELFcAID_DLD_1057.minus.bw"
#peakFile         <- "/projects/b1025/tango/TANGO-817/TANGO-817/peaks/293-DMSO-PolIIrep1.macsPeaks.bed" 
#assembly         <- "hg19"
#txTssDown        <- 500
#Cores            <- 10
#outName          <- "/projects/b1025/arw/analysis/yuki/degrons/DLD1/tables/PRO_0h_NELFcAID_DLD_1057.parentMax"

if (identical(txTssDown,character(0))){
   txTssDown <- 500
}else{
    txTssDown <- as.numeric(txTssDown)
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
   outName <- sub(".plus.bw", "", basename(plusBw))
}

cat("regions:", regions, sep="\n")
cat("plusBw:", plusBw, sep="\n")
cat("minusBw:", minusBw, sep="\n")
cat("txTssDown:", txTssDown, sep="\n")
cat("numCores:", Cores, sep="\n")
cat("outName:", outName, sep="\n")
cat("assembly:", assembly, sep="\n")

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
library(GenomicFeatures)
library(parallel)
library(biomaRt)

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
    print("regions are in bed format")
    print("columns are assumed as seqnames, start, end, gene_id, tx_name, strand")
    boo <-  read.table(regions)
    names(boo) <- c("seqnames", "start", "end", "gene_id", "tx_name", "strand")
    gnModel <- as(boo, "GRanges")
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
}else{
    print("regions are in GRanges object")
    print("throws out any columns with the word max or Max in it")
    gnModel          <- get(load(regions))
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
    mcols(gnModel)   <- mcols(gnModel)[,names(mcols(gnModel))[grep("Max|max", names(mcols(gnModel)), invert=TRUE)] ]
}


###########################################
## Find the max in the promoter region
###########################################
TSS            <- promoters(gnModel, upstream=1, downstream=txTssDown)
TSS$gene_id    <- names(TSS)
names(TSS)     <- NULL
head(ranges(TSS))


getMax <- function(model){
    plus      <- import.bw(plusBw,RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(minusBw,RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1
    max.cov  <- data.frame(t(
        with(as.data.frame(model),{
            mcmapply(function(seqname,start,end,strand,tx_name){
                if(strand == '-'){
                    bwMax       <- max(minus.cov[[seqname]][start:end])
                    maxPos      <- which.max(minus.cov[[seqname]][start:end])
                    inputMax    <- minus.cov[[seqname]][start+maxPos-1]@values
                    NewStart    <- as.integer(start+maxPos-1)
                }else{
                    bwMax       <- max(plus.cov[[seqname]][start:end])
                    maxPos      <- which.max(plus.cov[[seqname]][start:end])
                    inputMax    <- plus.cov[[seqname]][start+maxPos-1]@values
                    NewStart    <- start+maxPos-1
                }
                c(tx_name, bwMax, NewStart)
            }   
           ,mc.cores=Cores
           ,as.character(seqnames),start,end,as.character(strand),as.character(tx_name))
        })
    ))
    names(max.cov)      <- c("tx_name", "tssMaxCpm", "tssMaxStart")
    max.cov$tssMaxCpm   <- as.numeric(paste(max.cov$tssMaxCpm))
    max.cov$tssMaxStart <- as.numeric(paste(max.cov$tssMaxStart))
    max.cov$tx_name     <- as.character(paste(max.cov$tx_name))
    rownames(max.cov)   <- max.cov$tx_name
    max.cov
}

max.df <- getMax(model=TSS)

stopifnot( max.df[gnModel$tx_name,"tx_name"] == paste(gnModel$tx_name))
#[1] TRUE

elementMetadata(gnModel)[["tssMaxCpm"]]   <- max.df[gnModel$tx_name,"tssMaxCpm"]
elementMetadata(gnModel)[["tssMaxStart"]] <- max.df[gnModel$tx_name,"tssMaxStart"]

### to filter can not use the max cpm but need an average coverage in the first kb or so

## save the gnModel as a granges object and a tab delimited file
save(gnModel, file=paste0(outName, ".rda"))
write.table(as.data.frame(gnModel)
           ,file=paste0(outName, ".txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )
