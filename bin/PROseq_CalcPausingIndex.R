args <- commandArgs()

help <- function(){
    cat("PROseq_CalcPausingIndex.R :
- From a provided GRnages object of transcripts calcualte the average coverage in a defined promoter region and body region
\tfor in a strand aware manner .
- You pass an option to use a defined distance rather than the entire gene body.
- You can use either the annotated TSS or the tssMaxStart from PROseq_pickGenesTSS.R \n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest      [required]\n")    
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                      [default = hg19]\n")    
    cat("--promUp      : number of nucleotides upstream of the tss                  [required]\n")
    cat("--promDown    : number of nucleotides downstream the tss                   [required]\n")
    cat("--bodyLen     : if you want a smaller set body region (ex 2000)            [default = entire gene length w/o promoter]\n")
    cat("--bwFiles     : path to bigWig Files                                       [required]\n")
    cat("--bwPattern   : grep pattern for bigWigs to use quotes (ex. PolII.*293T)   [optional]
                         if not provided uses all bw in path\n")
    cat("--numCores    : number of cores to use                                     [default = 10 ]\n")
    cat("--outName     : prefix to your out file names                              [default = basename(bigWigFile) ]\n")
    cat("--tssMax      : use the max postion of coverage in the annotated tss (0/1) [default = 0; use annotated ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions   <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    promUp    <- as.numeric( sub( '--promUp=', '', args[grep('--promUp=', args)] ))
    promDown  <- as.numeric( sub( '--promDown=', '', args[grep('--promDown=', args)] ))
    bodyLen  <- sub( '--bodyLen=', '', args[grep('--bodyLen=', args)] )
    bwFiles   <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])
    tssMax    <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])
}

outName
bwPattern
bwFiles

#regions <- "tables/PRO_0h_NELFeAID_DLD_1092.filteredProteinCodingTx.rda"
#setwd("/projects/b1025/arw/analysis/yuki/degrons/DLD1/")
#promUp <- 25
#promDown <- 25
#numCores <- 10
#bodyLen=2000
#tssMax <- 1
#outName <- "tables/pausing_index/PRO_0h_NELFeAID_DLD_1092_Tss25bp"
#bwFiles <- "data_proseq"
#bwPattern <- "PRO.*NELFe.*1092"


## test files
#regions <- "tables/PRO_0h_Aux_PAF1AID_DLD1_846.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#bwFiles <- "data_proseq"
#bwPattern <- "Aux_PAF1AID_DLD1_846"
#promUp <- 5
#promDown <- 5
#numCores <- 10
#bodyLen <- 2000
#tssMax <- 1
#outName <- "tables/pausing_index/PROseq_Aux_PAF1AID_DLD1_846_Tss10bp"
#Cores <- 10

## select all bigWig files
bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
bws <- bws[grep(bwPattern,bws,invert=FALSE)]
bws

if (identical(tssMax,character(0))){
   tssMax <- 0
}else{
    tssMax <- as.numeric(tssMax)
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
   outName <- sub(".filteredProteinCodingTx.rda|.bed", "", basename(regions))
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
    ## this is for enseble version 75
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    anno <- get(load(annoFile))
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

###############################
## set up tss and body regions
###############################

## read in transcripts
if( length(grep(".bed$", regions)) > 0 ){
    boo <-  read.table(regions)
    names(boo) <- c("seqnames", "start", "end", "gene_id", "tx_name", "strand")
    gnModel <- as(boo, "GRanges")
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
}else{
    gnModel          <- get(load(regions))
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
}

## define TSS
if ( tssMax== 1 ){
    print("use max as Tss")
    maxTss                                                      <- gnModel
    print("fix start")
    start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart
    print("fix end")
    end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart
    print("define tss")
    Tss                                                         <- promoters(maxTss ,upstream=promUp ,downstream=promDown)
    if (identical(bodyLen,character(0))){
        Body      <- resize(maxTss, fix='end', width=width(maxTss)-promDown )
    }else{
        print(paste("using a body length of:", bodyLen))
        bodyLen   <- as.numeric(bodyLen)
        maxTss.rs <- resize(maxTss,fix='start',width= bodyLen + promDown)
        Body      <- resize(maxTss.rs,fix='end',width=width(maxTss.rs) - promDown )
    }
}else{
    Tss  <- promoters(gnModel ,upstream=promUp ,downstream=promDown)
    ranges(Tss)
    if (identical(bodyLen,character(0))){
        Body    <- resize(gnModel, fix='end', width=width(gnModel)-promDown )
    }else{
        print(paste("using a body length of:", bodyLen))
        bodyLen <- as.numeric(bodyLen)
        gnModel <- resize(gnModel,fix='start',width= bodyLen + promDown)
        Body    <- resize(gnModel,fix='end',width=width(gnModel) - promDown )
    }
    ranges(Body)
}
ranges(Tss)
ranges(Body)

names(Tss) <- NULL
names(Body) <- NULL

### get rid of zombie processes
library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

#for (iter in 1:T) {
#  res[iter] = mclapply(1:10, fun)
#  wait()
#}


###############################
## calc average coverage in tss and body regions
###############################
calcCov <- function(bw,model){
    print(basename(bw))
    print("import plus strand")
    plus      <- import.bw(paste0(bw, ".plus.bw"),  RangedData=FALSE,selection = BigWigSelection(model))
    print("import minus strand")
    minus     <- import.bw(paste0(bw, ".minus.bw"), RangedData=FALSE,selection = BigWigSelection(model))
    print("calc cov for plus strand")
    plus.cov  <- coverage(plus, weight='score')
    print("calc cov for minus strand")
    minus.cov <- coverage(minus, weight='score') *-1
    print("calc average cov for genes")
    mean.cov  <- with(as.data.frame(model),{
        mcmapply(function(Seqname,start,end,strand){
            if(strand == '-'){
                signif(mean(minus.cov[[Seqname]][start:end]))
            }else{
                signif(mean(plus.cov[[Seqname]][start:end]))
            }
        }   
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })
    print("convert to data frame")
    mean.cov           <- data.frame(mean.cov)
    print("add gene name")
    rownames(mean.cov) <- model$gene_id
    print("add column name")
    colnames(mean.cov) <- sub(".bw", "", basename(bw))
    mean.cov
}

bwBase            <- unique(sub(".minus.bw|.plus.bw", "", bws))
print("calc tssCov")
tssCov            <- do.call(cbind, lapply(bwBase, model=Tss, calcCov))
#tssCov            <- do.call(cbind, mclapply(bwBase, model=Tss, calcCov, mc.cores=1, mc.preschedule=FALSE))
print("calc bodyCov")
bodyCov           <- do.call(cbind, lapply(bwBase, model=Body, calcCov))
#bodyCov           <- do.call(cbind, mclapply(bwBase, model=Body, calcCov, mc.cores=1, mc.preschedule=FALSE))

print("fix the column names")
colnames(tssCov)  <- sub("$", "_tss",  colnames(tssCov))
colnames(bodyCov) <- sub("$", "_body", colnames(bodyCov))

stopifnot(rownames(tssCov)==rownames(bodyCov))

## combine tss and body dataframes
print("combine data frames")
df <- cbind(as.data.frame(gnModel), tssCov, bodyCov)

print("check the directories")
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

## save the gnModel as a granges object and a tab delimited file
print("write table")
write.table(df
           ,file=paste0(outName, ".pausingIndexAverageCoverages.txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )


print("done")
