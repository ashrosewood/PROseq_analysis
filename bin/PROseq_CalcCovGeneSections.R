args <- commandArgs()

help <- function(){
    cat("PROseq_CalcCovGeneSections.R :
-For a GRanges object of transcripts of interest section the gene body after removing the tss and defined downstream region
\tand calculate the transcript total coverage and total coverage with each of the specified number of sections.
-There is an option for PROseq to use the start site of the max coverage (tssMaxStart)from PROseq_pickGenesTssFromStart.R
\tto define the promoter region to be removed.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest          [required]\n")    
    cat("--assembly    : genome (hg19, mm9, mm10, dm3, sacCer3)                         [hg19]\n")
    cat("--promDown    : number of nucleotides downstream the maxTss position to remove [required]\n")
    cat("--numSections : number of sections to break the gene body into                 [default = 4 ]\n")    
    cat("--bwFiles     : path to bigWig Files                                           [required]\n")
    cat("--bwPattern   : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [optional]
                         if not provided uses all bw in path\n")
    cat("--outName     : prefix to your out file names                                  [required ]\n")
    cat("--numCores    : number of cores to use                                         [default = 6 ]\n")
    cat("--tssMax      : use the max postion of coverage in the annotated tss (0/1)     [default = 0; use annotated ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions     <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly    <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    promDown    <- as.numeric( sub( '--promDown=', '', args[grep('--promDown=', args)] ))
    numSections <- sub( '--numSections=', '', args[grep('--numSections=', args)])
    bwFiles     <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern   <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    outName     <- sub( '--outName=', '',args[grep('--outName=',args)])
    Cores       <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    tssMax    <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])

}

## for testing 
#regions     <- "tables/PRO_0h_Aux_PAF1AID_DLD1_846tss500.filteredProteinCodingTx.rda"
#assembly    <- "hg19"
#bwFiles     <- "data_proseq"
#bwPattern   <- "Aux_PAF1AID_DLD1_846"
#promDown    <- 500
#Cores    <- 10
#outName <- "tables/elongation_index/PROseq_Aux_PAF1AID_DLD1_846tss500_TssDown500sections4SumCov.txt"
#tssMax <- 1
#numSections <- 4

if (identical(Cores,character(0))){
   Cores <- 6
}else{
    Cores <- as.numeric(Cores)
}

## set defaults
if (identical(numSections,character(0))){
   numSections <- 4
}else{
    numSections <- as.numeric(numSections)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
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

## select all bigWig files
bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
bws <- bws[grep(bwPattern,bws,invert=FALSE)]
bws

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
    start(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']) <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='+']$tssMaxStart
    end(maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-'])   <- maxTss[paste(as.data.frame(maxTss)[,"strand"])=='-']$tssMaxStart
    Tss                                                         <- promoters(maxTss ,upstream=0 ,downstream=promDown)
    Body                                                        <- resize(maxTss, fix='end', width=width(maxTss)-promDown )
    Model <- maxTss
}else{
    Tss                                                         <- promoters(gnModel ,upstream=0 ,downstream=promDown)
    Body                                                        <- resize(gnModel, fix='end', width=width(gnModel)-promDown )
    Model <- gnModel
}
##
print("tss ranges")
ranges(Tss)
print("body ranges")
ranges(Body)

####################
## set up function
####################
Body$gene_id    <- names(Body)
names(Body)     <- NULL

sumCov <- function(bw,model){
    basename(bw)
    plus      <- import.bw(paste0(bw, ".plus.bw"),  RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(paste0(bw, ".minus.bw"), RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    sum.cov  <- with(as.data.frame(model),{
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
    sum.cov           <- data.frame(sum.cov)
    rownames(sum.cov) <- model$gene_id
    colnames(sum.cov) <- sub(".bw", "", basename(bw))
    sum.cov
}

bwBase            <- unique(sub(".minus.bw|.plus.bw", "", bws))
## entire tx coverage
entireTx <- do.call(cbind,mclapply(bwBase,model=Model,sumCov,mc.cores=1))
names(entireTx) <- sub("$", "_entireTx", names(entireTx))

########################
## section the gene
########################

for (i in 1:numSections){
    ## Set up GB sections, positive strand
    print(paste("Section",i,"positive"))
    gbSection                                        <- Body
    gbSection$lengths                                <- width(gbSection)
    ## plus strand
    start(gbSection[which(strand(gbSection) =="+")]) <- start(Body[which(strand(Body) =="+")]) + (i-1)*(round(gbSection[which(strand(Body) =="+")]$lengths/numSections))
    end(gbSection[which(strand(gbSection)=="+")])    <- start(Body[which(strand(Body)=="+")]) + i*round(gbSection[which(strand(Body) =="+")]$lengths/numSections)
    ## negative strand; in granges for minus strand the start is the end
    print(paste("Section",i,"negative"))
    start(gbSection[which(strand(Body) =="-")])   <- end(Body[which(strand(Body) =="-")]) - i*(gbSection[which(strand(Body) =="-")]$lengths/numSections)
    end(  gbSection[which(strand(Body)=="-") ])   <- end(Body[which(strand(Body) =="-")]) - (i-1)*(gbSection[which(strand(Body) =="-")]$lengths/numSections)
    ## calc the total cov in that region
    sectionSum                                       <- do.call(cbind,mclapply(bwBase,model=gbSection,sumCov,mc.cores=1))
    names(sectionSum)                                <- sub("$", paste0("_cut",i), names(sectionSum))
    assign(paste0("cut_", i), sectionSum) 
}

all <- cbind(entireTx, get("cut_1"), get("cut_2"), get("cut_3"), get("cut_4")) 

df <- as.data.frame(gnModel)

stopifnot(gnModel$gene_id == rownames(all))

df.all <- cbind(df, all)

write.table(df.all
           ,file=outName
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
