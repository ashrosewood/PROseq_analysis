args <- commandArgs()

help <- function(){
    cat("PROseq_pickGenesTSS.R :
- This script generates a filtered gene list from a Txdb based on PRO-seq in defined region around tss and filters for ensembl protein coding
  genes that are also refseq validated.
- Filter for a minimal gene length.
- Filter for a minimal distance to the nearest gene after selecting the best transcript.
- If a peakFile is provided, then take the transcripts only overlapping peaks.
- Require a max rpm in the tss region of at least 1.
- This script is set up for PolII ChIP-seq and does not take the stranded coverage into account\n")
    cat("Usage: \n")
    cat("--txdbFile    : transcript data base (.txdb)                             [required]\n")
    cat("--assembly    : genome (hg19, mm9, mm10, dm3, sacCer3)                   [default = hg19]\n")
    cat("--txTssDown   : number of nucleotides downstream tss for tx picking only [default = 50]\n")    
    cat("--minLength   : min transcript length                                    [default = 2000]\n")
    cat("--minDist     : min distance between genes                               [default = 2000]\n")
    cat("--plusBw      : bigWigFile for plus strand                               [required]\n")
    cat("--minusBw     : bigWigFile for minus strand                              [required]\n")
    cat("--peakFile    : bed file                                                 [optional]\n")
    cat("--numCores    : number of cores to use                                   [default = 10 ]\n")
    cat("--outName     : prefix to your out file names                            [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    txdbFile    <- sub( '--txdbFile=', '', args[grep('--txdbFile=', args)] )
    assembly    <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    txTssDown   <- sub( '--txTssDown=', '', args[grep('--txTssDown=', args)])
    minLength   <- sub( '--minLength=', '', args[grep('--minLength=', args)])
    minDist     <- sub( '--minDist=', '', args[grep('--minDist=', args)])
    plusBw      <- sub( '--plusBw=', '', args[grep('--plusBw=', args)])
    minusBw     <- sub( '--minusBw=', '', args[grep('--minusBw=', args)])
    peakFile    <- sub( '--peakFile=', '', args[grep('--peakFile=', args)])
    Cores       <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName     <- sub( '--outName=', '',args[grep('--outName=',args)])
}

#txdbFile         <- "/projects/b1025/arw/anno/hg19/hg19_ucsc_refGene.txdb"
#txdbFile         <- "/projects/p20742/anno/Txdb/hsapiens_gene_ensembl_Ens75.txdb"
#plusBw           <- "/projects/b1025/arw/analysis/fei/paf1/data_proseq/PRO_0h_Aux_PAF1AID_DLD1_846.plus.bw"
#minusBw          <- "/projects/b1025/arw/analysis/fei/paf1/data_proseq/PRO_0h_Aux_PAF1AID_DLD1_846.minus.bw"
#peakFile         <- "/projects/b1025/tango/TANGO-817/TANGO-817/peaks/293-DMSO-PolIIrep1.macsPeaks.bed" 
#assembly         <- "hg19"
#txTssDown        <- 300
#minLength        <- 3000
#minDist          <- 1000
#Cores            <- 10
#outName          <- "tables/PRO_0h_Aux_PAF1AID_DLD1_846"

if (identical(txTssDown,character(0))){
   txTssDown <- 100
}else{
    txTssDown <- as.numeric(txTssDown)
}

if (identical(minLength,character(0))){
   minLength <- 2000
}else{
    minLength <- as.numeric(minLength)
}

if (identical(minDist,character(0))){
   minDist <- 2000
}else{
    minDist <- as.numeric(minDist)
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
library(GenomicFeatures)
library(parallel)
library(biomaRt)

if (assembly == "hg19") {
    organism <- Hsapiens
    ## this is for enseble version 75
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    if(!(file.exists(annoFile))) {  
        bm <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
        anno <- getBM(mart=bm, attributes=c('ensembl_gene_id'
                                           ,'ensembl_transcript_id'
                                           ,'external_gene_name'
                                           ,'gene_biotype'
                                           ,'transcript_biotype'
                                           ,'refseq_mrna'
                                           ,'refseq_ncrna'
                                           ,'entrezgene'
                                           ,'description'))
        anno[anno$refseq_mrna=="",  "refseq_mrna"] <- NA
        anno[anno$refseq_ncrna=="", "refseq_ncrna"] <- NA
        save(anno, file=annoFile)
    }else{
        anno <- get(load(annoFile))
    }}
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

##############
## tss average coverage
##############

## load model
txdb <- loadDb(txdbFile)

## filter out chroms
seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$|^M$",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

## toss chrY and chrM
seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("chrY|chrM",seqlevels(txdb), invert=TRUE)]

if(length(grep("^chr",seqlevels(txdb)))==0){
    seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
}

## make gnModel to combine with counts for only looking at the browser
## rpkms if calculated with need to use the exons below
gnModel          <- transcriptsBy(txdb, 'gene')
gnModel          <- unlist(gnModel)## Gets the genomic region convered by transcripts
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]

fiveUTRs          <- fiveUTRsByTranscript(txdb, use.names=TRUE)
fiveUTRs          <- unlist(fiveUTRs)
seqinfo(fiveUTRs) <- seqinfo(organism)[seqlevels(fiveUTRs)]
fiveUTRs$tx_name  <- names(fiveUTRs)

## add gene info to gnModel
gnModel$gene_id            <- names(gnModel)
#iv                         <- match(gnModel$tx_name, anno$refseq_mrna)
iv                         <- match(gnModel$tx_name, anno$ensembl_transcript_id)
gnModel$external_gene_name <- anno[iv, "external_gene_name"]
gnModel$gene_biotype       <- anno[iv, "gene_biotype"]
gnModel$tx_biotype         <- anno[iv, "transcript_biotype"]
gnModel$refseq_mrna        <- anno[iv, "refseq_mrna"]

## filter for only protein coding with refseq mrna id
gnModel <- gnModel[!is.na(gnModel$refseq_mrna) & gnModel$tx_biotype=="protein_coding" & gnModel$gene_biotype=="protein_coding"]

fiveUTRs <- fiveUTRs[fiveUTRs$tx_name %in% gnModel$tx_name]
fiveUTRs <- fiveUTRs[fiveUTRs$exon_rank==1]

length(gnModel)
length(fiveUTRs)
## not the same

###########################################
## pick the mosth highly occupied transcript tss
########################################### 
TSS            <- promoters(gnModel, upstream=0, downstream=txTssDown)
TSS$gene_id    <- names(TSS)
names(TSS)     <- NULL

calcCov <- function(model){
    plus      <- import.bw(plusBw,RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(minusBw,RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    mean.cov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                mean(minus.cov[[seqname]][start:end])
            }else{
                mean(plus.cov[[seqname]][start:end])
            }
        }   
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })       
    mean.cov           <- data.frame(mean.cov)
    rownames(mean.cov) <- model$tx_name
    mean.cov$mean.cov
}

elementMetadata(gnModel)[["tssAvgCov"]]      <- calcCov(model=TSS)
elementMetadata(fiveUTRs)[["fiveUtrAvgCov"]] <- calcCov(model = resize(fiveUTRs, width=width(fiveUTRs)+150, fix='start') )

gnModel$fiveUtrAvgCov                         <- NA
names(fiveUTRs)                               <- paste(fiveUTRs$tx_name)

boo                                           <- as.data.frame(fiveUTRs)
iv                                            <- match(gnModel$tx_name, boo$tx_name)
gnModel$fiveUtrAvgCov                         <- boo[iv, "fiveUtrAvgCov"]
gnModel$AvgCov                                <- gnModel$fiveUtrAvgCov

## sort and take most highly occupied tss
or         <- order( elementMetadata(gnModel)$AvgCov,decreasing=TRUE )
gnModel.or <- gnModel[or]
reps       <- duplicated( gnModel.or$gene_id )
gnModel    <- gnModel.or[!reps]

## calculate distance to the nearest gene that does not directly overlap
dists                     <- as.data.frame( distanceToNearest(gnModel, ignore.strand=TRUE) )
gnModel$distanceToNearest <- dists$distance

dists                     <- as.data.frame( distanceToNearest(gnModel, ignore.strand=FALSE) )
gnModel$distanceToNearestStranded <- dists$distance

###########################################
## Filter for tss's with peaks if provided
###########################################

if ( ! identical( peakFile, character(0)) ){
    peaks                                               <- import.bed(peakFile)
    ol                                                  <- as.data.frame(findOverlaps(promoters(gnModel, upstream=0, downstream=txTssDown)
                                                                                     ,peaks)
                                                                         )
    gnModel$peakOverLap                                 <- 0
    mcols(gnModel[unique(ol$queryHits)])["peakOverLap"] <- 1
    gnModel                                             <- gnModel[gnModel$peakOverLap == 1 ]
    gnModel$peakOverLap                                 <- NULL
}

###########################################
## Filter for minimal length
###########################################
gnModel <- gnModel[width(gnModel) >= minLength]

###########################################
## Filter for minimal distance between genes
###########################################
gnModel <- gnModel[gnModel$distanceToNearestStranded >= minDist]

###########################################
## Find the max in the promoter region
###########################################

UTR <- fiveUTRs[fiveUTRs$tx_name %in% gnModel$tx_name]
mcols(UTR) <- mcols(UTR)[,c("tx_name","fiveUtrAvgCov"  )]
names(UTR) <- NULL
UTR$group <- "UTR"

noUTR <- gnModel[is.na(gnModel$fiveUtrAvgCov)]
mcols(noUTR) <- mcols(noUTR)[,c("tx_name","fiveUtrAvgCov"  )]
names(noUTR) <- NULL 
noUTR$group <- "noUTR"

summary(width(resize(UTR, width=width(UTR)+150, fix='start')))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  151.0   222.0   287.0   342.4   396.0  3870.0

regions <- c(resize(UTR, width=width(UTR)+150, fix='start')
            ,promoters(noUTR, upstream=0, downstream=200)
             )

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

max.df <- getMax(model=regions)

all.equal(max.df[gnModel$tx_name,"tx_name"], paste(gnModel$tx_name))
#[1] TRUE

elementMetadata(gnModel)[["tssMaxCpm"]]   <- max.df[gnModel$tx_name,"tssMaxCpm"]
elementMetadata(gnModel)[["tssMaxStart"]] <- max.df[gnModel$tx_name,"tssMaxStart"]

### to filter can not use the max cpm but need an average coverage in the first kb or so
maxWin <- gnModel

start(maxWin[paste(as.data.frame(maxWin)[,"strand"])=='+']) <- maxWin[paste(as.data.frame(maxWin)[,"strand"])=='+']$tssMaxStart

end(maxWin[paste(as.data.frame(maxWin)[,"strand"])=='-']) <- maxWin[paste(as.data.frame(maxWin)[,"strand"])=='-']$tssMaxStart

calcCov <- function(model, Fun){
    plus      <- import.bw(plusBw,RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(minusBw,RangedData=FALSE,selection = BigWigSelection(model))
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    mean.cov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                if(Fun=="mean"){
                    mean(minus.cov[[seqname]][start:end])
                }else{
                    sum(minus.cov[[seqname]][start:end])
                }
            }else{
                if(Fun=="mean"){
                    mean(plus.cov[[seqname]][start:end])
                }else{
                    sum(plus.cov[[seqname]][start:end])
                }
            }
        }   
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })       
    mean.cov           <- data.frame(mean.cov)
    rownames(mean.cov) <- model$tx_name
    mean.cov$mean.cov
}

## try sum instead of mean
elementMetadata(maxWin)[["up50down1000SumCov"]]      <- calcCov(model=promoters(maxWin, upstream=50, downstream=1000), Fun="sum")

stopifnot(names(maxWin) == names(gnModel))

#test <- maxWin[maxWin$up50down1000SumCov > 250 & maxWin$tssMaxCpm >= 9]
#test[order( elementMetadata(test)$win20AvgCov,decreasing=TRUE )]
#length(test)

gnModel$maxUp50down1000SumCov <- maxWin$up50down1000SumCov

summary(gnModel$maxUp50down1000SumCov)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.0     0.0   111.0   184.8   274.6  5247.0
summary(gnModel$tssMaxCpm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   3.661   8.229  10.980 392.900 
summary(gnModel[gnModel$maxUp50down1000SumCov > 250]$tssMaxCpm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.220   9.762  14.640  20.510  23.180 392.900 

gnModel <- gnModel[gnModel$maxUp50down1000SumCov > 250 & maxWin$tssMaxCpm >= 9 ]

gnModel$tx_id <- NULL
gnModel$tx_biotype <- NULL
gnModel$gene_biotype <- NULL

## save the gnModel as a granges object and a tab delimited file
if ( ! identical( outName, character(0)) ){
    save(gnModel, file=paste0(outName, ".filteredProteinCodingTx.rda"))
    write.table(as.data.frame(gnModel)
               ,file=paste0(outName, ".filteredProteinCodingTx.txt")
               ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
                )
}else{
    save(gnModel, file=paste0(sub(".bw", "", basename(bwFile)), ".filteredProteinCodingTx.rda"))
    write.table(as.data.frame(gnModel)
               ,file=paste0(sub(".bw", "", basename(bwFile)), ".filteredProteinCodingTx.txt")
               ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
                )
}
