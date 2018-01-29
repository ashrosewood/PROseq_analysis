args <- commandArgs()

help <- function(){
    cat("PROseq_heatMatrix.R :
- From a provided GRnages obnject make binned heatmap matrix of defined windows around peaks, tss, or tes.
- For smaller windows 5-10kb, I recommend binning in 25bp. For widows > 20kb, I recommend 50bp bins.
- Tables for meta plots exproted as a .rda file.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest         [required]\n")
    cat("--type        : Peaks, Tss, or Tes                                            [required]\n")
    cat("--upStream    : distance upstream of the region to take                       [default = 2.5kb]\n")
    cat("--downStream  : distance downstream of the region to take                     [default = 2.5kb]\n")
    cat("--bins        : number to bin coverage the window should be divisable by this [default = 25bp]\n")
    cat("--plusBw      : bigWigFile for plus strand                                    [required]\n")
    cat("--minusBw     : bigWigFile for minus strand                                   [required]\n")    
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                         [default = hg19]\n")    
    cat("--bwFile      : path to bigWig File                                           [required]\n")
    cat("--numCores    : number of cores to use                                        [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                 [default = basename(bigWigFile) ]\n")
    cat("--tssMax      : use the max postion of coverage in the annotated tss (0/1)     [default = 0; use annotated ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    Type       <- sub( '--type=', '', args[grep('--type=', args)] )
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    Bins       <- sub( '--bins=', '', args[grep('--bins=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    plusBw     <- sub( '--plusBw=', '', args[grep('--plusBw=', args)])
    minusBw    <- sub( '--minusBw=', '', args[grep('--minusBw=', args)])
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    tssMax     <- sub( '--tssMax=', '',args[grep('--tssMax=',args)])
}

## for debugging
#regions <- "tables/Pol2_293T_DMSO_817_rep1.filteredProteinCodingTx.rda"
#regions <- "tables/PRO_0h_Aux_PAF1AID_DLD1_846tss500.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#upStream <- 1000
#downStream <- 1000
#Bins <- 0
#Type <- "Tss"
#tssMax <- 1
#plusBw <- "/projects/b1025/arw/analysis/kevin/SEC/data_proseq/PRO_468_293T_836.plus.bw"
#minusBw <- "/projects/b1025/arw/analysis/kevin/SEC/data_proseq/PRO_468_293T_836.minus.bw"
#plusBw           <- "/projects/b1025/arw/analysis/fei/paf1/data_proseq/PRO_1h_Aux_PAF1AID_DLD1_846.plus.bw"
#minusBw          <- "/projects/b1025/arw/analysis/fei/paf1/data_proseq/PRO_1h_Aux_PAF1AID_DLD1_846.minus.bw"
#Cores <- 10
#outName <-paste("tables/heatmapsMaxTss", sub(".plus.bw", "", basename(plusBw)), sep="/")
#setwd("../")
#setwd("../../../kevin/SEC/")


## make output directory if needed
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bw", "", basename(bwFile))
}

if (identical(upStream,character(0))){
    upStream <- 1000
}else{
    upStream <- as.numeric(upStream)
}

if (identical(downStream,character(0))){
    downStream <- 1000
}else{
    downStream <- as.numeric(downStream)
}

if (identical(Bins,character(0))){
    Bins <- 25
}else{
    Bins <- as.numeric(Bins)
}

if (identical(tssMax,character(0))){
    tssMax <- 0
}else{
    tssMax <- as.numeric(tssMax)
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
}else{
    print("use annotated Tss")
    Model          <- get(load(regions))
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
    Model$length   <- width(Model)
}

if( Type=="Tss" ){
    print("model is Tss")
    Model.win      <- promoters(Model, upstream=upStream, downstream=downStream)
    Model.win$name <- Model.win$gene_id
    head(ranges(Model.win))
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    Model.win      <- promoters(tes, upstream=upStream, downstream=upStream)
    Model.win$name <- Model.win$gene_id
    head(ranges(Model.win))
}else if( Type=="Peaks" ){
    print("model is Peaks")
    Model.win      <- promoters(resize(Model, width=1, fix="center"), upstream=upStream, downstream=downStream)
    Model.win$name <- Model.win$name
    head(ranges(Model.win))
}

idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
if (length(idx) != 0L){
    print( paste("out of bounds removed", idx, names(Model.win[idx])) )
    Model.win <- Model.win[-idx]
}

seqlevels(Model.win,force=TRUE) <- seqlevels(Model.win)[grep("chrM",seqlevels(Model.win), invert=TRUE)]


matBin <- function(model){
    ## need to make a matrix for the sense and antisense strands for each region
    fname <- sub("$", paste0("_", Type, "up",upStream,"down", downStream, "bins",Bins, ".rda"), outName)
    plus      <- import.bw(plusBw,RangedData=FALSE,selection = BigWigSelection(model))
    minus     <- import.bw(minusBw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    plus.cov  <- coverage(plus, weight='score')
    minus.cov <- coverage(minus, weight='score') *-1    
    ## get coverage for the sense in defined window
    sense.cov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                r <- minus.cov[[seqname]][start:end]
                r <- rev(r)
                return(r)
            }else{
                r <- plus.cov[[seqname]][start:end]
                return(r)
            }
        }           
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })       
    sense.mat <- do.call(rbind, mclapply(sense.cov, as.numeric,mc.cores=Cores))
    sense.mat <- data.frame(sense.mat)
    ## get coverage for the sense in defined window
    anti.cov  <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            if(strand == '-'){
                r <- plus.cov[[seqname]][start:end]
                r <- rev(r)
                return(r)
            }else{
                r <- minus.cov[[seqname]][start:end]
                return(r)
            }
        }           
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })       
    anti.mat <- do.call(rbind, mclapply(anti.cov, as.numeric,mc.cores=Cores))
    anti.mat <- data.frame(anti.mat)
    ## bin the matrix if needed
    if( Bins > 0){
        binMatrix <- function(cov){
            window.cov <- function(row){
                window <- as.integer(ncol(cov)/Bins)
                window.coverage <- lapply(0:(window-1), function(jump)
                    rowMeans(row[as.numeric(jump*Bins+1):as.numeric((jump*Bins)+Bins)])
                    )
                t(as.matrix(unlist(window.coverage)))
            }
            win <- mclapply(1:nrow(cov), function(i)
                window.cov(cov[i,]),mc.cores=Cores)    
            bin.mat <- do.call(rbind, mclapply(win, as.numeric, mc.cores=Cores))        
            df <- data.frame(bin.mat)
            rownames(df) <- model$name
            df
        }
        sense.df <- binMatrix(cov=sense.mat)
        anti.df  <- binMatrix(cov=anti.mat)
        save(sense.df, file=sub(".rda", ".sense.rda", fname))
        save(anti.df,  file=sub(".rda", ".anti.rda", fname))
    }else{
        rownames(sense.mat) <- model$name
        rownames(anti.mat)  <- model$name
        sense.df            <- sense.mat
        anti.df             <- anti.mat
        save(sense.df, file=sub(".rda", ".sense.rda", fname))
        save(anti.df,  file=sub(".rda", ".anti.rda", fname))
    }
    print("done")
}


matBin(model=Model.win)
