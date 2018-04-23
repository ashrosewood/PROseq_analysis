args <- commandArgs()

help <- function(){
    cat("ChIP_callDifferentialPeaksHyper.R :
- This is set up for ChIP seq data with one replicate to call peaks based on a hypergeometirc p-value.
- It will calculated the number of reads that overlap the entire peak regions in the experiment and control.\n")
    cat("Usage: \n")
    cat("--Control    : Path to control bam file                                 [required]\n")
    cat("--Experiment : Path to experiment bam file                              [required]\n")
    cat("--assembly   : genome (hg19, mm9, mm10, dm3, sacCer3)                   [hg19]\n")
    cat("--PeakFile   : Path to peaks (rda or bed) must have a name in column    [required]\n")
    cat("--outName    : Path to output file (no file extension ie .txt)          [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Control     <- sub( '--Control=', '', args[grep('--Control=', args)] )
    Experiment  <- sub( '--Experiment=', '', args[grep('--Experiment=', args)] )
    assembly    <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    PeakFile    <- sub( '--PeakFile=', '', args[grep('--PeakFile=', args)] )
    outName     <- sub( '--outName=', '',args[grep('--outName=',args)] )
}

## set cores to 8 if not set
if (identical(assembly,character(0))){
    assembly <- "hg19"
}

print(paste("Control:", Control))
print(paste("Experiment:", Experiment))
print(paste("Peaks:", PeakFile))

# example for debugging
##setwd("/projects/b1025/arw/analysis/fei/paf1")
##Control    <- "data_chipseq/H3K4me3_0h_Aux_PAF1AID_DLD1_1015.bam"
##Experiment <- "data_chipseq/H3K4me3_1D_Aux_PAF1AID_DLD1_1015.bam"
##PeakFile   <- "data_chipseq/H3K4me3_0h_Aux_PAF1AID_DLD1_1015.macsPeaks.bed"
##outName    <- "tables/peak_calling/H3K4me3_1D_vs_0h_Aux_PAF1AID_DLD1_1015"

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)

print(assembly)
if (assembly == "hg19") organismStr <- "Hsapiens"
if (assembly == "mm9") organismStr <- "Mmusculus"
if (assembly == "mm10") organismStr <- "Mmusculus"
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

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

param           <- ScanBamParam(what='mapq')

if( length(grep(".bed$", PeakFile)) > 0 ){
    peaks           <- import(PeakFile)
    seqinfo(peaks)  <- seqinfo(organism)[seqlevels(peaks)]
}else{
    peaks           <- get(load(PeakFile))
    seqinfo(peaks)  <- seqinfo(organism)[seqlevels(peaks)]
    peaks
}

hyperTable <- function(Control,Experiment){
    ##-------------read in control and experiment bam files-----------##
    con             <- readGAlignments(Control,param=param)  
    seqinfo(con)    <- seqinfo(organism)[seqlevels(con)]
    exp             <- readGAlignments(Experiment,param=param)  
    seqinfo(exp)    <- seqinfo(organism)[seqlevels(exp)]
    con.all         <- length(con)
    exp.all         <- length(exp)
    print(paste("Control total:", con.all))
    print(paste("Experiment total:", exp.all))
    ##-------------add counts to peaks-----------##
    mcols(peaks, level="within")[,sub(".bam", ".counts", basename(Control))]    <- countOverlaps(peaks, con)
    mcols(peaks, level="within")[,sub(".bam", ".counts", basename(Experiment))] <- countOverlaps(peaks, exp)
    mcols(peaks, level="within")[,sub(".bam", ".total", basename(Control))]     <- con.all
    mcols(peaks, level="within")[,sub(".bam", ".total", basename(Experiment))]  <- exp.all    
    ##-----------run phyper and make tables for plots---------##
    Dir  <- dirname(outName)
    if(!(file.exists( Dir ))) {
        dir.create(Dir,FALSE,TRUE)  
    }
    df                                              <- as.data.frame(peaks) 
    ## calc library normilization factor
    NORM                                            <- con.all / exp.all
    ## calc fold chage with pseudo count
    Con                                             <- df[,sub(".bam", ".counts", basename(Control))]
    Exp                                             <- df[,sub(".bam", ".counts", basename(Experiment))]
    df[,sub(".bam", ".rpkm", basename(Control))]    <- Con/( df$width/1000 * con.all/10^6 )
    df[,sub(".bam", ".rpkm", basename(Experiment))] <- Exp/( df$width/1000 * con.all/10^6 )
    df$log2FC                                       <- log2(NORM*(Exp + 1/NORM) / (Con + 1)) 
    ## set up variables
    m <- as.numeric(Exp)
    n <- as.numeric(Exp+Con)
    M <- as.numeric(exp.all)
    N <- as.numeric(exp.all+con.all)
    ## calculate enriched and depleted hypergeometric pvals
    df$p.over <- signif(phyper(m-1,n,N-n,M,lower.tail=FALSE))
    df$p.over.adj <- p.adjust(df$p.over, method="BH")
    ## calculate enriched and depleted hypergeometric pvals
    df$p.under <- signif(phyper(m,n,N-n,M,lower.tail=TRUE))
    df$p.under.adj <- p.adjust(df$p.under, method="BH")  
    ## write table
    write.table(df, file=paste(outName, ".phyper.txt", sep=""), sep="\t", ,col.names=TRUE,row.names=FALSE, quote=FALSE)
    print("done")
}

hyperTable(Control=Control, Experiment=Experiment)
