args <- commandArgs()

help <- function(){
    cat("createROtracks.R :
- Create tracks for run on sequencing such as GROseq and PROseq or TT-seq\n")
    cat("Usage: \n")
    cat("--bamFile    : Sample bam file (must be a full path)                           [required]\n")    
    cat("--outName    : Prefix for your output file                                     [default = sample]\n")    
    cat("--assembly   : Genome build hg19, mm9, mm10, dm3 ect.                          [default = hg19]\n")
    cat("--sepStrands : If the strands should be separated or kept in one track (0/1)   [default = 1]\n")
    cat("--flipStrand : If the strand should be flipped. Depends where the primer is.   [default = 1]\n")
    cat("               (typicllay GROseq no ;TTseq, RNAseq, PROseq yes) (0/1)                       \n")
    cat("--extLen     : number of bases to extend the reads to                          [default = 0]\n")
    cat("--noStrand   : replace strand with * (ChIP-seq) (0/1)                          [default = 0]\n")
    cat("--threePrime : Only report 1 position at the 3' end of the read (PRO-seq)(0/1) [default = 0]\n")
    cat("--fivePrime  : Only report 1 position at the 5' end of the read (0/1)          [default = 0]\n")
    cat("--SpikeIn    : Fraction of control spike over experiment to normalize reads    [default = 1]
                         if your control spikein total reads are 3549426 and experiment
                         are 4812918, your fraction is 0.737479 \n")    
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    bamFile     <- sub('--bamFile=', '',    args[grep('--bamFile=', args)])
    outName     <- sub('--outName=', '',    args[grep('--outName=', args)])
    assembly    <- sub('--assembly=', '',   args[grep('--assembly=', args)])
    sepStrands  <- sub('--sepStrands=', '', args[grep('--sepStrands=', args)])
    flipStrand  <- sub('--flipStrand=', '', args[grep('--flipStrand=', args)])
    extLen      <- sub('--extLen=', '',     args[grep('--extLen=', args)])
    noStrand    <- sub('--noStrand=', '',   args[grep('--noStrand=', args)])
    threePrime  <- sub('--threePrime=', '', args[grep('--threePrime=', args)])
    fivePrime   <- sub('--fivePrime=', '',  args[grep('--fivePrime=', args)])
    SpikeIn     <- sub('--SpikeIn=', '',  args[grep('--SpikeIn=', args)])
}

print(bamFile)

if (identical(SpikeIn,character(0))){
   SpikeIn <- 1
}else{
    SpikeIn <- as.numeric(SpikeIn)
}

if (identical(outName,character(0))){
   outName <- sub(".bam", "", bamFile)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(sepStrands,character(0))){
   sepStrands <- 1
}else{
   sepStrands <- as.numeric(sepStrands)
}

if (identical(flipStrand,character(0))){
    flipStrand <- 1
}else{
    flipStrand <- as.numeric(flipStrand)
}

if (identical(extLen,character(0))){
   extLen <- 0
}else{
   extLen <- as.numeric(extLen)
}

if (identical(threePrime,character(0))){
   threePrime <- 0
}else{
   threePrime <- as.numeric(threePrime)
}

if (identical(fivePrime,character(0))){
   fivePrime <- 0
}else{
   fivePrime <- as.numeric(fivePrime)
}

if (identical(fivePrime,character(0))){
   fivePrime <- 0
}else{
   fivePrime <- as.numeric(fivePrime)
}

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster

bam2bw <- function(BF,organism){
    cat("opening:", BF, sep="\n")
    bd                       <- readGAlignments(BF)
    cat("remove unassembled and viral chromosomes:", BF, sep="\n") 
    seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_|EBV",seqlevels(bd), invert=TRUE)]
    print(seqlevels(bd))
    if( flipStrand > 0 ){
        cat("change the strand:", BF, sep="\n")
        strand(bd)           <- ifelse(strand(bd) == '+', '-', '+')
    }
    if( noStrand > 0 ){
        strand(bd)           <- '*'
    }
    cat("convert to GRanges\n")
    mygr                    <- as(bd,"GRanges")
    if (extLen > 0){
        cat("extend reads:", BF, sep="\n")
        mygr                 <- resize(mygr, extLen, fix='start')
    }
    if (threePrime > 0){
        cat("take three prime position of the reads:", BF, sep="\n")
        mygr                 <- resize(mygr, width=1, fix='end')
        outName              <- sub("$", ".3prime", outName)
    }
    if ( fivePrime > 0 ){
        cat("take five prime position of the reads:", BF, sep="\n")
        mygr                 <- resize(mygr, width=1, fix='start')
        outName              <- sub("$", ".5prime", outName)
    }
    if (sepStrands > 0){        
        cat("getting coverage for separate strands\n")
        ## get plus coverage                                                             
        plus                  <- coverage(mygr[strand(mygr) == "+"])
        plus.rpm              <- plus*(1e6/length(bd))
        seqlengths(plus.rpm)  <- seqlengths(organism)[names(plus.rpm)]
        ## get minus coverage
        minus                 <- coverage(mygr[strand(mygr) == "-"])
        minus.rpm             <- minus*(-1e6/length(bd))
        seqlengths(minus.rpm) <- seqlengths(organism)[names(minus.rpm)]
        ## set the outfile name
        plus.outfile          <- sub("$", ".rpm.plus.bw", outName)
        minus.outfile         <- sub("$", ".rpm.minus.bw", outName)
        ## export rpm to bigWig
        cat(paste("exporting r.p.m. plus bigwig", plus.outfile, "\n", sep="\n"))
        export.bw( plus.rpm, plus.outfile )
        cat(paste("exporting r.p.m. minus bigwig", minus.outfile, "\n", sep="\n"))
        export.bw( minus.rpm, minus.outfile )      
        if( SpikeIn != 1 ){
            print("making SpikeIn normalized tracks")
            ## get plus coverage
            plus.si               <- plus*SpikeIn
            seqlengths(plus.si)   <- seqlengths(organism)[names(plus.si)]
            ## get minus coverage
            minus.si             <- minus*(-SpikeIn)
            seqlengths(minus.si) <- seqlengths(organism)[names(minus.si)]
            ## set the outfile name
            plus.si.outfile      <- sub("$", ".spikeNorm.plus.bw", outName)
            minus.si.outfile     <- sub("$", ".spikeNorm.minus.bw", outName)
            ## export rpm to bigWig
            cat(paste("exporting spikeNorm plus bigwig", plus.si.outfile, "\n", sep="\n"))
            export.bw( plus.si, plus.si.outfile )
            cat(paste("exporting spikeNorm minus bigwig", minus.si.outfile, "\n", sep="\n"))
            export.bw( minus.si, minus.si.outfile ) 
        }
        cat("export complete:", BF, sep="\n")
    }else{       
        cat("getting coverage for both strands\n")
        cov                  <- coverage(mygr)
        rpm                  <- cov*(1e6/length(bd))
        seqlengths(rpm)      <- seqlengths(organism)[names(rpm)]
        ## export rpm to bigWig
        outfile              <- sub("$", ".bw", outName)
        cat( paste("exporting r.p.m. bigwig", outfile, "\n", sep="\t") )
        export.bw( rpm, outfile )
        if( SpikeIn != 1 ){
            print("making SpikeIn normalized tracks")
            cov.si             <- cov*SpikeIn
            seqlengths(cov.si) <- seqlengths(organism)[names(cov.si)]
            ## export rpm to bigWig
            outfile.si         <- sub("$", ".spikeNorm.bw", outName)
            cat( paste("exporting spikeNorm bigwig", outfile.si, "\n", sep="\n") )
            export.bw( cov.si, outfile.si )
        }
        cat( "export complete:", BF, sep="\n" )          
    }
}

# for each element of our vector, call the bam2bw function
mclapply(bamFile,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
