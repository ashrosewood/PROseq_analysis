args <- commandArgs()

help <- function(){
    cat("createROtracks.R :
- Author: Ashley Woodfin
- Create tracks for run on sequencing such as GROseq and PROseq or TT-seq\n")
    cat("Usage: \n")
    cat("--bamFile    : Sample bam file (must be a full path)                         [required]\n")    
    cat("--outName    : Prefix for your output file                                   [default = sample]\n")    
    cat("--assembly   : Genome build hg19, mm9, mm10, dm3 ect.                        [default = hg19]\n")
    cat("--sepStrands : If the strands should be separated or kept in one track (0/1) [default = 1]\n")
    cat("--flipStrand : If the strand should be flipped. Depends where the primer is. [default = 1]\n")
    cat("               (typicllay GROseq no ;PROseq yes) (0/1)                        \n")
    cat("--extLen     : number of bases to extend the reads to                        [default = 0]\n")
    cat("--threePrime : For PROseq make tracks of the 3' end of the reads (0/1)       [default = 0]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    bamFile     <- sub('--bamFile=', '', args[grep('--bamFile=', args)])
    outName     <- sub('--outName=', '', args[grep('--outName=', args)])
    assembly    <- sub('--assembly=', '', args[grep('--assembly=', args)])
    sepStrands  <- sub('--sepStrands=', '', args[grep('--sepStrands=', args)])
    flipStrand  <- sub('--flipStrand=', '', args[grep('--flipStrand=', args)])
    extLen      <- sub('--extLen=', '', args[grep('--extLen=', args)])
    threePrime  <- sub('--threePrime=', '', args[grep('--threePrime=', args)])
}

print(bamFile)

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
    seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_|chrM",seqlevels(bd), invert=TRUE)]
    print(seqlevels(bd))
    if( flipStrand > 0 ){
        cat("change the strand:", BF, sep="\n")
        strand(bd) <- ifelse(strand(bd) == '+', '-', '+')
    }
    cat("convert to GRanges\n")
    mygr <- as(bd,"GRanges")
    if (extLen > 0){
        cat("extend reads:", BF, sep="\n")
        mygr <- resize(mygr, extLen, fix='start')
    }
    if (sepStrands > 0){        
        cat("getting coverage for separate strands\n")
        ## get plus coverage                                                             
        plus                 <- coverage(mygr[strand(mygr) == "+"])
        plus.rpm             <- plus*(1e6/length(bd))
        seqlengths(plus.rpm) <- seqlengths(organism)[names(plus.rpm)]
        ## get minus coverage
        minus                 <- coverage(mygr[strand(mygr) == "-"])
        minus.rpm             <- minus*(-1e6/length(bd))
        seqlengths(minus.rpm) <- seqlengths(organism)[names(minus.rpm)]
        ## set the outfile name
        plus.outfile  <- sub("$", ".plus.bw", outName)
        minus.outfile <- sub("$", ".minus.bw", outName)
        ## export rpm to bigWig
        cat(paste("exporting to plus bigwig", plus.outfile, "\n", sep="\t"))
        export.bw(plus.rpm, plus.outfile)
        cat(paste("exporting to minus bigwig", minus.outfile, "\n", sep="\t"))
        export.bw(minus.rpm, minus.outfile)      
        cat("export complete:", BF, sep="\n")
        if (threePrime > 0){
            cat( "take last position of read:", BF, sep="\n" )     
            mygr                                                    <- as(bd,"GRanges")
            start(mygr[paste(as.data.frame(mygr)[,"strand"])=='+']) <- end(mygr[paste(as.data.frame(mygr)[,"strand"])=='+'])
            end(mygr[paste(as.data.frame(mygr)[,"strand"])=='-'])   <- start(mygr[paste(as.data.frame(mygr)[,"strand"])=='-'])
            plus                                                    <- coverage(mygr[strand(mygr) == "+"])
            plus.rpm                                                <- plus*(1e6/length(bd))
            seqlengths(plus.rpm)                                    <- seqlengths(organism)[names(plus.rpm)]
            ## get minus coverage
            minus                                                   <- coverage(mygr[strand(mygr) == "-"])
            minus.rpm                                               <- minus*(-1e6/length(bd))
            seqlengths(minus.rpm)                                   <- seqlengths(organism)[names(minus.rpm)]
            ## set the outfile name
            plus.outfile                                            <- sub("$", ".3prime.plus.bw", outName)
            minus.outfile                                           <- sub("$", ".3prime.minus.bw", outName)
            ## export rpm to bigWig
            cat(paste("exporting plus bigwig", plus.outfile, "\n", sep="\t"))
            export.bw(plus.rpm, plus.outfile)
            cat(paste("exporting to minus bigwig", minus.outfile, "\n", sep="\t"))
            export.bw(minus.rpm, minus.outfile)      
            cat("export complete:", BF, sep="\n")
        }
    }else{
        cat("getting coverage for both strands\n")
        cov <- coverage(mygr)
        rpm <- cov*(1e6/length(bd))
        seqlengths(rpm) <- seqlengths(organism)[names(rpm)]
        ## export rpm to bigWig
        outfile <- sub("$", ".bw", outName)
        cat(paste("exporting to bigwig", outfile, "\n", sep="\t"))
        export.bw(rpm, outfile)
        cat("export complete:", BF, sep="\n")
          
    }
}

# for each element of our vector, call the bam2bw function
mclapply(bamFile,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
