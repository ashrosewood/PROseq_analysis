args <- commandArgs()

help <- function(){
    cat("PolII_makeMetaGenePlot.R :
- Make meta gene plot from metaGeneMatrix.R and save as pdf\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                          [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T) [optional]\n")    
    cat("--outName    : prefix for plot name. Type will be appended to this in script   [required]\n")
    cat("--cols       : need the same number as samples separated by comma              [default = viridis palette]\n")
    cat("--upStream   : number of bp upstream of Tss in matrix                          [required]\n")
    cat("--downStream : number of bp upstream of Tss in matrix                          [required]\n")
    cat("--Height     : plot height if you do not want to be set the max                [default = max of colMeans]\n")
    cat("--Label      : Alternate plot label if not Type (example: Pause_site)          [default = Type]
                         Type is used to grep the files. If file name match use Type.                   \n")    
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir  <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern    <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols       <- sub( '--cols=', '',args[grep('--cols=',args)])
    upStream   <- sub( '--upStream=', '',args[grep('--upStream=',args)])
    downStream <- sub( '--downStream=', '',args[grep('--downStream=',args)])
    Height     <- sub( '--Height=', '',args[grep('--Height=',args)])
    Label      <- sub( '--Label=', '',args[grep('--Label=',args)])
}

if (identical(Label,character(0))){
   Label <- "Tss"
}

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

Dir
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda|.tssMax.metaGene.rda",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
    colnames(df) <- paste(sub(".df", "", oname), 1:ncol(df), sep=".")
    assign(oname, df)
}

SAMPLES <- ls(pattern=".df$")
SAMPLES

## make a data frame with matrix column averages
df <- do.call(cbind, lapply(SAMPLES, function(x){
    avg           <- as.data.frame( colMeans(get(x) ))
    colnames(avg) <- sub(".df", "", x)
    avg
  }
))
df$x <- seq(1,nrow(df))
print(head(df))
df.long <- melt(df, id.vars = "x")

mround <- function(x,base){ 
    base*ceiling(x/base) 
} 

if( max(df[,grep("^x$",names(df),invert=TRUE)]) > 0.5 ){
    Ymax <- mround(ceiling(max(df[,grep("^x$",names(df),invert=TRUE)])),2)
    Min <- 0
}else{
    Ymax <- mround(max(df[,grep("^x$",names(df),invert=TRUE)]),0.05)
    Min <- mround(min(df[,grep("^x$",names(df),invert=TRUE)]),-0.05)
}

if (identical(Height, character(0))){
    Height <- Ymax
}else{
    Height <- as.numeric(Height)
}

if (identical(cols,character(0))){
    Cols <- viridis(length(foo))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[gsub(".df|.3prime|.36bp", "", SAMPLES), "color"])
    Cols
}

pdf(file=sub("$", ".metaGenePlot.pdf", outName),width=8,height=5)
print({
    p <-

    ggplot(df.long, aes(x=x, y=value, color=variable)) +
    geom_line(size=0.5)+
    scale_color_manual("sample", values = Cols )+
    scale_x_continuous(breaks=c(0
                               ,100
                               ,100+(400/3)
                               ,100+(400/3*2)
                               ,100+(400/3*3)
                               ,600
                                )
                      ,labels=c(paste0("-",upStream)
                               ,Label
                               ,"33%"
                               ,"66%"
                               ,"TES"
                               ,paste0("+", downStream)
                                ))+
    ylab("signal")+
    xlab("")+
    ylim(0,Height)+
    ggtitle(basename(outName))+
    theme(panel.grid.major = element_blank()
         ,panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key=element_blank())+
    theme(text = element_text(size=12),
          axis.text.x = element_text(vjust=1,color="black",size=12),
          axis.text.y = element_text(color="black",size=12),
          plot.title=element_text(size=12))

})
dev.off()
