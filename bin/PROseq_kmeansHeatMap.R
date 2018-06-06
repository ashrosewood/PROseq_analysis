args <- commandArgs()

help <- function(){
    cat("PROseq_kmeansHeatMap.R :
- Use kmeans clustering to cluster heatmaps cdts.
- This script is meant for multiple files that need to be compared to a control.
- They will be separated by cluster than sorted by rowSums in each cluster
- A fold change cdt is also made.\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                           [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T)  [optional]\n")    
    cat("--clusNum    : number of clusters for heatmap                                   [default = 4]\n")
    cat("--Control    : number of bp upstream of Tss in matrix                           [required]\n")
    cat("--Cluster    : Cluster by coverage heatmap (Coverage) or fold change (FC)       [default = FC]\n")
    cat("--outName    : prefix to your out file names (No .extention)                    [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir        <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern    <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    clusNum    <- sub( '--clusNum=', '',args[grep('--clusNum=',args)])
    Control    <- sub( '--Control=', '',args[grep('--Control=',args)])
    Cluster    <- sub( '--Cluster=', '',args[grep('--Cluster=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
}

setwd("/projects/b1025/arw/analysis/yuki/degrons/DLD1")
Dir <- "tables/heatmaps/parental_genes"
Pattern <- "PRO_4h_NELF.*Tssup50down100bins0.sense.log2FC.rda|PRO_2h_NELFeAID.*Tssup50down100bins0.sense.log2FC.rda"
clusNum=3
Cluster <- "FC"
outName <- "test"
#Control <- "Pol2_CT_12Mio_HCT116_1087"

if (identical(Cluster,character(0))){
    Cluster <- "FC"
}

if (identical(clusNum,character(0))){
   clusNum <- 4
}else{
   clusNum <- as.numeric(clusNum)
}

## show params
print(paste("clusNum", clusNum))
print(paste("Dir", Dir))
print(paste("Pattern", Pattern))
print(paste("outName", outName))
#print(paste("Control", Control))
print(paste("Cluster by", Cluster))

## select rda files for pattern of interest
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda|_Peaks.*rda|_Tss.*rda|_Tes.*rda|.log2FC.rda",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
    colnames(df) <- paste(sub(".df", "", oname), 1:ncol(df), sep=".")
    assign(oname, df)
}

SAMPLES <- ls(pattern=".df$")
SAMPLES

## sanity check must be the same genes in the same order
stopifnot(rownames(get(SAMPLES[1]))==rownames(get(SAMPLES[2])))

## put the control first and set experiments
exps <- SAMPLES[SAMPLES!=paste0(Control, ".df")]

## report the order in the data frame
df.all <- get(SAMPLES[
print(paste("control sample 1", paste0(Control, ".df")))
for( i in 1:length(exps) ){
    print(paste("sample", i+1, exps[i]))
    df.all <- cbind(df.all, get(exps[i]))
}

## make the fold change heat map
fc <- do.call(cbind, lapply(1:length(exps), function(i){
    print(paste("FC sample", i, exps[i]))
    exp         <- as.matrix(get(exps[i]))
    con         <- get(SAMPLES[SAMPLES==paste0(Control, ".df")])
    pseudo      <- min(c(min(exp[!exp == 0]), min(con[!con == 0])))
    exp         <- exp + pseudo
    con         <- con + pseudo
    FC          <- log2(exp/con)
    FC
}))

df.all <- do.call(cbind, lapply(1:length(SAMPLES), function(i){
    get(SAMPLES[i])
}))






##################
## define clusters
##################
if( Cluster == "Coverage"){
    ## cluster by coverage heatmap
    print( paste("cluster by kmeans for", clusNum, "clsuters") )
    ## kmeans parameters
    go.paras       <- list(knc=clusNum, max.iter=20, nrs=30)
    ## force kmeans to use the same random seed everytime
    set.seed(20)    
    km             <- kmeans(df.all, centers=go.paras$knc, 
                             iter.max=go.paras$max.iter, nstart=go.paras$nrs)
    ## get cluster numbers and sort by Control occupancy within each cluster
    clus              <- data.frame(km$cluster)
    clus$totalCov     <- rowSums( get(SAMPLES[SAMPLES==paste0(Control, ".df")]) )
    clus.or           <- clus[with(clus, order(km.cluster, -totalCov)), ]
    clus.or$gene_name <- rownames(clus.or)
    write.table(clus.or
               ,file=sub("$", paste(".kmeans", clusNum, "clusterIds.txt", sep="."), outName)
               ,sep="\t",row.names=FALSE, quote=FALSE)
}else{
    ## cluster by FC heatmap
    print( paste("cluster by kmeans for", clusNum, "clsuters") )
    ## kmeans parameters
    go.paras       <- list(knc=clusNum, max.iter=20, nrs=30)
    ## force kmeans to use the same random seed everytime
    set.seed(20)    

    km             <- kmeans(df.all, centers=go.paras$knc, 
                             iter.max=go.paras$max.iter, nstart=go.paras$nrs)

    ## get cluster numbers and sort by Control occupancy within each cluster
    clus              <- data.frame(km$cluster)
    clus$totalCov     <- rowSums( get(SAMPLES[1]) )
    clus.or           <- clus[with(clus, order(km.cluster, -totalCov)), ]
    clus.or$gene_name <- rownames(clus.or)
    write.table(clus.or
               ,file=sub("$", paste(".kmeans", clusNum, "clusterIds.txt", sep="."), outName)
               ,sep="\t",row.names=FALSE, quote=FALSE)


}

## report the number in each cluster
for ( x in 1:clusNum){
    print(paste("cluster", x, sum(clus.or$km.cluster==x), "genes")) 
}    

write.table(cbind(UID=rownames(clus.or[clus.or$km.cluster==1,])
                 ,NAME=rownames(clus.or[clus.or$km.cluster==1,])
                 ,df.all[rownames(clus.or[clus.or$km.cluster==1,]),]
                  ), ,file="test.clus.cdt"
           ,sep="\t",row.names=FALSE)

write.table(cbind(UID=rownames(clus.or[clus.or$km.cluster==2,])
                 ,NAME=rownames(clus.or[clus.or$km.cluster==2,])
                 ,df.all[rownames(clus.or[clus.or$km.cluster==2,]),]
                  ), ,file="test.clus2.cdt"
           ,sep="\t",row.names=FALSE)

write.table(cbind(UID=rownames(clus.or[clus.or$km.cluster==3,])
                 ,NAME=rownames(clus.or[clus.or$km.cluster==3,])
                 ,df.all[rownames(clus.or[clus.or$km.cluster==3,]),]
                  ), ,file="test.clus3.cdt"
           ,sep="\t",row.names=FALSE)

plot(colMeans(get(SAMPLES[1])[rownames(clus.or[clus.or$km.cluster==2,]),]))


##############################
## write combined coverage cdt
##############################
all.cdt <- data.frame()
## add arbitrary divider for the combined heatmap
for ( x in 1:clusNum){
    print(paste("Combined heat cluster",  x))
    cdt <- cbind(UID=rownames(clus.or[clus.or$km.cluster==x,])
                ,NAME=rownames(clus.or[clus.or$km.cluster==x,])
                ,df.all[rownames(clus.or[clus.or$km.cluster==x,]),]
                 )
    fake <- data.frame(cbind(UID=paste0("clus",x)
                            ,NAME=paste0("clus",x)
                                        #,t(data.frame(num=rep(ceiling(max(df.all)),ncol(df.all))))
                             ,data.frame( matrix(0, ncol = ncol(df.all), nrow = 25) )
                             ))
    colnames(fake) <- colnames(cdt)
    rownames(fake) <- sub("^", paste0("clus",x, "."), rownames(fake))
    all.cdt <- rbind(all.cdt, cdt, fake)   
}
## write the cdt for each sample
write.table(all.cdt
           ,file=sub("$", paste(".kmeans", clusNum, "cdt", sep="."), outName)
           ,sep="\t",row.names=FALSE)

print(paste("coverage heatmap saved as", sub("$", paste(".kmeans", clusNum, "cdt", sep="."), outName)))













########################
## write combined FC cdt
########################
write.table(
           ,file=sub("$", paste(".kmeans", clusNum, "log2FC.cdt", sep="."), outName)
           ,sep="\t",row.names=FALSE)    



fc.cdt                  <- data.frame()
for( x in 1:clusNum){
    print(paste("Combined FC heat cluster",  x))
    cdt <- cbind(UID=rownames(clus.or[clus.or$km.cluster==x,])
                ,NAME=rownames(clus.or[clus.or$km.cluster==x,])
                ,fc[rownames(clus.or[clus.or$km.cluster==x,]),]
                 )
    fake <- data.frame(cbind(UID=paste0("clus",x)
                            ,NAME=paste0("clus",x)
                            ,data.frame( matrix(0, ncol = ncol(fc), nrow = 25) )
                             ))
    colnames(fake) <- colnames(cdt)
    rownames(fake) <- sub("^", paste0("clus",x, "."), rownames(fake))
    fc.cdt <- rbind(fc.cdt, cdt, fake)
}
write.table(fc.cdt
           ,file=sub("$", paste(".kmeans", clusNum, "log2FC.cdt", sep="."), outName)
           ,sep="\t",row.names=FALSE)    

print(paste("coverage heatmap saved as", sub("$", paste(".kmeans", clusNum, "cdt", sep="."), outName)))

print("Done")
