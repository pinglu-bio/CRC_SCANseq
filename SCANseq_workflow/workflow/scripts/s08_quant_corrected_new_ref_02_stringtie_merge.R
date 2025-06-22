suppressMessages({
    options(stringsAsFactors = F)
})

filelist <- snakemake@input[["gtf"]]
samplelist <- snakemake@params[["sample"]]

for (i in c(1:length(filelist))){
    print(i)
    sample <- samplelist[[i]]
    gtf <- read.table(file=filelist[[i]], header = F,sep="\t")
    gtf <- gtf[which(gtf$V3 %in% "transcript"),]
    cov <- data.frame(
        TranID = sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x){x[4]}),
        cov = as.numeric(sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x){x[8]}))
    )
    colnames(cov) <- c("TranID", sample)
    
    if (i==1) {
        merge_cov <- cov
    } else {
        merge_cov <- merge(merge_cov, cov, by="TranID", all=T)
    }
}
write.table(merge_cov, file=snakemake@output[["counts"]], quote=FALSE, row.names=FALSE, sep = "\t")