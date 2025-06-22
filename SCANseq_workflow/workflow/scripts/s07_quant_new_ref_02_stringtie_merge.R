# conda activate renv
# R version 4.2.2
suppressMessages({
	library(data.table)
	library(parallel)
	library(dplyr)
	options(stringsAsFactors = FALSE)
})

#samplelist <- fread(file='config/samples_5400.txt', header=F)$V1
#filelist <- file.path("results/07_quant_new_ref", samplelist, "quant_new_ref_stringtie.gtf")
samplelist <- snakemake@params[["sample"]]
filelist <- snakemake@input[["gtf"]]
# samplelist <- samplelist[1:10]
# filelist <- filelist[1:10]

# Function to process each file and return a data.table
process_file <- function(file, sample) {
	gtf <- read.table(file, header = FALSE, sep = "\t", comment.char="#")
	gtf <- gtf[which(gtf$V3 %in% "transcript"),]
	cov <- data.table(
		TranID = sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x) x[4]), 
		cov = as.numeric(sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x) x[6]))
	)
	print(sample)
	setnames(cov, c("TranID", sample))
	return(cov)
}

# Parallel processing setup
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("filelist", "samplelist", "process_file"))
clusterEvalQ(cl, library(data.table))

# Process all files in parallel 
data_list <- parLapply(cl, seq_along(filelist), function(i) {
	process_file(filelist[[i]], samplelist[[i]])
})
stopCluster(cl)

# Merge all data tables
merge_cov <- Reduce(function(x, y) merge(x, y, by = "TranID", all = TRUE), data_list)

# Write the merged data table to file
fwrite(merge_cov, file = snakemake@output[["counts"]], quote = FALSE, row.names = FALSE, sep = "\t")


## old version
# for (i in c(1:length(filelist))){
#     #print(i)
#     sample <- samplelist[[i]]
#     print(sample)
#     gtf <- read.table(file=filelist[[i]], header = F,sep="\t", comment.char="#")
#     gtf <- gtf[which(gtf$V3 %in% "transcript"),]
#     cov <- data.frame(
#         TranID = sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x){x[4]}),
#         cov = as.numeric(sapply(strsplit(gsub(";", "", gtf$V9), split = " "), function(x){x[6]}))
#     )
#     colnames(cov) <- c("TranID", sample)
    
#     if (i==1) {
#         merge_cov <- cov
#     } else {
#         merge_cov <- merge(merge_cov, cov, by="TranID", all=T)
#     }
# }
# write.table(merge_cov, file=snakemake@output[["counts"]], quote=FALSE, row.names=FALSE, sep = "\t")
