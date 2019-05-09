library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(readr)


myArgs <- commandArgs(trailingOnly = T)
# path to your idat files
path <- myArgs[1]
pathBValue <- myArgs[2]

if(length(myArgs) != 2){
    stop("Usage: Rscript [path_to_450K_idat] [path_to_beta_value]")
}

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")

cat("Start reading in ", path, "\n")
rgSet <- read.metharray.exp(base = path, recursive = F, force=TRUE)
# rgSet <- read.metharray.exp(base = path, recursive = F)
sampleName <- sapply(colnames(rgSet), FUN = function(X){unlist(strsplit(X, split = "_"))[1]})
colnames(rgSet) <- sampleName

# detect P value
detP <- detectionP(rgSet)

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
detP <- detP[,keep]
mSetSq <- preprocessQuantile(rgSet)
mSetRaw <- preprocessRaw(rgSet)

detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

cat("Filter out snps\n")
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

cat("Filter out sex chromosomes\n")
# filter out pos on sex chr
xReactiveProbes <- read.csv(file=paste(dataDirectory, "48639-non-specific-probes-Illumina450k.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
bVals <- getBeta(mSetSqFlt)

cat("Start writing to ", pathBValue, "\n")
bValsOut <- cbind(data.frame(pos = rownames(bVals)), bVals)
write_delim(x = as.data.frame(bValsOut), path = pathBValue, delim = "\t", col_names = T)
