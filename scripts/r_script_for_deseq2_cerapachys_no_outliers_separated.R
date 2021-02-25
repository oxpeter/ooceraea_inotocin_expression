# load DESeq2 library and functions:
library("DESeq2")

##########################################################################
########## Oocerea subjects (genes) ###################################
##########################################################################


# set location of HTSeq count files:
htseq_dir <- "../htseq"

# extract sample files:
sampleFilesP <- grep("P.*(_SP|_FL).*htseq.gene.out", list.files(htseq_dir), value=TRUE)
sampleFilesR <- grep("R.*(_SP|_FL).*htseq.gene.out", list.files(htseq_dir), value=TRUE)
sampleFilesP
sampleFilesR

names <- c("P", "R") # for saving results of each loop consistently with HTSeq file names
i <- 0
for(sampleFiles in c(sampleFilesP, sampleFilesR)){

i <- i+1
sampleset <- names[i]

# assign phase and batch to each sample, based on filename:
sampleCondition <- sub(".*(_SP|_FL).*", "\\1", sampleFiles)
sampleCondition

sampleBatch <- sub("[PM]_(B[0-9]).*", "\\1", sampleFiles)
sampleBatch


# construct table, including design (incorporates batch effect in model)
sampleTable <- data.frame(sampleName = sampleFiles,
                           fileName = sampleFiles,
                           condition = sampleCondition,
			   batch = sampleBatch)

### PERFORM TESTS ###                          
# standard
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = htseq_dir,
                                        design= ~ batch + condition)
dds
dds_std <- DESeq(dds, minReplicatesForReplace=5)

# look at results of tests:
resQF <- results(dds_std, alpha=0.05, contrast=c("condition", "_SP", "_FL"))


# produce and save the variance adjusted count values for each gene and save to file:
vsd <- varianceStabilizingTransformation(dds_std)
vstMat <- assay(vsd)
write.table(as.data.frame(vstMat), file=paste("../DESeqResults/Cerapachys.genes.all_samples", sampleset, "vst.tbl", sep = "."))


# order results on adjusted pvalue:
resQFOrdered <- resQF[order(resQF$padj),]
head(resQFOrdered)
summary(resQF)



# count number of genes with adjusted pvalue less than 0.05
SigQF <- subset(resQFOrdered, padj < 0.05)
length(SigQF$padj)

# graph all variance adjusted genes:
plotMA(resQF, main="statary and foraging Cerapachys Samples (Genes)", ylim=c(-2,2))

# save pvalue and log fold change results for all genes to file:
write.table(as.data.frame(resQFOrdered), file=paste("../DESeqResults/Cerapachys.genes.all_samples", sampleset, "pvalues.out", sep = ".") )


}
