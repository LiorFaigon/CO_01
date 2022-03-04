# Champions Oncology - Exercise for bioinformatics candidates
# Submitted by Lior Faigon

# Load the data into R and make sure the count and annotation data are consistent with each other.
  #noteforself : consider removing header and rownames at reading to facilitate deseq2
raw_counts = read.csv("https://raw.githubusercontent.com/LiorFaigon/CO_01/master/counts.txt",
                      header = T,
                      row.names = 1,
                      sep = "\t")
gene_annotations = read.csv("https://raw.githubusercontent.com/LiorFaigon/CO_01/master/gene-annotation.txt",
                            header = T,
                            sep = "\t")
colnames(gene_annotations)[4] = "NAME" # fill empty colname in gene_annotations
sample_annotations = read.csv("https://raw.githubusercontent.com/LiorFaigon/CO_01/master/sample-annotation.txt",
                            header = T,
                            row.names = 1,
                            sep = "\t")

raw_counts_samples = as.list(colnames(raw_counts))
raw_counts_genes = as.list(rownames(raw_counts))
gene_annotations_entries = as.list(gene_annotations$ENSEMBL)
sample_annotations_entries = as.list(rownames(sample_annotations))

length(setdiff(raw_counts_samples, sample_annotations_entries)) #0
length(setdiff(sample_annotations_entries, raw_counts_samples)) #0
length(setdiff(raw_counts_genes, gene_annotations_entries)) #32489
length(setdiff(gene_annotations_entries, raw_counts_genes)) #0

genes_missing_annotation = c(raw_counts_genes[!(raw_counts_genes %in% gene_annotations_entries)])
gene_annotations_supplemented = data.frame("SYMBOL" = vector(mode="character", length=length(genes_missing_annotation)),
                                           "NAME" = vector(mode="character", length=length(genes_missing_annotation))
                                           )
gene_annotations_supplemented$ENSEMBL = genes_missing_annotation
gene_annotations_supplemented$ENTREZID = NA
gene_annotations_supplemented = gene_annotations_supplemented[,c(3,4,1,2)]
gene_annotations_supplemented = rbind(gene_annotations, gene_annotations_supplemented)

# Filter the count data for lowly-expressed genes, for example, only keep genes with a CPM >= 1 in at least 75% samples, in at least one of the groups.
# After consultation, normalizing by DESeq2's normalized counts >= 5 instead of CPM

library(DESeq2)
library("dplyr")

raw_counts_rowid = tibble::rownames_to_column(raw_counts, "ENSEMBL") #generate tale without header or rownames for DESeq2
sample_annotations_rowid = tibble::rownames_to_column(sample_annotations, "SAMPLE_ID") #generate tale without header or rownames for DESeq2
dds = DESeqDataSetFromMatrix(countData=raw_counts_rowid, 
                              colData=sample_annotations_rowid, 
                              design=~type, tidy = TRUE)
dds = estimateSizeFactors(dds)
normalized_counts = counts(dds, normalized=TRUE)

counts_threshold = 5 # assign threshold for normalized counts
normalized_counts_threshold = normalized_counts > counts_threshold # filter by counts threshold

sample_threshold = 0.75 # percent of samples above threshold in condition group
sample_threshold_normal = floor(sum(sample_annotations$type == "normal", na.rm = TRUE) * sample_threshold) # calculate threshold samples for normal
sample_threshold_lesional = floor(sum(sample_annotations$type == "lesional", na.rm = TRUE) * sample_threshold) # calculate threshold samples for lesional

normal_list = sample_annotations_rowid[which(sample_annotations_rowid$type == "normal"), ][1] #get indx of normal samples
lesional_list = sample_annotations_rowid[which(sample_annotations_rowid$type == "lesional"), ][1] #get indx of lesional samples

normalized_counts_threshold_normal = normalized_counts_threshold[, which(colnames(normalized_counts_threshold) %in% normal_list$SAMPLE_ID)] #generate matrix of normal samples
genes_passing_normal = rowSums(normalized_counts_threshold_normal) >= sample_threshold_normal # get boolean for genes passing sample threshold in normal
sum(genes_passing_normal, na.rm = TRUE) #23977

normalized_counts_threshold_lesional = normalized_counts_threshold[, which(colnames(normalized_counts_threshold) %in% lesional_list$SAMPLE_ID)] #generate matrix of lesional samples
genes_passing_lesional = rowSums(normalized_counts_threshold_lesional) >= sample_threshold_lesional # get boolean for genes passing sample threshold in lesional
sum(genes_passing_lesional, na.rm = TRUE) #23346

normalized_counts_filtered = normalized_counts[which(genes_passing_normal==TRUE | genes_passing_lesional==TRUE), ] # filter normalized counts by genes passing sample threshold in normal or lesional



