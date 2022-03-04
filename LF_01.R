# Champions Oncology - Exercise for bioinformatics candidates
# Submitted by Lior Faigon

#Load required libraries:
#library(Biobase) # for ExpressionSet


# Load the data into R and make sure the count and annotation data are consistent with each other.

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


