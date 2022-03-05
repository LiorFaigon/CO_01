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

# Generate an object that contains the library-size normalized log-CPM data. Save it as a binary file (.rda or .rds).
# The equivalent to library-size normalized log-CPM data using DESeq2 is the vst (variance stabilizing transformations) data.
# unclear if instructions refer to filtered or unfiltered data, doing both.

vst_normalized_counts = assay(vst(dds)) # vst for unfiltered data
vst_normalized_counts_filtered = vst_normalized_counts[which(rownames(vst_normalized_counts) %in% rownames(normalized_counts_filtered)), ] # filter vst data genes by the previously filtered data

save(vst_normalized_counts, file = "vst_normalized_counts_unfiltered.rda")
save(vst_normalized_counts_filtered, file = "vst_normalized_counts_filtered.rda")

#Generate basic plots of your choice to investigate its main properties (library sizes, densities, PCA coloured by group, etc...).
  # library sizes: 
    # test variance between normal and lesional samples: significantly unequal
    # t.test: p-value = 3.403e-09

library(ggplot2)
library(ggsignif)

lib_size_vst_norm_filt = data.frame("library_size" = colSums(vst_normalized_counts_filtered),
                                    "type" = sample_annotations$type)
var.test(lib_size_vst_norm_filt$library_size[which(lib_size_vst_norm_filt$type == "normal")], 
         lib_size_vst_norm_filt$library_size[which(lib_size_vst_norm_filt$type == "lesional")], 
         alternative = "two.sided")
lib_size_ttest = t.test(lib_size_vst_norm_filt$library_size ~ lib_size_vst_norm_filt$type)

ggplot(lib_size_vst_norm_filt, aes(type, library_size)) +
  geom_violin(scale = "area") +
  labs(title = "Library Sizes for Filtered Data with T.Test", x = "sample type", y = "library size (vst normalized)") +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  geom_signif(comparisons = list(c("normal", "lesional")),   
              map_signif_level=TRUE, annotation = c(signif(lib_size_ttest$p.value, digits=3)))

  # densities:

plot(density(vst_normalized_counts_filtered),main="Expression Distribution",xlab="VST normalized counts (filtered)", ylab="Density")
for (s in 1:length(colnames(vst_normalized_counts_filtered))){
  d <- density(vst_normalized_counts_filtered[,s])
  lines(d)
}

  # PCA:

library(ggbiplot)

PCA_01 = prcomp(t(vst_normalized_counts_filtered))
summary(PCA_01)

ggbiplot(PCA_01, 
         choices = c(1,2), 
         var.axes = F,
         obs.scale = 1, 
         var.scale = 1,
         groups = sample_annotations$type
         ) +
  ggtitle("Principal Component Analysis")

PC1genenames_pos = rownames(data.frame(sort(PCA_01$rotation[,"PC1"], decreasing=TRUE)[1:10])) # get 10 most positively influential genes for PC1
PC1genenames_neg = rownames(data.frame(sort(PCA_01$rotation[,"PC1"], decreasing=FALSE)[1:10])) # get 10 most negatively influential genes for PC1

# The PCA plot may suggest the presence of outlier/mis-labeled samples in this dataset. Try to identify them and remove them from the downstream analysis.
  #PCA suggests 1 mislabeled sample, no obvious outliers. checking for more mislabeled samples along suggested cluster coordinates

PC_dim_12 = data.frame(PCA_01$x[,1:2])
PC_dim_12$type = sample_annotations$type

PC_dim_12[which(PC_dim_12$PC1 > 0 & PC_dim_12$type == "normal"), ] # no mislabeled normal samples
PC_dim_12[which(PC_dim_12$PC1 < 0 & PC_dim_12$type == "lesional"), ] # found 1 mislabeled lesional sample: SRR1146216 
vst_normalized_counts_filtered = vst_normalized_counts_filtered[ , ! colnames(vst_normalized_counts_filtered) %in% c("SRR1146216")]
save(vst_normalized_counts_filtered, file = "vst_normalized_counts_filtered.rda") #update filtered file

#Run a differential expression analysis comparing lesional vs normal samples. This can be done according to your preference either on the count data or the normalized log-CPM data, using appropriate statistical method.

raw_counts_rowid_filtered = raw_counts_rowid[ raw_counts_rowid$ENSEMBL %in% rownames(vst_normalized_counts_filtered), 
                                              ! colnames(raw_counts_rowid) %in% c("SRR1146216")]
sample_annotations_rowid_filtered = sample_annotations_rowid[ sample_annotations_rowid$SAMPLE_ID %in% colnames(vst_normalized_counts_filtered), ]
dds_filt = DESeqDataSetFromMatrix(countData=raw_counts_rowid_filtered, 
                                  colData=sample_annotations_rowid_filtered, 
                                  design=~type, tidy = TRUE)
dds_filt = DESeq(dds_filt)

res = results(dds_filt)
summary(res)
plotMA(res) #most points on y=0 intercept and skewed lower as expected.
plotDispEsts(dds_filt) #data scattered around curve as expected.

#Export the results in a tab-separated text/CSV file: a table with genes in rows along with gene annotations and any relevant statistic.

res_tab = head(results(dds_filt, tidy=TRUE), n = length(rownames(vst_normalized_counts_filtered)))
colnames(res_tab)[1] = "ENSEMBL" #change colname from default
res_tab$"padj<0.01" = res_tab$padj<0.01 # add column for significantly DE genes

res_tab$"log2FC>|1|"[res_tab$padj<0.01 & res_tab$log2FoldChange > 1] = "UP" # add column for significantly DE genes passing logFC threshold upregulated
res_tab$"log2FC>|1|"[res_tab$padj<0.01 & res_tab$log2FoldChange < -1] = "DOWN" # add column for significantly DE genes passing logFC threshold downregulated

res_tab = merge(res_tab,gene_annotations_supplemented,by="ENSEMBL") # merge DE results with gene annotations to include genes with missing annotations
write.csv(res_tab,"differential_expression_res.csv", row.names = FALSE)

#Select the top 100 most significant *annotated* genes and generate a heatmap of the log-CPM data, with samples in columns, annotated with the group variable.
  #(DESeq2 VST normalized counts data)
library(genefilter)
library(RColorBrewer)

sig_genes = res_tab[(res_tab$`log2FC>|1|` == "UP" | res_tab$`log2FC>|1|` == "DOWN") & res_tab$NAME != "" , 1] #filter significantly DE above threshold with annotation
sig_genes = as.list(sig_genes[!is.na(sig_genes)]) #make list and remove NAs

sum(gene_annotations_supplemented[ gene_annotations_supplemented$ENSEMBL %in% sig_genes,]$NAME =="") # "Trust But Verify"

vst_normalized_counts_filtered_sig = vst_normalized_counts_filtered[rownames(vst_normalized_counts_filtered) %in% sig_genes, ]

topVarGenes = head( order( rowVars( vst_normalized_counts_filtered_sig ), decreasing=TRUE ), 100 ) #find 100 most significant annotated genes

sample_type_4hm = as.numeric(as.factor(sample_annotations[rownames(sample_annotations) %in% colnames(vst_normalized_counts_filtered_sig), 1])) #generate sample type numeric vector for heatmap

heatmap(vst_normalized_counts_filtered_sig[topVarGenes, ], 
        ColSideColors = brewer.pal(9, "Set1")[sample_type_4hm],
        scale="column")

# Generate a volcano plot (x-axis is the effect size and y-axis is the p-value) for this analysis. The selected 100 most significant genes should be colored.

most_sig = rownames(vst_normalized_counts_filtered_sig[topVarGenes,])

par(mfrow=c(1,1)) #reset par
with(res_tab, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot")) # Make the volcano plot
with(subset(res_tab, padj<.01 & abs(log2FoldChange)>1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue")) # overlay with significant genes
with(res_tab[ res_tab$ENSEMBL %in% most_sig, ], points(log2FoldChange, -log10(pvalue), pch=20, col="red")) #overlay with 100 most significant genes

# Export list of significantly upregulated genes for gene set enrichment analysis in DAVID

sig_up = res_tab[res_tab$`log2FC>|1|` == "UP", 1]
sig_up = sig_up[!is.na(sig_up)]
write.csv(sig_up,"significantly_upregulated_genes.csv", row.names = FALSE)
write.csv(as.character(gene_annotations_entries),"background_genes.csv", row.names = FALSE) #add background genes for Gene ontology enrichment analysis
