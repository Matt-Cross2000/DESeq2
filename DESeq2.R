

## Package Installation and library loading

install.packages(c("readr", "jsonlite", "ggpubr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.18")
BiocManager::install(c("tximport", "DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db"))

library(readr)
library(jsonlite)
library(ggpubr)
library(tximport)
library(DESeq2)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define file paths
base_path <- "C:/Users/mattcross/Documents/PGT Precision Medicine/High Throughput Methods/Practical"
quant_path <- file.path(base_path, "Quant Files")

# Import transcript-to-gene mapping
tx2gene <- read_csv(file.path(base_path, "tx2gene.gencode.v44.csv"))

# List Salmon quant files
sample_files <- list.files(quant_path, full.names = TRUE)

# Import expression data
txi <- tximport(sample_files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# Define sample metadata
sample_conditions <- c("Derived", "Derived", "Derived", "Parental", "Parental", "Parental")
sample_table <- data.frame(
  run = basename(sample_files),
  condition = factor(sample_conditions)
)
rownames(sample_table) <- colnames(txi$counts)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds)
res_lfc <- lfcShrink(dds, coef = "condition_Parental_vs_Derived", type = "apeglm")

# Filter results (optional)
res_ordered <- res[order(res$pvalue), ]
res_sig <- results(dds, alpha = 0.05)

# Summaries
summary(res)
summary(res_sig)

# MA Plots
par(mfrow = c(1, 3))
plotMA(res, ylim = c(-2, 2), main = "Raw")
plotMA(res_lfc, ylim = c(-2, 2), main = "LFC Shrinkage")
plotMA(res_sig, ylim = c(-2, 2), main = "Adjusted p < 0.05")

# Prepare gene list for enrichment (sorted named vector)
geneList <- res$log2FoldChange
names(geneList) <- sub("\\.[0-9]*$", "", rownames(res))
geneList <- sort(geneList, decreasing = TRUE)

# Bitr conversion
gene.df <- bitr(names(geneList), fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

# Group GO
ggo <- groupGO(
  gene = gene.df$ENSEMBL,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  level = 3,
  readable = TRUE
)

# GSEA
ego_gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  minGSSize = 100,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Enrichment Analysis
ego_enrich <- enrichGO(
  gene = names(geneList),
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot results
gse_plot <- dotplot(ego_gsea, showCategory = 50)
enrich_plot <- dotplot(ego_enrich, showCategory = 50)

# Combine plots
ggarrange(gse_plot, enrich_plot)




