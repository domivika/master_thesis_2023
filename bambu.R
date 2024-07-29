
# data preparation
library(bambu)
library(ggbio)
fa.file <- "/path/to/Homo_sapiens_assembly38_masked.fasta"
gtf.file <-  "/path/to/gencode.v46.basic.annotation.gtf"
annotations <- prepareAnnotations(gtf.file) 
samples.bam <- list.files("/path/to/bams", pattern = ".primary.bam$", full.names = TRUE)

# running Bambu 
se <- bambu(reads = samples.bam, annotations = annotations, genome = fa.file, ncore = 17)

# save output
writeBambuOutput(se, path = "./output")

# plot
plotBambu(se, type = "annotation", transcript_id = "ENSG00000175061.20")
plotBambu(se, type = "heatmap") # heatmap 
plotBambu(se, type = "pca") # PCA

# Load necessary libraries
library(DESeq2)
library(pheatmap)

# Set working directory to where your files are located
setwd("/path/to/your/directory")

# Read the counts data
counts <- read.table("counts_gene.txt", header=TRUE, row.names=1, sep="\t")

# Read the metadata
sample_metadata <- data.frame(row.names = samples)

# Ensure the counts are integers
counts <- round(counts)


# Create DESeq2 dataset with the condition in the design formula
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_metadata,
                              design = ~ 1)  # Adjust as per your experimental design

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get the DESeq2 results
results <- results(dds)

# Remove rows with NA in pvalue
results <- results[complete.cases(results$pvalue), ]

# Filter out genes with not significant p-value (e.g., p-value >= 0.05)
significant_genes <- results[results$pvalue < 0.05, ]

# Perform variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Get the variance-stabilized data
vsd_data <- assay(vsd)

# Filter vsd_data for significant genes
vsd_data_filtered <- vsd_data[rownames(significant_genes), ]

# Calculate row variances
row_vars <- rowVars(vsd_data_filtered)

# Sort genes by variance
sorted_genes <- rownames(vsd_data_filtered)[order(row_vars, decreasing=TRUE)]

# Get the top variable genes
top_variable_genes <- sorted_genes[1:25]  # Adjust the number as needed

# Plot the top variable genes
pheatmap(vsd_data_filtered[top_variable_genes, ])