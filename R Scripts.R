Count Normalization using DESeq2
# Load libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# Loading Files

# 1 hour Control replicate 1
mock1hA <- read_tsv("mock1hA_ReadsPerGene.tsv",
                           col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                           skip = 4) # skip the first 4 lines

# 1 hour Control replicate 2
mock1hB <- read_tsv("mock1hB_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines
# 6 hour Control replicate 1
mock6hA <- read_tsv("mock6hA_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 6 hour Control replicate 2
mock6hB <- read_tsv("mock6hB_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 12 hour Control replicate 1
mock12hA <- read_tsv("mock12hA_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 12 hour Control replicate 2
mock12hB <- read_tsv("mock12hB_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 1 hour Treatment replicate 1
vir1hA <- read_tsv("vir1hA_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 1 hour Treatment replicate 2
vir1hB <- read_tsv("vir1hB_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines
# 6 hour Treatment replicate 1
vir6hA <- read_tsv("vir6hA_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 6 hour Treatment replicate 2
vir6hB <- read_tsv("vir6hB_ReadsPerGene.tsv",
                    col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                    skip = 4) # skip the first 4 lines

# 12 hour Treatment replicate 1
vir12hA <- read_tsv("vir12hA_ReadsPerGene.tsv",
                     col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                     skip = 4) # skip the first 4 lines

# 12 hour Treatment replicate 2
vir12hB <- read_tsv("vir12hB_ReadsPerGene.tsv",
                     col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                     skip = 4) # skip the first 4 lines

# Compiling Data into Dataframe
counts_1hr <- data.frame(row.names = mock1hA$gene_id,
                  mock1hA = mock1hA$sense,
                  mock1hB = mock1hB$sense,
                  vir1hA = vir1hA$sense,
                  vir1hB = vir1hB$sense)

counts_6hr <- data.frame(row.names = mock6hA$gene_id,
                         mock6hA = mock6hA$sense,
                         mock6hB = mock6hB$sense,
                         vir6hA = vir6hA$sense,
                         vir6hB = vir6hB$sense)

counts_12hr <- data.frame(row.names = mock12hA$gene_id,
                         mock12hA = mock12hA$sense,
                         mock12hB = mock12hB$sense,
                         vir12hA = vir12hA$sense,
                         vir12hB = vir12hB$sense)

# Transforming counts into a Matrix
counts_1hr_matrix<- as.matrix(counts_1hr)
class(counts_1hr_matrix)

counts_6hr_matrix<- as.matrix(counts_6hr)
class(counts_6hr_matrix)

counts_12hr_matrix<- as.matrix(counts_12hr)
class(counts_12hr_matrix)

# Making Metadata file that contains column information for matrix
metadata_1hr <- data.frame(row.names = colnames(counts_1hr_matrix),
                       condition = c("1_hour_c", "1_hour_c", "1_hour_t", "1_hour_t"))

metadata_6hr <- data.frame(row.names = colnames(counts_6hr_matrix),
                           condition = c("6_hour_c", "6_hour_c", "6_hour_t", "6_hour_t"))

metadata_12hr <- data.frame(row.names = colnames(counts_12hr_matrix),
                           condition = c("12_hour_c", "12_hour_c", "12_hour_t", "12_hour_t"))

# Confirming Metadata row names match Matrix Column Names
colnames(counts_1hr_matrix) == rownames(metadata_1hr)
colnames(counts_6hr_matrix) == rownames(metadata_6hr)
colnames(counts_12hr_matrix) == rownames(metadata_12hr)

# Creating DESeq2 object
dds_1hr_matrix <- DESeqDataSetFromMatrix(countData = counts_1hr_matrix, #matrix
                                     colData = metadata_1hr, #metadata file
                                     design = ~condition)

dds_6hr_matrix <- DESeqDataSetFromMatrix(countData = counts_6hr_matrix, #matrix
                                         colData = metadata_6hr, #metadata file
                                         design = ~condition)

dds_12hr_matrix <- DESeqDataSetFromMatrix(countData = counts_12hr_matrix, #matrix
                                         colData = metadata_12hr, #metadata file
                                         design = ~condition)

## Set control condition using the relevel function
dds_1hr_matrix$condition <- relevel(dds_1hr_matrix$condition, ref = "1_hour_c")
dds_6hr_matrix$condition <- relevel(dds_6hr_matrix$condition, ref = "6_hour_c")
dds_12hr_matrix$condition <- relevel(dds_12hr_matrix$condition, ref = "12_hour_c")

# Check the levels of dds_matrix$condition
levels(dds_1hr_matrix$condition)
levels(dds_6hr_matrix$condition)
levels(dds_12hr_matrix$condition)

# Running DESeq2 on the datasets
dds_1hr <- DESeq(dds_1hr_matrix)
dds_6hr <- DESeq(dds_6hr_matrix)
dds_12hr <- DESeq(dds_12hr_matrix)

# Sanity Check by Sample Clustering
# Perform log transformation on our count data
rld_1hr <- rlog(dds_1hr)
rld_6hr <- rlog(dds_6hr)
rld_12hr <- rlog(dds_12hr)

# Generate a PCA plot with DESeq2's plotPCA function
plotPCA(rld_1hr, intgroup = "condition")
plotPCA(rld_6hr, intgroup = "condition")
plotPCA(rld_12hr, intgroup = "condition")

# Generate Heat Map
sample_dists_1hr <- dist(t(assay(rld_1hr)))
sample_1hr_dist_matrix <- as.matrix(sample_dists_1hr)
colnames(sample_1hr_dist_matrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_1hr_dist_matrix,
         clustering_distance_rows = sample_dists_1hr,
         clustering_distance_cols = sample_dists_1hr,
         col = colours)

sample_dists_6hr <- dist(t(assay(rld_6hr)))
sample_6hr_dist_matrix <- as.matrix(sample_dists_6hr)
colnames(sample_6hr_dist_matrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_6hr_dist_matrix,
         clustering_distance_rows = sample_dists_6hr,
         clustering_distance_cols = sample_dists_6hr,
         col = colours)

sample_dists_12hr <- dist(t(assay(rld_12hr)))
sample_12hr_dist_matrix <- as.matrix(sample_dists_12hr)
colnames(sample_12hr_dist_matrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_12hr_dist_matrix,
         clustering_distance_rows = sample_dists_12hr,
         clustering_distance_cols = sample_dists_12hr,
         col = colours)

# Saving DESeqDataSet object
saveRDS(dds_1hr, "dds_1hr.rds")
saveRDS(dds_6hr, "dds_6hr.rds")
saveRDS(dds_12hr, "dds_12hr.rds")

# Names of the results that DESeq2 calculated
resultsNames(dds_1hr)
resultsNames(dds_6hr)
resultsNames(dds_12hr)

# Extracting results for comparison between the Treatment and Control groups
res_1hr <- results(dds_1hr, name = "condition_1_hour_t_vs_1_hour_c") %>% as.data.frame()
head(res_1hr)

res_6hr <- results(dds_6hr, name = "condition_6_hour_t_vs_6_hour_c") %>% as.data.frame()
head(res_6hr)

res_12hr <- results(dds_12hr, name = "condition_12_hour_t_vs_12_hour_c") %>% as.data.frame()
head(res_12hr)

# Filtering out NA
res_1hr_no_NA <- res_1hr %>%
  drop_na()

res_6hr_no_NA <- res_6hr %>%
  drop_na()

res_12hr_no_NA <- res_12hr %>%
  drop_na()

# Filtering results with adjusted p-value < 0.05
res_1hr_filtered <- res_1hr_no_NA %>%
  filter(padj <= 0.05) %>%
  rownames_to_column("gene_id")

res_6hr_filtered <- res_6hr_no_NA %>%
  filter(padj <= 0.05) %>%
  rownames_to_column("gene_id")

res_12hr_filtered <- res_12hr_no_NA %>%
  filter(padj <= 0.05) %>%
  rownames_to_column("gene_id")

# Filtering biologically relevant for 2X higher/lower expression
res_1hr_filtered_final <- res_1hr_filtered %>%
  filter(log2FoldChange <= 0 | log2FoldChange >= 0) %>%
  rownames_to_column("gene_id")

res_6hr_filtered_final <- res_6hr_filtered %>%
  filter(log2FoldChange <= 0 | log2FoldChange >= 0) %>%
  rownames_to_column("gene_id")

res_12hr_filtered_final <- res_12hr_filtered %>%
  filter(log2FoldChange <= 0 | log2FoldChange >= 0) %>%
  rownames_to_column("gene_id")

# Filtering for the top 10 genes
top10_1hr_genes <- res_1hr_filtered_final %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 10)

top10_6hr_genes <- res_6hr_filtered_final %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 10)

top10_12hr_genes <- res_12hr_filtered_final %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 10)

# Filtering for the bottom 10 genes
bottom10_1hr_genes <- res_1hr_filtered_final %>%
  arrange(log2FoldChange) %>%
  head(n = 10)

bottom10_6hr_genes <- res_6hr_filtered_final %>%
  arrange(log2FoldChange) %>%
  head(n = 10)

bottom10_12hr_genes <- res_12hr_filtered_final %>%
  arrange(log2FoldChange) %>%
  head(n = 10)

# Saving biologically relevant .csv
write_csv(res_1hr_filtered, "arabidopsis_1hr_results_no_log.csv")
write_csv(res_6hr_filtered, "arabidopsis_6hr_results_no_log.csv")
write_csv(res_12hr_filtered, "arabidopsis_12hr_results_no_log.csv")
write_csv(res_1hr_filtered_final, "arabidopsis_1hr_results.csv")
write_csv(res_6hr_filtered_final, "arabidopsis_6hr_results.csv")
write_csv(res_12hr_filtered_final, "arabidopsis_12hr_results.csv")

# Saving top 10 genes .csv
write_csv(top10_1hr_genes, "arabidopsis_1hr_filtered_up.csv")
write_csv(top10_6hr_genes, "arabidopsis_6hr_filtered_up.csv")
write_csv(top10_12hr_genes, "arabidopsis_12hr_filtered_up.csv")

# Saving bottom 10 genes .csv
write_csv(bottom10_1hr_genes, "arabidopsis_1hr_filtered_down.csv")
write_csv(bottom10_6hr_genes, "arabidopsis_6hr_filtered_down.csv")
write_csv(bottom10_12hr_genes, "arabidopsis_12hr_filtered_down.csv")

GO Enrichment Analysis

# Running topGo

# Installing topGo
BiocManager::install("topGo")

# Load topGo
library("topGO")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(topGO))

# Read in "gene universe" file
geneID2GO <- readMappings("arabidopsis_GOIDs.tsv")
geneUniverse <- names(geneID2GO)

# Load the differential expression data for 1hr
arabidopsis_1hr <- read.csv("arabidopsis_1hr_results.csv")

# Filter for statistically significant upregulated genes
arabidopsis_1hr_up_genes <- arabidopsis_1hr %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)

# Filter for statistically signficant downregulated genes
arabidopsis_1hr_down_genes <- arabidopsis_1hr %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)

# process files
arabidopsis_1hr_upregulated_genes <- as.character(arabidopsis_1hr_up_genes$gene_id)
arabidopsis_1hr_downregulated_genes <- as.character(arabidopsis_1hr_down_genes$gene_id)

# factor the names
arabidopsis_1hr_up_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_1hr_upregulated_genes))
arabidopsis_1hr_down_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_1hr_downregulated_genes))

# SEt names for the gene list
names(arabidopsis_1hr_up_gene_list) <- geneUniverse
names(arabidopsis_1hr_down_gene_list) <- geneUniverse

# build the GOdata object in topGO for upregulated
arabidopsis_1hr_up_GO_data <- new("topGOdata",
                  description = "Arabidopsis_1hr",
                  ontology = "BP",
                  allGenes = arabidopsis_1hr_up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

# build the GOdata object in topGO for downregulated
arabidopsis_1hr_down_GO_data <- new("topGOdata",
                    description = "Arabidopsis_1hr",
                    ontology = "BP",
                    allGenes = arabidopsis_1hr_down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# run the Fisher's exact tests with the weight01 algorithm
# resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
arabidopsis_1hr_up_result <- runTest(arabidopsis_1hr_up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")
arabidopsis_1hr_down_result <- runTest(arabidopsis_1hr_down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# summarize the results
arabidopsis_1hr_up_summary <- GenTable(arabidopsis_1hr_up_GO_data,
                       weight01 = arabidopsis_1hr_up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 20)
arabidopsis_1hr_down_summary <- GenTable(arabidopsis_1hr_down_GO_data,
                         weight01 = arabidopsis_1hr_down_result,
                         orderBy = "down_result",
                         ranksOf = "down_result",
                         topNodes = 20)



# Load the differential expression data for 6hr
arabidopsis_6hr <- read.csv("arabidopsis_6hr_results.csv")

# Filter for statistically significant upregulated genes
arabidopsis_6hr_up_genes <- arabidopsis_6hr %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)

# Filter for statistically signficant downregulated genes
arabidopsis_6hr_down_genes <- arabidopsis_6hr %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)

# process files
arabidopsis_6hr_upregulated_genes <- as.character(arabidopsis_6hr_up_genes$gene_id)
arabidopsis_6hr_downregulated_genes <- as.character(arabidopsis_6hr_down_genes$gene_id)

# factor the names
arabidopsis_6hr_up_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_6hr_upregulated_genes))
arabidopsis_6hr_down_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_6hr_downregulated_genes))

# SEt names for the gene list
names(arabidopsis_6hr_up_gene_list) <- geneUniverse
names(arabidopsis_6hr_down_gene_list) <- geneUniverse

# build the GOdata object in topGO for upregulated
arabidopsis_6hr_up_GO_data <- new("topGOdata",
                                  description = "Arabidopsis_6hr",
                                  ontology = "BP",
                                  allGenes = arabidopsis_6hr_up_gene_list,
                                  annot = annFUN.gene2GO,
                                  gene2GO = geneID2GO)

# build the GOdata object in topGO for downregulated
arabidopsis_6hr_down_GO_data <- new("topGOdata",
                                    description = "Arabidopsis_6hr",
                                    ontology = "BP",
                                    allGenes = arabidopsis_6hr_down_gene_list,
                                    annot = annFUN.gene2GO,
                                    gene2GO = geneID2GO)

# run the Fisher's exact tests with the weight01 algorithm
# resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
arabidopsis_6hr_up_result <- runTest(arabidopsis_6hr_up_GO_data,
                                     algorithm = "weight01",
                                     statistic = "fisher")
arabidopsis_6hr_down_result <- runTest(arabidopsis_6hr_down_GO_data,
                                       algorithm = "weight01",
                                       statistic = "fisher")

# summarize the results
arabidopsis_6hr_up_summary <- GenTable(arabidopsis_6hr_up_GO_data,
                                       weight01 = arabidopsis_6hr_up_result,
                                       orderBy = "up_result",
                                       ranksOf = "up_result",
                                       topNodes = 20)
arabidopsis_6hr_down_summary <- GenTable(arabidopsis_6hr_down_GO_data,
                                         weight01 = arabidopsis_6hr_down_result,
                                         orderBy = "down_result",
                                         ranksOf = "down_result",
                                         topNodes = 20)



# Load the differential expression data for 12hr
arabidopsis_12hr <- read.csv("arabidopsis_12hr_results.csv")

# Filter for statistically significant upregulated genes
arabidopsis_12hr_up_genes <- arabidopsis_12hr %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)

# Filter for statistically signficant downregulated genes
arabidopsis_12hr_down_genes <- arabidopsis_12hr %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)

# process files
arabidopsis_12hr_upregulated_genes <- as.character(arabidopsis_12hr_up_genes$gene_id)
arabidopsis_12hr_downregulated_genes <- as.character(arabidopsis_12hr_down_genes$gene_id)

# factor the names
arabidopsis_12hr_up_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_12hr_upregulated_genes))
arabidopsis_12hr_down_gene_list <- factor(as.integer(geneUniverse %in% arabidopsis_12hr_downregulated_genes))

# Set names for the gene list
names(arabidopsis_12hr_up_gene_list) <- geneUniverse
names(arabidopsis_12hr_down_gene_list) <- geneUniverse

# build the GOdata object in topGO for upregulated
arabidopsis_12hr_up_GO_data <- new("topGOdata",
                                  description = "Arabidopsis_12hr",
                                  ontology = "BP",
                                  allGenes = arabidopsis_12hr_up_gene_list,
                                  annot = annFUN.gene2GO,
                                  gene2GO = geneID2GO)

# build the GOdata object in topGO for downregulated
arabidopsis_12hr_down_GO_data <- new("topGOdata",
                                    description = "Arabidopsis_12hr",
                                    ontology = "BP",
                                    allGenes = arabidopsis_12hr_down_gene_list,
                                    annot = annFUN.gene2GO,
                                    gene2GO = geneID2GO)

# run the Fisher's exact tests with the weight01 algorithm
# resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
arabidopsis_12hr_up_result <- runTest(arabidopsis_12hr_up_GO_data,
                                     algorithm = "weight01",
                                     statistic = "fisher")
arabidopsis_12hr_down_result <- runTest(arabidopsis_12hr_down_GO_data,
                                       algorithm = "weight01",
                                       statistic = "fisher")

# summarize the results
arabidopsis_12hr_up_summary <- GenTable(arabidopsis_12hr_up_GO_data,
                                       weight01 = arabidopsis_12hr_up_result,
                                       orderBy = "up_result",
                                       ranksOf = "up_result",
                                       topNodes = 20)
arabidopsis_12hr_down_summary <- GenTable(arabidopsis_12hr_down_GO_data,
                                         weight01 = arabidopsis_12hr_down_result,
                                         orderBy = "down_result",
                                         ranksOf = "down_result",
                                         topNodes = 20)

