# ---- 0) Install libraries if needed ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma", "Biobase"))
install.packages("ggplot2")     
install.packages("ggfortify")   

install.packages("ggfortify")

# ---- Libraries ----
library(limma)
library(Biobase)
library(ggplot2)
library(ggfortify)

# ---- Set WD ----
setwd("/Users/favourigwezeke/Personal_System/Research/Dr. Charles Nnadi/Project 1")

# ---- 1) Load expression matrix (probes in rows, GSM columns) ----
# NOTE: first column has no header, contains probe IDs → use row.names = 1
exprs_matrix <- read.csv("Project1_expression_matrix.csv", 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = FALSE)
exprs_matrix <- as.matrix(exprs_matrix)

# ---- 2) Confirm sample order ----
expected_cols <- c(
  "GSM253220.CEL.gz","GSM253221.CEL.gz","GSM253222.CEL.gz","GSM253223.CEL.gz",
  "GSM253224.CEL.gz","GSM253225.CEL.gz","GSM253226.CEL.gz","GSM253227.CEL.gz",
  "GSM253228.CEL.gz","GSM253229.CEL.gz","GSM253230.CEL.gz","GSM253231.CEL.gz",
  "GSM253232.CEL.gz","GSM253233.CEL.gz","GSM253234.CEL.gz","GSM253235.CEL.gz",
  "GSM253236.CEL.gz","GSM253237.CEL.gz","GSM253238.CEL.gz","GSM253239.CEL.gz",
  "GSM253240.CEL.gz","GSM253241.CEL.gz","GSM253242.CEL.gz","GSM253243.CEL.gz"
)
if(!all(expected_cols %in% colnames(exprs_matrix))) {
  stop("ERROR: not all expected GSM columns were found in your expression matrix. Check column names exactly.")
}

# Reorder columns safely
exprs_matrix <- exprs_matrix[, expected_cols]

# ---- 3) Create sample metadata ----
treatment_full <- c(
  "Control","Control","Control","CQ","CQ","CQ",
  "Control","Control","Control","CQ","CQ","CQ",
  "Control","Control","Control","CQ","CQ","CQ",
  rep("gDNA", 6)
)
genotype_full <- c(
  rep("106_1", 6),
  rep("106_1_76I", 6),
  rep("106_1_76I_352K", 6),
  rep("gDNA", 6)
)
sample_df <- data.frame(
  Sample = colnames(exprs_matrix),
  Genotype = genotype_full,
  Treatment = treatment_full,
  stringsAsFactors = FALSE
)
cat("Full sample table (first 20 rows):\n")
print(sample_df)

# ---- 4) Exclude gDNA samples ----
rna_idx <- which(sample_df$Treatment != "gDNA")
exprs_rna <- exprs_matrix[, rna_idx]
sample_df_rna <- sample_df[rna_idx, ]
rownames(sample_df_rna) <- NULL
cat("After excluding gDNA, RNA samples dims:\n")
print(dim(exprs_rna))
print(sample_df_rna)


# ---- 5) Build grouping factor and design ----
group <- paste(sample_df_rna$Genotype, sample_df_rna$Treatment, sep = "_")
group_factor <- factor(group, levels = unique(group))
cat("Groups used:\n")
print(table(group_factor))
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
cat("Design matrix columns:\n")
print(colnames(design))
print(design)


# ---- 6) Quick QC: boxplot of samples and PCA ----
# Boxplot (sample distributions)
pdf("QC_boxplot_samples.pdf", width = 10, height = 6)
boxplot(as.data.frame(exprs_rna), las = 2, main = "Sample distributions (expression)", cex.axis=0.8)
dev.off()
cat("Saved QC_boxplot_samples.pdf\n")


pdf("QC_PCA_samples.pdf", width = 7, height = 6)
pca_res <- prcomp(t(exprs_rna), scale. = TRUE)
autoplot(pca_res, data = sample_df_rna, colour = 'Genotype', shape = 'Treatment', label = FALSE) +
  ggtitle("PCA of samples (RNA only)")
dev.off()
cat("Saved QC_PCA_samples.pdf\n")

# ---- 7) Fit limma model and contrasts ----

# Rename columns to valid R names
colnames(design) <- make.names(colnames(design))

# Check what they became
print(colnames(design))

# Define contrasts that match your experimental questions:
#  - CQ_effect: average CQ vs average Control across the three genotypes
#  - Genotype_effect:  (averaged Controls of mutant genotypes) vs WT Control
#  - Per-genotype CQ contrasts as well
contrast_matrix <- makeContrasts(
  CQ_effect = (X106_1_CQ + X106_1_76I_CQ + X106_1_76I_352K_CQ)/3 -
    (X106_1_Control + X106_1_76I_Control + X106_1_76I_352K_Control)/3,
  Genotype_effect = ((X106_1_76I_Control + X106_1_76I_352K_Control)/2) - X106_1_Control,
  CQ_106_1 = X106_1_CQ - X106_1_Control,
  CQ_76I   = X106_1_76I_CQ - X106_1_76I_Control,
  CQ_352K  = X106_1_76I_352K_CQ - X106_1_76I_352K_Control,
  levels = design
)

fit <- lmFit(exprs_rna, design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


# ---- 8) Extract DE results ----
de_CQ <- topTable(fit2, coef = "CQ_effect", number = Inf, sort.by = "P", adjust.method = "BH")
de_Genotype <- topTable(fit2, coef = "Genotype_effect", number = Inf, sort.by = "P", adjust.method = "BH")
de_CQ_106_1 <- topTable(fit2, coef = "CQ_106_1", number = Inf, sort.by = "P", adjust.method = "BH")

write.csv(de_CQ, "CQ_differential_expression_all.csv", row.names = TRUE)
write.csv(de_Genotype, "Genotype_differential_expression_all.csv", row.names = TRUE)
write.csv(de_CQ_106_1, "CQ_106_1_differential_expression_all.csv", row.names = TRUE)
cat("Saved differential expression CSVs.\n")

write.csv(head(de_CQ, 200), "CQ_top200.csv", row.names = TRUE)


# ---- 9) Top variable genes (signature) ----
gene_vars <- apply(exprs_rna, 1, var, na.rm = TRUE)
topN <- 978
top_variable_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(topN, length(gene_vars))]
signature_matrix <- exprs_rna[top_variable_genes, , drop = FALSE]
rownames(signature_matrix) <- top_variable_genes
write.csv(signature_matrix, "PF_expression_signatures_top978.csv", row.names = TRUE)
cat("Saved PF_expression_signatures_top978.csv\n")


# ---- 10) Group-average profiles ----
group_levels <- levels(group_factor)
group_profiles <- sapply(group_levels, function(g) {
  cols <- which(group_factor == g)
  rowMeans(signature_matrix[, cols, drop = FALSE])
})
colnames(group_profiles) <- group_levels
write.csv(group_profiles, "group_expression_profiles_top978.csv", row.names = TRUE)
cat("Saved group_expression_profiles_top978.csv\n")


# ---- 11) Save workspace ----
save.image("GSE10022_processed_workspace.RData")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("Saved workspace and session info.\n")

cat("\n✓ Pipeline complete. Key outputs:\n")
cat(" - QC_boxplot_samples.pdf\n - QC_PCA_samples.pdf\n - CQ_differential_expression_all.csv\n - Genotype_differential_expression_all.csv\n - PF_expression_signatures_top978.csv\n - group_expression_profiles_top978.csv\n - GSE10022_processed_workspace.RData\n")






