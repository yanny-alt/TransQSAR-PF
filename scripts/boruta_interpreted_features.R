# BIOLOGICAL INTERPRETATION OF BORUTA-SELECTED TRANSCRIPTOMIC FEATURES
# This script maps the 13 Boruta features to actual genes and their biological functions

library(data.table)
library(dplyr)


setwd("/Users/favourigwezeke/Personal_System/Research/Dr. Charles Nnadi/Project 1")


cat("=== BIOLOGICAL INTERPRETATION OF BORUTA FEATURES ===\n")

# Your 13 Boruta-selected transcriptomic features with importance scores
boruta_features <- data.frame(
  Feature = c("CQ_DE.CQ_DE_128", "CQ_DE.CQ_DE_169", 
              "Genotype_DE.Genotype_DE_88", "Genotype_DE.Genotype_DE_92", 
              "Genotype_DE.Genotype_DE_98", "Genotype_DE.Genotype_DE_184",
              "CQ_106_1_DE.CQ_106_1_DE_40", "CQ_106_1_DE.CQ_106_1_DE_79", 
              "CQ_106_1_DE.CQ_106_1_DE_135",
              "Variability.Variability_Pf.12.198.0_CDS_at", 
              "Variability.Variability_Pf.2.13.0_CDS_at",
              "Variability.Variability_Pf.13_1.443.0_CDS_at", 
              "Variability.Variability_Pf.5.281.0_CDS_at"),
  Importance = c(5.308, 4.542, 4.172, 2.327, 1.969, 1.973, 
                 9.428, 3.288, 0.998, 4.612, 3.226, 4.106, 1.829),
  Feature_Type = rep("Transcriptomic_Boruta", 13)
)

# Sort by importance
boruta_features <- boruta_features[order(boruta_features$Importance, decreasing = TRUE), ]

cat("13 Boruta-Selected Transcriptomic Features (by importance):\n")
print(boruta_features)


cat("\n=== STEP 1: EXTRACT GENE IDs FROM FEATURE NAMES ===\n")

# Extract the row indices to map back to original data
extract_gene_info <- function(feature_name) {
  if (grepl("Variability", feature_name)) {
    # Extract probe ID from Variability features
    probe_id <- gsub("Variability\\.Variability_", "", feature_name)
    probe_id <- gsub("_CDS_at", "", probe_id)
    return(list(type = "Variability", probe_id = probe_id, row_index = NA))
  } else if (grepl("CQ_DE\\.", feature_name)) {
    # Extract row index from differential expression features
    row_index <- as.numeric(gsub("CQ_DE\\.CQ_DE_", "", feature_name))
    return(list(type = "CQ_DE", probe_id = NA, row_index = row_index))
  } else if (grepl("Genotype_DE\\.", feature_name)) {
    row_index <- as.numeric(gsub("Genotype_DE\\.Genotype_DE_", "", feature_name))
    return(list(type = "Genotype_DE", probe_id = NA, row_index = row_index))
  } else if (grepl("CQ_106_1_DE\\.", feature_name)) {
    row_index <- as.numeric(gsub("CQ_106_1_DE\\.CQ_106_1_DE_", "", feature_name))
    return(list(type = "CQ_106_1_DE", probe_id = NA, row_index = row_index))
  }
  return(list(type = "Unknown", probe_id = NA, row_index = NA))
}

# Apply extraction to all features
boruta_features$gene_info <- lapply(boruta_features$Feature, extract_gene_info)
boruta_features$feature_type <- sapply(boruta_features$gene_info, function(x) x$type)
boruta_features$row_index <- sapply(boruta_features$gene_info, function(x) x$row_index)
boruta_features$probe_id <- sapply(boruta_features$gene_info, function(x) x$probe_id)

cat("Feature breakdown by type:\n")
table(boruta_features$feature_type)


cat("\n=== STEP 2: LOAD YOUR ORIGINAL DATA TO MAP GENES ===\n")

# Load the differential expression results to map row indices to genes
if (file.exists("CQ_differential_expression_all.csv")) {
  de_cq_all <- fread("CQ_differential_expression_all.csv")
  cat("Loaded CQ DE data with", nrow(de_cq_all), "genes\n")
} else {
  cat("CQ_differential_expression_all.csv not found - will need to create mapping manually\n")
}

if (file.exists("Genotype_differential_expression_all.csv")) {
  de_genotype_all <- fread("Genotype_differential_expression_all.csv")
  cat("Loaded Genotype DE data with", nrow(de_genotype_all), "genes\n")
} else {
  cat("Will attempt to use existing de_genotype object\n")
}

# Load GPL annotation file for probe-to-gene mapping
if (file.exists("GPL1321-15512.txt")) {
  lines <- readLines("GPL1321-15512.txt")
  header_line <- which(!grepl("^#", lines))[1]
  gpl <- fread("GPL1321-15512.txt", skip = header_line - 1, sep = "\t", header = TRUE, fill = TRUE)
  
  # Clean up the mapping
  probe_to_gene <- gpl[, .(ID, ORF, `Gene Symbol`, `Gene Title`)]
  probe_to_gene <- probe_to_gene[ORF != "" & !is.na(ORF)]
  colnames(probe_to_gene) <- c("probe_id", "gene_id", "gene_symbol", "gene_title")
  
  cat("Loaded GPL annotation with", nrow(probe_to_gene), "probe-gene mappings\n")
} else {
  cat("GPL annotation file not found - will use available data\n")
  probe_to_gene <- NULL
}


cat("\n=== STEP 3: MAP BORUTA FEATURES TO ACTUAL GENES ===\n")

# Initialize results dataframe
boruta_gene_mapping <- data.frame(
  Feature = boruta_features$Feature,
  Importance = boruta_features$Importance,
  Feature_Type = boruta_features$feature_type,
  Gene_ID = NA,
  Gene_Symbol = NA,
  Gene_Description = NA,
  Biological_Context = NA,
  stringsAsFactors = FALSE
)

# Map each feature type
for (i in 1:nrow(boruta_features)) {
  feature_info <- boruta_features$gene_info[[i]]
  
  if (feature_info$type == "Variability" && !is.null(probe_to_gene)) {
    # Map variability features using probe IDs
    probe_match <- probe_to_gene[probe_to_gene$probe_id == feature_info$probe_id, ]
    if (nrow(probe_match) > 0) {
      boruta_gene_mapping$Gene_ID[i] <- probe_match$gene_id[1]
      boruta_gene_mapping$Gene_Symbol[i] <- probe_match$gene_symbol[1]
      boruta_gene_mapping$Gene_Description[i] <- probe_match$gene_title[1]
      boruta_gene_mapping$Biological_Context[i] <- "Expression variability across strains"
    }
  } else if (feature_info$type %in% c("CQ_DE", "Genotype_DE", "CQ_106_1_DE")) {
    # For DE features, we need to map row indices back to the original data
    if (feature_info$type == "CQ_DE" && exists("de_cq_all")) {
      if (feature_info$row_index <= nrow(de_cq_all)) {
        gene_row <- de_cq_all[feature_info$row_index, ]
        if (!is.null(probe_to_gene) && "V1" %in% names(gene_row)) {
          probe_match <- probe_to_gene[probe_to_gene$probe_id == gene_row$V1, ]
          if (nrow(probe_match) > 0) {
            boruta_gene_mapping$Gene_ID[i] <- probe_match$gene_id[1]
            boruta_gene_mapping$Gene_Symbol[i] <- probe_match$gene_symbol[1]
            boruta_gene_mapping$Gene_Description[i] <- probe_match$gene_title[1]
            boruta_gene_mapping$Biological_Context[i] <- paste("Chloroquine response,", 
                                                               "logFC =", round(gene_row$logFC, 3))
          }
        }
      }
    } else if (feature_info$type == "Genotype_DE" && exists("de_genotype_all")) {
      if (feature_info$row_index <= nrow(de_genotype_all)) {
        gene_row <- de_genotype_all[feature_info$row_index, ]
        if (!is.null(probe_to_gene) && "V1" %in% names(gene_row)) {
          probe_match <- probe_to_gene[probe_to_gene$probe_id == gene_row$V1, ]
          if (nrow(probe_match) > 0) {
            boruta_gene_mapping$Gene_ID[i] <- probe_match$gene_id[1]
            boruta_gene_mapping$Gene_Symbol[i] <- probe_match$gene_symbol[1]
            boruta_gene_mapping$Gene_Description[i] <- probe_match$gene_title[1]
            boruta_gene_mapping$Biological_Context[i] <- paste("Genotype difference,", 
                                                               "logFC =", round(gene_row$logFC, 3))
          }
        }
      }
    }
  }
}

# Fill in missing information with feature names
boruta_gene_mapping$Gene_ID[is.na(boruta_gene_mapping$Gene_ID)] <- "Unknown"
boruta_gene_mapping$Gene_Symbol[is.na(boruta_gene_mapping$Gene_Symbol)] <- "Unknown"
boruta_gene_mapping$Gene_Description[is.na(boruta_gene_mapping$Gene_Description)] <- "Unknown"
boruta_gene_mapping$Biological_Context[is.na(boruta_gene_mapping$Biological_Context)] <- "Unknown"

cat("Biological mapping of Boruta features:\n")
print(boruta_gene_mapping)

cat("\n=== STEP 4: CATEGORIZE BIOLOGICAL FUNCTIONS ===\n")

# Categorize the biological functions based on your pathway analysis
categorize_biological_function <- function(gene_desc, gene_id, context) {
  gene_desc_lower <- tolower(gene_desc)
  gene_id_lower <- tolower(gene_id)
  context_lower <- tolower(context)
  
  # Based on GSEA results, categorize into major functional groups
  if (grepl("pfemp|var|erythrocyte membrane", gene_desc_lower) || 
      grepl("pf3d7_04|pf3d7_06|pf3d7_12", gene_id_lower)) {
    return("Virulence & Antigenic Variation")
  } else if (grepl("ribosom|translation|rna|protein synthesis", gene_desc_lower)) {
    return("Protein Synthesis & RNA Processing")
  } else if (grepl("conserved|hypothetical|unknown", gene_desc_lower)) {
    return("Conserved Unknown Function")
  } else if (grepl("metabol|enzyme|kinase|phosphatase", gene_desc_lower)) {
    return("Metabolism & Signaling")
  } else if (grepl("transport|membrane|channel", gene_desc_lower)) {
    return("Transport & Membrane")
  } else if (grepl("dna|replication|repair|transcription", gene_desc_lower)) {
    return("DNA/RNA Processing")
  } else if (grepl("chloroquine|drug", context_lower)) {
    return("Drug Response")
  } else if (grepl("variability", context_lower)) {
    return("Strain Variability")
  } else {
    return("Other/Unclassified")
  }
}

# Apply categorization
boruta_gene_mapping$Functional_Category <- mapply(
  categorize_biological_function,
  boruta_gene_mapping$Gene_Description,
  boruta_gene_mapping$Gene_ID,
  boruta_gene_mapping$Biological_Context
)

# Summary by functional category
cat("Functional categorization of Boruta-selected features:\n")
category_summary <- boruta_gene_mapping %>%
  group_by(Functional_Category) %>%
  summarise(
    Count = n(),
    Total_Importance = round(sum(Importance), 2),
    Avg_Importance = round(mean(Importance), 2),
    Features = paste(Feature, collapse = "; ")
  )

print(category_summary)


cat("\n=== STEP 5: CREATE PUBLICATION-READY RESULTS ===\n")

# Create a clean table for publication
publication_table <- boruta_gene_mapping %>%
  arrange(desc(Importance)) %>%
  mutate(
    Rank = 1:n(),
    Importance_Rounded = round(Importance, 2)
  ) %>%
  select(Rank, Feature, Importance_Rounded, Functional_Category, 
         Gene_ID, Gene_Symbol, Gene_Description, Biological_Context)

# Save the results
write.csv(publication_table, "boruta_biological_interpretation_publication.csv", row.names = FALSE)
write.csv(category_summary, "boruta_functional_categories_summary.csv", row.names = FALSE)

cat("=== PUBLICATION-READY BIOLOGICAL INTERPRETATION ===\n")
cat("Top 5 most important Boruta-selected transcriptomic features:\n\n")

for (i in 1:min(5, nrow(publication_table))) {
  row <- publication_table[i, ]
  cat(sprintf("%d. %s (Importance: %.2f)\n", 
              row$Rank, row$Functional_Category, row$Importance_Rounded))
  cat(sprintf("   Gene: %s (%s)\n", row$Gene_Symbol, row$Gene_ID))
  cat(sprintf("   Function: %s\n", row$Gene_Description))
  cat(sprintf("   Context: %s\n\n", row$Biological_Context))
}

cat("\nKey Biological Insights from Boruta Selection:\n")
cat("1. Drug Response Signatures: Features capturing chloroquine response\n")
cat("2. Strain Variability: Genes with high expression variability across P. falciparum strains\n")
cat("3. Genotype Differences: Genetic background effects on gene expression\n")
cat("4. Functional Categories: Focus on", length(unique(boruta_gene_mapping$Functional_Category)), "major functional groups\n")

cat("\nFiles created for publication:\n")
cat("✓ boruta_biological_interpretation_publication.csv - Detailed gene mapping\n")
cat("✓ boruta_functional_categories_summary.csv - Functional category summary\n")


cat("\n=== STEP 6: CONNECT TO YOUR PATHWAY ANALYSIS ===\n")

# Load your previous GSEA results if available
if (file.exists("fgsea_results_clean.csv")) {
  fgsea_results <- fread("fgsea_results_clean.csv")
  cat("Connecting Boruta features to GSEA pathways...\n")
  
  # Show how Boruta features relate to your significant pathways
  if (nrow(fgsea_results[pval < 0.05]) > 0) {
    cat("Your significant GSEA pathways that may relate to Boruta features:\n")
    sig_pathways <- head(fgsea_results[pval < 0.05], 5)
    for (i in 1:nrow(sig_pathways)) {
      cat(sprintf("- %s (p = %.4f, NES = %.2f)\n", 
                  sig_pathways$pathway[i], sig_pathways$pval[i], sig_pathways$NES[i]))
    }
  }
}

cat("\n=== BIOLOGICAL INTERPRETATION COMPLETE ===\n")
cat("You now have the biological context for your 13 Boruta features!\n")
cat("This directly addresses the 'Which features and why' question for publication.\n")

# Create a summary statement for your paper
cat("\n=== SUMMARY FOR YOUR PAPER ===\n")
cat("Boruta feature selection identified 13 critical transcriptomic features from", 
    "an original set of 764 features (98.3% reduction). These features represent:\n")

for (i in 1:nrow(category_summary)) {
  cat(sprintf("- %s: %d features (%.1f%% total importance)\n",
              category_summary$Functional_Category[i],
              category_summary$Count[i],
              100 * category_summary$Total_Importance[i] / sum(boruta_features$Importance)))
}

cat("\nThis biologically-guided feature selection improved QSAR model performance\n")
cat("from R² = 0.719 (QSAR-only) to R² = 0.762 (Boruta-integrated), representing\n")
cat("a 6.1% improvement while dramatically reducing model complexity.\n")


library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# === Figure 1: Bar Plot of Feature Importance by Functional Category ===
importance_by_category <- publication_table %>%
  group_by(Functional_Category) %>%
  summarise(Total_Importance = sum(Importance_Rounded)) %>%
  arrange(desc(Total_Importance))

ggplot(importance_by_category, aes(x = reorder(Functional_Category, Total_Importance), 
                                   y = Total_Importance, 
                                   fill = Functional_Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Total Feature Importance by Functional Category",
       x = "Functional Category",
       y = "Sum of Feature Importance") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Save Figure 1
ggsave("figure1_feature_importance_by_category.png", width = 8, height = 5, dpi = 300)


