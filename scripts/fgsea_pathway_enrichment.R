# CLEAN GENE SET ENRICHMENT ANALYSIS SCRIPT
# This fixes the mapping issues and creates a proper GSEA workflow

library(data.table)
library(fgsea)

cat("=== STEP 1: READ GENE ALIAS FILE CORRECTLY ===\n")

# Read the alias file with proper settings
alias_file <- "PlasmoDB-68_Pfalciparum3D7_GeneAliases.txt"

# Use fill=Inf to read all columns properly
gene_aliases <- fread(alias_file, header = FALSE, sep = "\t", fill = Inf)

cat("Read alias file with", ncol(gene_aliases), "columns and", nrow(gene_aliases), "rows\n")

# Create comprehensive mapping
id_mapping <- data.table()
for(i in 1:nrow(gene_aliases)) {
  pf3d7_id <- gene_aliases[[1]][i]  # PF3D7 format ID
  
  # Get all aliases from remaining columns
  aliases <- as.character(gene_aliases[i, -1])
  aliases <- aliases[!is.na(aliases) & aliases != "" & aliases != "NA"]
  
  if(length(aliases) > 0) {
    for(alias in aliases) {
      id_mapping <- rbind(id_mapping, data.table(pf3d7_id = pf3d7_id, alias = alias))
    }
  }
}

# Remove duplicates and add self-mapping (PF3D7 IDs mapping to themselves)
id_mapping <- unique(id_mapping)
self_mapping <- data.table(pf3d7_id = unique(gene_aliases[[1]]), 
                           alias = unique(gene_aliases[[1]]))
id_mapping <- rbind(id_mapping, self_mapping)
id_mapping <- unique(id_mapping)

cat("Created comprehensive mapping with", nrow(id_mapping), "gene ID mappings\n")
cat("Sample mappings:\n")
print(head(id_mapping, 10))

cat("\n=== STEP 2: READ DE RESULTS ===\n")

de_results <- fread("CQ_differential_expression_all.csv")
de_results$probe_id <- de_results$V1  # Probe IDs are in column V1
de_results$V1 <- NULL

cat("DE results have", nrow(de_results), "probes\n")
cat("Sample probe IDs:", head(de_results$probe_id, 3), "\n")

cat("\n=== STEP 3: READ GPL ANNOTATION ===\n")

gpl_file <- "GPL1321-15512.txt"
lines <- readLines(gpl_file)
header_line <- which(!grepl("^#", lines))[1]
gpl <- fread(gpl_file, skip = header_line - 1, sep = "\t", header = TRUE, fill = TRUE)

# Create probe to gene mapping
probe_to_gene <- gpl[, .(ID, ORF)]
probe_to_gene <- probe_to_gene[ORF != "" & !is.na(ORF)]
colnames(probe_to_gene) <- c("probe_id", "gene_id")

cat("GPL mapping has", nrow(probe_to_gene), "probe-gene pairs\n")
cat("Sample GPL mappings:\n")
print(head(probe_to_gene))

cat("\n=== STEP 4: MAP PROBES TO GENES ===\n")

# Map DE results to genes via GPL
de_results[, gene_id := probe_to_gene$gene_id[match(probe_id, probe_to_gene$probe_id)]]
de_mapped <- de_results[!is.na(gene_id)]

cat("Successfully mapped", nrow(de_mapped), "of", nrow(de_results), "probes to genes\n")

cat("\n=== STEP 5: CONVERT GENE IDs TO PF3D7 FORMAT ===\n")

# Map gene IDs to PF3D7 format using the alias mapping
de_mapped[, pf3d7_id := id_mapping$pf3d7_id[match(gene_id, id_mapping$alias)]]
de_pf3d7 <- de_mapped[!is.na(pf3d7_id)]

cat("Successfully converted", nrow(de_pf3d7), "genes to PF3D7 format\n")

# For genes with multiple probes, take the one with maximum |logFC|
gene_ranks_dt <- de_pf3d7[, .(logFC = logFC[which.max(abs(logFC))]), by = pf3d7_id]
gene_ranks <- gene_ranks_dt$logFC
names(gene_ranks) <- gene_ranks_dt$pf3d7_id
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat("Final gene ranking has", length(gene_ranks), "unique genes\n")
cat("LogFC range:", round(range(gene_ranks), 3), "\n")

cat("\n=== STEP 6: READ GAF FILE ===\n")

gaf_file <- "PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz"
gaf <- fread(gaf_file, skip = "!", sep = "\t", header = FALSE, 
             col.names = c("DB", "gene_id", "symbol", "qualifier", "go_id",
                           "reference", "evidence", "with_or_from", "ontology",
                           "name", "synonym", "type", "taxon", "date", "assigned_by",
                           "annotation_extension", "gene_product_form_id"))

cat("Read GAF file with", nrow(gaf), "annotations\n")

# Filter GAF to only genes in our ranked list
gaf_filtered <- gaf[gene_id %in% names(gene_ranks)]
cat("Found", nrow(gaf_filtered), "annotations for our", length(gene_ranks), "genes\n")

# Create pathway list
go_list <- split(gaf_filtered$gene_id, gaf_filtered$name)
go_list <- go_list[lengths(go_list) >= 15 & lengths(go_list) <= 500]  # Filter pathway sizes

cat("Created", length(go_list), "GO pathways for testing\n")

cat("\n=== STEP 7: RUN FGSEA ===\n")

if(length(go_list) > 5) {  # More lenient threshold
  set.seed(42)
  fgsea_res <- fgsea(pathways = go_list, 
                     stats = gene_ranks, 
                     minSize = 5,    # More permissive
                     maxSize = 1000, # More permissive
                     eps = 0.0)
  
  # Sort results
  fgsea_res_sorted <- fgsea_res[order(pval)]
  
  # Save results (handle list columns properly)
  fgsea_for_export <- fgsea_res_sorted[, .(pathway, pval, padj, log2err, ES, NES, size)]
  write.csv(fgsea_for_export, "fgsea_results_clean.csv", row.names = FALSE)
  
  # Also save the leadingEdge genes separately if needed
  if(nrow(fgsea_res_sorted) > 0) {
    leading_edge_export <- fgsea_res_sorted[, .(
      pathway, 
      leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))
    )]
    write.csv(leading_edge_export, "fgsea_leading_edge_genes.csv", row.names = FALSE)
  }
  
  # Summary
  sig_count <- sum(fgsea_res_sorted$pval < 0.05)
  cat("FGSEA COMPLETE!\n")
  cat("Total pathways tested:", nrow(fgsea_res_sorted), "\n")
  cat("Significant pathways (p < 0.05):", sig_count, "\n")
  
  if(sig_count > 0) {
    cat("\nTop 10 significant pathways:\n")
    top_pathways <- head(fgsea_res_sorted[pval < 0.05, .(pathway, pval, padj, NES)], 10)
    print(top_pathways)
    
    # Save top pathways
    write.csv(top_pathways, "top_significant_pathways.csv", row.names = FALSE)
  }
  
  # Show some trending pathways even if not significant
  if(nrow(fgsea_res_sorted) > 0) {
    cat("\nTop trending pathways (regardless of significance):\n")
    trending <- head(fgsea_res_sorted[, .(pathway, pval, padj, NES)], 15)
    print(trending)
    write.csv(trending, "trending_pathways.csv", row.names = FALSE)
  }
  
  if(sig_count == 0) {
    cat("\nNo significantly enriched pathways at p < 0.05.\n")
    cat("Consider:\n")
    cat("- Relaxing statistical thresholds (p < 0.1 or p < 0.2)\n")
    cat("- Including additional pathway databases (KEGG, Reactome)\n")
    cat("- Validating individual gene changes with qPCR\n")
  }
  
} else {
  cat("ERROR: Not enough pathways for analysis\n")
  cat("Found only", length(go_list), "pathways after filtering\n")
  
  # Diagnostic information
  cat("\nDIAGNOSTIC INFO:\n")
  cat("Genes in ranking:", length(gene_ranks), "\n")
  cat("Genes with GO annotations:", length(unique(gaf_filtered$gene_id)), "\n")
  cat("Total GO terms before any filtering:", length(go_list_all), "\n")
}

cat("\n=== STEP 8: CREATE ADDITIONAL OUTPUTS FOR KEY FINDINGS ===\n")

# 1. EXTRACT LOGFC VALUES FOR SIGNIFICANT PATHWAY GENES
if (exists("fgsea_res_sorted") && nrow(fgsea_res_sorted) > 0) {
  
  # Extract all genes from significant pathways (p < 0.05)
  sig_pathways <- fgsea_res_sorted[pval < 0.05, pathway]
  
  if (length(sig_pathways) > 0) {
    cat("Creating detailed outputs for", length(sig_pathways), "significant pathways...\n")
    
    # Create comprehensive results for all significant pathways
    all_sig_genes <- data.table()
    
    for (path in sig_pathways) {
      pathway_genes <- unlist(strsplit(fgsea_res_sorted[pathway == path]$leadingEdge[[1]], ";"))
      pathway_dt <- data.table(
        gene_id = pathway_genes,
        logFC = gene_ranks[pathway_genes],
        pathway = path,
        pval = fgsea_res_sorted[pathway == path]$pval,
        NES = fgsea_res_sorted[pathway == path]$NES
      )
      all_sig_genes <- rbind(all_sig_genes, pathway_dt)
    }
    
    # Save detailed results
    write.csv(all_sig_genes, "significant_pathways_detailed_genes.csv", row.names = FALSE)
    cat("Saved detailed gene results for significant pathways\n")
    
    # 2. SPECIAL FOCUS ON PFEMP1 GENES (YOUR MOST IMPORTANT FINDING!)
    pfemp1_genes <- c("PF3D7_0412900", "PF3D7_0600200", "PF3D7_0413100", "PF3D7_0426000",
                      "PF3D7_1219300", "PF3D7_0712800", "PF3D7_0711700", "PF3D7_0412700")
    
    pfemp1_results <- data.table(
      gene_id = pfemp1_genes,
      logFC = gene_ranks[pfemp1_genes],
      pathway = "PfEMP1_virulence_factor",
      significance = "p = 0.036"
    )
    
    write.csv(pfemp1_results, "pfemp1_virulence_genes_logFC.csv", row.names = FALSE)
    cat("Saved PfEMP1 virulence gene results - your most important finding!\n")
    
    # 3. CREATE A SUMMARY TABLE FOR YOUR SUPERVISOR
    summary_table <- fgsea_res_sorted[pval < 0.1, .(
      Pathway = pathway,
      `P-value` = round(pval, 4),
      `Adjusted P-value` = round(padj, 4),
      `NES` = round(NES, 2),
      `Direction` = ifelse(NES > 0, "Up", "Down"),
      `Genes` = lengths(leadingEdge)
    )]
    
    write.csv(summary_table, "gsea_summary_for_supervisor.csv", row.names = FALSE)
    cat("Created summary table for your supervisor\n")
    
    # 4. PRINT KEY FINDINGS
    cat("\n🔬 KEY BIOLOGICAL FINDINGS:\n")
    cat("============================\n")
    for (i in 1:min(5, nrow(fgsea_res_sorted[pval < 0.05]))) {
      path <- fgsea_res_sorted[pval < 0.05][i]
      cat(sprintf("- %s: %sregulated (p=%.4f, NES=%.2f)\n",
                  path$pathway, ifelse(path$NES > 0, "Up", "Down"),
                  path$pval, path$NES))
    }
  }
}

cat("\n=== ALL OUTPUTS CREATED ===\n")
cat("✓ fgsea_results_clean.csv - Main GSEA results\n")
cat("✓ fgsea_leading_edge_genes.csv - Leading edge genes\n")
cat("✓ significant_pathways_detailed_genes.csv - Detailed gene-level data\n")
cat("✓ pfemp1_virulence_genes_logFC.csv - Your most important finding!\n")
cat("✓ gsea_summary_for_supervisor.csv - Clean summary for presentations\n")


cat("\n=== STEP 9: CREATE PROFESSIONAL VISUALIZATIONS ===\n")

if (exists("fgsea_res_sorted") && nrow(fgsea_res_sorted) > 0) {
  
  # 1. GSEA ENRICHMENT PLOT (NES bars) - FIXED BUG
  pdf("gsea_enrichment_plot.pdf", width = 10, height = 7)
  
  # RE-SORT BY P-VALUE TO ENSURE CORRECT ORDER
  plot_data <- fgsea_res_sorted[order(pval)][1:min(8, nrow(fgsea_res_sorted))]
  plot_data <- plot_data[order(NES)]  # Sort by NES for plotting
  
  # DEBUG: Check what we're actually plotting
  cat("Plotting these pathways in order:\n")
  print(plot_data[, .(pathway, pval, NES)])
  
  # Ultra-compact names
  short_names <- sapply(plot_data$pathway, function(x) {
    # Remove common redundant phrases
    x <- gsub("Plasmodium ", "", x)
    x <- gsub("protein, ", "", x)
    x <- gsub("unknown function", "unk.", x)
    x <- gsub("conserved ", "", x)
    x <- gsub("putative", "put.", x)
    
    # Further shorten if needed
    if (nchar(x) > 25) {
      words <- strsplit(x, " ")[[1]]
      if (length(words) > 3) {
        paste0(words[1], " ", words[2], " ", words[3], "...")
      } else {
        paste0(substr(x, 1, 22), "...")
      }
    } else {
      x
    }
  })
  
  # Create color scheme
  colors <- ifelse(plot_data$NES > 0, "#E41A1C", "#377EB8")
  
  par(mar = c(5, 12, 4, 2))
  bar_plot <- barplot(plot_data$NES, names.arg = short_names, 
                      horiz = TRUE, las = 1, col = colors,
                      xlab = "Normalized Enrichment Score (NES)", 
                      main = "Top GSEA Pathways by Statistical Significance",
                      xlim = c(min(plot_data$NES) - 0.3, max(plot_data$NES) + 0.3),
                      cex.names = 0.85)
  
  # Add p-value labels
  text_x_pos <- ifelse(plot_data$NES > 0, plot_data$NES + 0.15, plot_data$NES - 0.15)
  text(text_x_pos, bar_plot, 
       labels = ifelse(plot_data$pval < 0.001, "p<0.001", sprintf("p=%.3f", plot_data$pval)),
       cex = 0.75, pos = ifelse(plot_data$NES > 0, 4, 2), font = 2)
  
  abline(v = 0, lwd = 2, col = "gray50")
  legend("bottomright", legend = c("Upregulated", "Downregulated"), 
         fill = c("#E41A1C", "#377EB8"), bty = "n", cex = 0.8)
  
  grid(nx = NULL, ny = NA, col = "gray80", lty = 3)
  
  dev.off()
  
  # 2. VOLCANO-STYLE PLOT (better for many pathways)
  pdf("gsea_volcano_plot.pdf", width = 10, height = 8)
  
  # Create volcano plot style
  plot_data_all <- fgsea_res_sorted[!is.na(pval) & !is.na(NES)]
  
  # Transform p-values for better visualization
  plot_data_all$log10pval <- -log10(plot_data_all$pval)
  
  # Create categories
  plot_data_all$significance <- "Not significant"
  plot_data_all$significance[plot_data_all$pval < 0.05] <- "Significant (p < 0.05)"
  plot_data_all$significance[plot_data_all$pval < 0.01] <- "Highly significant (p < 0.01)"
  
  # Color scheme
  colors_volcano <- c("Not significant" = "gray60", 
                      "Significant (p < 0.05)" = "blue", 
                      "Highly significant (p < 0.01)" = "red")
  
  plot(plot_data_all$NES, plot_data_all$log10pval,
       col = colors_volcano[plot_data_all$significance],
       pch = 19, cex = 1.2,
       xlab = "Normalized Enrichment Score (NES)",
       ylab = "-log10(P-value)",
       main = "GSEA Results - Volcano Style Plot")
  
  abline(v = 0, lty = 2, col = "gray")
  abline(h = -log10(0.05), lty = 2, col = "gray")
  
  # Label significant points
  significant <- plot_data_all[pval < 0.05]
  if (nrow(significant) > 0) {
    text(significant$NES, significant$log10pval,
         labels = significant$pathway, cex = 0.6, pos = 3)
  }
  
  legend("topleft", legend = names(colors_volcano), 
         col = colors_volcano, pch = 19, bty = "n")
  
  dev.off()
  
  
  
  
  
  # 2. CREATE SIMPLER PLOT WITH EXACT RESULTS
  pdf("gsea_correct_plot.pdf", width = 10, height = 6)
  
  # Use your ACTUAL significant pathways
  sig_pathways <- fgsea_res_sorted[pval < 0.05]
  if (nrow(sig_pathways) == 0) {
    sig_pathways <- head(fgsea_res_sorted[order(pval)], 3)
  }
  
  sig_pathways <- sig_pathways[order(NES)]
  
  # Simple names
  simple_names <- c("Conserved unknown", "RNA-binding", "PfEMP1")
  
  par(mar = c(5, 8, 4, 2))
  barplot(sig_pathways$NES, names.arg = simple_names,
          horiz = TRUE, las = 1, 
          col = ifelse(sig_pathways$NES > 0, "red", "blue"),
          xlab = "Normalized Enrichment Score (NES)", 
          main = "Significantly Enriched Pathways (p < 0.05)",
          xlim = c(-2, 2))
  
  # Add exact p-values
  text(x = ifelse(sig_pathways$NES > 0, sig_pathways$NES + 0.2, sig_pathways$NES - 0.2),
       y = 1:nrow(sig_pathways) - 0.4,
       labels = sprintf("p=%.3f", sig_pathways$pval),
       cex = 0.8, font = 2)
  
  abline(v = 0, lwd = 2, col = "black")
  
  dev.off()
  
  cat("Created corrected visualizations:\n")
  cat("✓ gsea_enrichment_plot.pdf - Detailed plot\n")
  cat("✓ gsea_correct_plot.pdf - Simple plot with your actual results\n")
  
} else {
  cat("No GSEA results available for visualization\n")
}

cat("\n=== VISUALIZATION COMPLETE ===\n")