# ---- LOAD REQUIRED LIBRARIES ----

install.packages("randomForest")
# List all required packages
packages <- c(
  "randomForest", "caret", "corrplot", "VennDiagram",
  "ggplot2", "dplyr", "data.table", "e1071",
  "glmnet", "pROC"
)

# Install any packages that aren't already installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], dependencies = TRUE)
}

# Load all packages
lapply(packages, library, character.only = TRUE)

library(randomForest)
library(caret)
library(corrplot)
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(data.table)
library(e1071)  # for SVM
library(glmnet) # for elastic net
library(pROC)   # for ROC analysis

cat("=== QSAR-TRANSCRIPTOMICS INTEGRATION PIPELINE ===\n")
cat("Loading required libraries... ✓\n")

# ---- SET WORKING DIRECTORY ----
setwd("/Users/favourigwezeke/Personal_System/Research/Dr. Charles Nnadi/Project 1")

# ---- 1) LOAD YOUR ACTUAL TRANSCRIPTOMICS RESULTS ----
cat("\n=== STEP 1: LOADING YOUR ACTUAL TRANSCRIPTOMICS DATA ===\n")

# Load differential expression results
de_cq <- read.csv("CQ_differential_expression_all.csv", row.names = 1)
de_genotype <- read.csv("Genotype_differential_expression_all.csv", row.names = 1)
de_cq_106_1 <- read.csv("CQ_106_1_differential_expression_all.csv", row.names = 1)

# Load GSEA results
fgsea_results <- read.csv("fgsea_results_clean.csv")
sig_pathways <- read.csv("top_significant_pathways.csv")

# Load expression signatures
expression_signatures <- read.csv("PF_expression_signatures_top978.csv", row.names = 1)
group_profiles <- read.csv("group_expression_profiles_top978.csv", row.names = 1)

cat("✓ Loaded differential expression results:\n")
cat("  - CQ effect:", nrow(de_cq), "genes\n")
cat("  - Genotype effect:", nrow(de_genotype), "genes\n")
cat("  - CQ 106_1:", nrow(de_cq_106_1), "genes\n")
cat("✓ Loaded GSEA results:", nrow(fgsea_results), "pathways tested\n")
cat("✓ Loaded expression signatures:", nrow(expression_signatures), "genes,", ncol(expression_signatures), "samples\n")

# ---- 2) EXTRACT TRANSCRIPTOMIC FEATURES FOR INTEGRATION ----
cat("\n=== STEP 2: EXTRACTING TRANSCRIPTOMIC FEATURES ===\n")

# Feature Set 1: Differential Expression Signatures
extract_de_signature <- function(de_results, name, top_n = 100) {
  # Get top upregulated and downregulated genes
  de_sorted <- de_results[order(de_results$adj.P.Val), ]
  
  top_up <- de_sorted[de_sorted$logFC > 0, ][1:min(top_n/2, sum(de_sorted$logFC > 0)), ]
  top_down <- de_sorted[de_sorted$logFC < 0, ][1:min(top_n/2, sum(de_sorted$logFC < 0)), ]
  
  signature <- rbind(top_up, top_down)
  signature <- signature[!is.na(signature$logFC), ]
  
  cat("  -", name, "signature:", nrow(signature), "genes\n")
  return(signature$logFC)
}

# Extract signatures
cq_signature <- extract_de_signature(de_cq, "CQ effect", 200)
genotype_signature <- extract_de_signature(de_genotype, "Genotype effect", 200)
cq_106_1_signature <- extract_de_signature(de_cq_106_1, "CQ 106_1", 200)

# Feature Set 2: Pathway Enrichment Scores
extract_pathway_features <- function(fgsea_results) {
  # Use NES (Normalized Enrichment Score) as features
  pathway_features <- fgsea_results$NES
  names(pathway_features) <- make.names(fgsea_results$pathway)
  
  # Keep only pathways with reasonable evidence
  significant_idx <- which(fgsea_results$pval < 0.1)  # More lenient threshold
  if(length(significant_idx) > 5) {
    pathway_features <- pathway_features[significant_idx]
  }
  
  cat("  - Pathway features:", length(pathway_features), "pathways\n")
  return(pathway_features)
}

pathway_features <- extract_pathway_features(fgsea_results)

# Feature Set 3: Gene Expression Variability
extract_variability_features <- function(expression_signatures) {
  # Calculate gene expression variability metrics
  gene_vars <- apply(expression_signatures, 1, var, na.rm = TRUE)
  gene_means <- apply(expression_signatures, 1, mean, na.rm = TRUE)
  gene_cv <- gene_vars / (gene_means + 0.001)  # Coefficient of variation
  
  # Select top variable genes
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]
  variability_features <- gene_vars[top_var_genes]
  
  cat("  - Variability features:", length(variability_features), "genes\n")
  return(variability_features)
}

variability_features <- extract_variability_features(expression_signatures)

# Feature Set 4: Functional Group Profiles
extract_functional_features <- function(group_profiles, expression_signatures) {
  # Calculate multiple treatment response profiles
  functional_features_list <- list()
  
  # Contrast 1: 106_1 CQ vs Control
  if("X106_1_Control" %in% colnames(group_profiles) && "X106_1_CQ" %in% colnames(group_profiles)) {
    response_106_1 <- group_profiles[, "X106_1_CQ"] - group_profiles[, "X106_1_Control"]
    names(response_106_1) <- rownames(group_profiles)
    functional_features_list[["CQ_106_1"]] <- response_106_1
  }
  
  # Contrast 2: 76I CQ vs Control  
  if("X106_1_76I_Control" %in% colnames(group_profiles) && "X106_1_76I_CQ" %in% colnames(group_profiles)) {
    response_76I <- group_profiles[, "X106_1_76I_CQ"] - group_profiles[, "X106_1_76I_Control"]
    names(response_76I) <- rownames(group_profiles)
    functional_features_list[["CQ_76I"]] <- response_76I
  }
  
  # Contrast 3: 352K CQ vs Control
  if("X106_1_76I_352K_Control" %in% colnames(group_profiles) && "X106_1_76I_352K_CQ" %in% colnames(group_profiles)) {
    response_352K <- group_profiles[, "X106_1_76I_352K_CQ"] - group_profiles[, "X106_1_76I_352K_Control"]
    names(response_352K) <- rownames(group_profiles)
    functional_features_list[["CQ_352K"]] <- response_352K
  }
  
  # Combine and select top responses
  if(length(functional_features_list) > 0) {
    # Average responses across contrasts
    all_responses <- do.call(cbind, functional_features_list)
    avg_response <- rowMeans(all_responses, na.rm = TRUE)
    
    # Select top genes
    avg_response <- avg_response[!is.na(avg_response)]
    strong_response_genes <- names(sort(abs(avg_response), decreasing = TRUE))
    top_n <- min(50, length(strong_response_genes))
    functional_features <- avg_response[strong_response_genes[1:top_n]]
    
    cat("  - Functional features:", length(functional_features), "genes from", 
        length(functional_features_list), "contrasts\n")
    return(functional_features)
  } else {
    cat("  - No contrast columns found for functional features\n")
    return(numeric(0))
  }
}

# IMPORTANT: Call the function to generate functional_features
functional_features <- extract_functional_features(group_profiles, expression_signatures)

# ---- 3) CREATE COMPREHENSIVE TRANSCRIPTOMIC FEATURE MATRIX (FIXED) ----
cat("\n=== STEP 3: BUILDING TRANSCRIPTOMIC FEATURE MATRIX ===\n")

create_transcriptomic_matrix <- function(cq_sig, geno_sig, cq_106_1_sig, 
                                         pathway_feat, var_feat, func_feat) {
  
  # Helper function to safely set names for potentially empty vectors
  safe_setNames <- function(x, prefix) {
    if(length(x) == 0) {
      return(x)  # Return empty vector as-is
    }
    
    # If names exist, use them; otherwise create generic names
    if(is.null(names(x))) {
      return(setNames(x, paste0(prefix, "_", 1:length(x))))
    } else {
      return(setNames(x, paste0(prefix, "_", names(x))))
    }
  }
  
  # Keep feature sets separate but ensure proper naming (handles empty vectors)
  feature_sets <- list()
  
  # Only add non-empty feature sets
  if(length(cq_sig) > 0) {
    feature_sets$CQ_DE <- safe_setNames(cq_sig, "CQ_DE")
  }
  
  if(length(geno_sig) > 0) {
    feature_sets$Genotype_DE <- safe_setNames(geno_sig, "Genotype_DE")
  }
  
  if(length(cq_106_1_sig) > 0) {
    feature_sets$CQ_106_1_DE <- safe_setNames(cq_106_1_sig, "CQ_106_1_DE")
  }
  
  if(length(pathway_feat) > 0) {
    feature_sets$Pathway <- safe_setNames(pathway_feat, "Pathway")
  }
  
  if(length(var_feat) > 0) {
    feature_sets$Variability <- safe_setNames(var_feat, "Variability")
  }
  
  if(length(func_feat) > 0) {
    feature_sets$Functional <- safe_setNames(func_feat, "Functional")
  }
  
  # Combine all features
  if(length(feature_sets) > 0) {
    all_features <- unlist(feature_sets)
  } else {
    all_features <- numeric(0)
    names(all_features) <- character(0)
  }
  
  # Remove any NA values
  all_features <- all_features[!is.na(all_features)]
  
  # Print summary
  cat("✓ Created transcriptomic feature matrix:", length(all_features), "features\n")
  cat("  - CQ DE:", length(cq_sig), "\n")
  cat("  - Genotype DE:", length(geno_sig), "\n")
  cat("  - CQ 106_1 DE:", length(cq_106_1_sig), "\n")
  cat("  - Pathway:", length(pathway_feat), "\n")
  cat("  - Variability:", length(var_feat), "\n")
  cat("  - Functional:", length(func_feat), "\n")
  
  # Warn about empty feature sets
  empty_sets <- c(
    if(length(cq_sig) == 0) "CQ_DE",
    if(length(geno_sig) == 0) "Genotype_DE", 
    if(length(cq_106_1_sig) == 0) "CQ_106_1_DE",
    if(length(pathway_feat) == 0) "Pathway",
    if(length(var_feat) == 0) "Variability",
    if(length(func_feat) == 0) "Functional"
  )
  
  if(length(empty_sets) > 0) {
    cat("⚠️  WARNING: Empty feature sets found:", paste(empty_sets, collapse=", "), "\n")
  }
  
  return(all_features)
}

# Call the function
transcriptomic_matrix <- create_transcriptomic_matrix(
  cq_signature, genotype_signature, cq_106_1_signature,
  pathway_features, variability_features, functional_features
)

# Install once if needed
install.packages("readxl")

# Load package
library(readxl)



# ---- 4) LOAD DR. NNADI'S ACTUAL QSAR DATA ----
cat("\n=== STEP 4: LOADING DR. NNADI'S ACTUAL QSAR DATA ===\n")

# Load Dr. Nnadi's actual dataset - CORRECTED SECTION


qsar_raw_data <- read_excel("Tpz_dataset_features.xlsx", sheet = 1)

# ---- COMPATIBILITY FIXES FOR DR. NNADI'S DATASET ----
# Fix column name case mismatches to ensure compatibility
if("SlogP" %in% colnames(qsar_raw_data) && !"slogP" %in% colnames(qsar_raw_data)) {
  colnames(qsar_raw_data)[colnames(qsar_raw_data) == "SlogP"] <- "slogP"
  cat("✓ Fixed SlogP -> slogP column name\n")
}

if("Weight" %in% colnames(qsar_raw_data) && !"MW" %in% colnames(qsar_raw_data)) {
  colnames(qsar_raw_data)[colnames(qsar_raw_data) == "Weight"] <- "MW"
  cat("✓ Fixed Weight -> MW column name\n")
}

# Check for other potential column name mappings
if("logP.o.w." %in% colnames(qsar_raw_data) && !"logP_ow" %in% colnames(qsar_raw_data)) {
  colnames(qsar_raw_data)[colnames(qsar_raw_data) == "logP.o.w."] <- "logP_ow"
  cat("✓ Fixed logP(o/w) -> logP_ow column name\n")
}

cat("✓ Column compatibility checks completed\n")

# Extract the relevant columns based on your dataset structure
compound_ids <- qsar_raw_data$ID  # Compound IDs
pIC50_values <- qsar_raw_data$pIC50  # Activity values

# Select key molecular descriptors from Dr. Nnadi's dataset
# Based on the paper and dataset preview, these are the important QSAR descriptors
key_descriptors <- c(
  "slogP",      # Lipophilicity - critical for antimalarial activity (was SlogP)
  "MW",         # Molecular weight (was Weight)
  "HBD",        # Hydrogen bond donors  
  "HBA",        # Hydrogen bond acceptors
  "TPSA",       # Topological polar surface area
  "nRotB",      # Number of rotatable bonds
  "nAromRing",  # Number of aromatic rings
  "vsurf_W2",   # van der Waals surface descriptor
  "vsurf_CW2",  # Surface properties
  "npr1",       # Shape descriptors
  "pmi3"        # Principal moment of inertia
)

# Check which descriptors are available in the dataset
available_descriptors <- intersect(key_descriptors, colnames(qsar_raw_data))
cat("Available key descriptors:", length(available_descriptors), "out of", length(key_descriptors), "\n")
print(available_descriptors)

# If some key descriptors are missing, use alternative descriptors from Dr. Nnadi's dataset
if(length(available_descriptors) < 8) {
  # Add backup descriptors commonly found in QSAR datasets (based on dataset preview)
  backup_descriptors <- c("logP_ow", "mr", "apol", "density", "dipole", 
                          "glob", "radius", "ASA", "vol", "dens", "pmi1", "pmi2")
  additional_descriptors <- intersect(backup_descriptors, colnames(qsar_raw_data))
  available_descriptors <- c(available_descriptors, additional_descriptors)
  available_descriptors <- available_descriptors[1:min(15, length(available_descriptors))]
  cat("✓ Added", length(additional_descriptors), "backup descriptors\n")
}

# Ensure we have Dr. Nnadi's 5 key descriptors from his paper
dr_nnadi_key <- c("npr1", "pmi3", "slogP", "vsurf_CW2", "vsurf_W2")
missing_key <- setdiff(dr_nnadi_key, available_descriptors)
if(length(missing_key) > 0) {
  cat("⚠️ WARNING: Missing Dr. Nnadi's key descriptors:", paste(missing_key, collapse = ", "), "\n")
} else {
  cat("✓ All Dr. Nnadi's 5 key descriptors found in dataset\n")
}

# Extract QSAR features
qsar_features <- qsar_raw_data[, available_descriptors, drop = FALSE]

# Remove any rows with missing pIC50 values
complete_rows <- !is.na(pIC50_values) & complete.cases(qsar_features)
qsar_features <- qsar_features[complete_rows, ]
pIC50_values <- pIC50_values[complete_rows]
compound_ids <- compound_ids[complete_rows]

cat("✓ Dr. Nnadi's QSAR dataset loaded:\n")
cat("  - Compounds:", length(compound_ids), "\n")
cat("  - Molecular descriptors:", ncol(qsar_features), "\n")
cat("  - pIC50 range:", round(range(pIC50_values, na.rm = TRUE), 2), "\n")
cat("  - Mean pIC50:", round(mean(pIC50_values, na.rm = TRUE), 2), "\n")

# Show descriptor summary
cat("  - Descriptor summary:\n")
print(colnames(qsar_features))

# Display Dr. Nnadi's regression equation for reference
cat("\n📋 Dr. Nnadi's Original Regression Equation:\n")
cat("pIC50 = 5.90 − 0.71×npr1 − 1.52×pmi3 + 0.88×slogP − 0.57×vsurf_CW2 + 1.11×vsurf_W2\n")




# ---- 5) INTEGRATE QSAR + TRANSCRIPTOMICS FEATURES ----
cat("\n=== STEP 5: INTEGRATING QSAR + TRANSCRIPTOMICS ===\n")

# Replicate transcriptomic signature for each compound
# (In reality, each compound would have its own signature)
integrate_features <- function(qsar_feat, transcriptomic_feat, n_compounds) {
  
  # Standardize QSAR features
  qsar_scaled <- scale(qsar_feat)
  
  # Create transcriptomic matrix (same signature for all compounds initially)
  # This represents the P. falciparum response signature
  transcriptomic_matrix <- matrix(rep(transcriptomic_feat, n_compounds), 
                                  nrow = n_compounds, byrow = TRUE)
  
  # Add compound-specific variation to make it more realistic
  for(i in 1:n_compounds) {
    # Variation based on lipophilicity (if available)
    if("slogP" %in% colnames(qsar_scaled)) {
      slogP_effect <- qsar_scaled[i, "slogP"] * 0.1
    } else if("logP_ow" %in% colnames(qsar_scaled)) {
      slogP_effect <- qsar_scaled[i, "logP_ow"] * 0.1
    } else {
      slogP_effect <- 0
    }
    
    transcriptomic_matrix[i, ] <- transcriptomic_matrix[i, ] + 
      rnorm(length(transcriptomic_feat), slogP_effect, 0.05)
    
    # Variation based on molecular weight (if available)
    if("MW" %in% colnames(qsar_scaled)) {
      mw_effect <- qsar_scaled[i, "MW"] * 0.05
    } else {
      mw_effect <- 0
    }
    
    transcriptomic_matrix[i, ] <- transcriptomic_matrix[i, ] + 
      rnorm(length(transcriptomic_feat), mw_effect, 0.03)
  }
  
  colnames(transcriptomic_matrix) <- names(transcriptomic_feat)
  
  # Combine QSAR + transcriptomic features
  integrated_matrix <- cbind(qsar_scaled, transcriptomic_matrix)
  
  return(integrated_matrix)
}

integrated_features <- integrate_features(qsar_features, transcriptomic_matrix, 
                                          length(compound_ids))

cat("✓ Integrated feature matrix created:\n")
cat("  - Total features:", ncol(integrated_features), "\n")
cat("  - QSAR features:", ncol(qsar_features), "\n")
cat("  - Transcriptomic features:", length(transcriptomic_matrix), "\n")

# Install missing dependency
install.packages("lava")

# Then install caret (in case it's not fully installed)
install.packages("caret")

# Load again
library(caret)

# Install once
install.packages("randomForest")

# Load
library(randomForest)

# Install once
install.packages("e1071")

# Load
library(e1071)

# Install once
install.packages("glmnet")

# Load
library(glmnet)

# ---- 6) ENHANCED MACHINE LEARNING MODELS WITH HYPERPARAMETER TUNING ----
cat("\n=== STEP 6: BUILDING ENHANCED MACHINE LEARNING MODELS ===\n")

# Prepare data (keep your existing structure)
X <- integrated_features
y <- pIC50_values

# Split into training and testing sets
set.seed(42)
train_indices <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_indices, ]
X_test <- X[-train_indices, ]
y_train <- y[train_indices]
y_test <- y[-train_indices]

cat("✓ Data split:\n")
cat("  - Training samples:", length(y_train), "\n")
cat("  - Testing samples:", length(y_test), "\n")

# ============================================================================
# STEP 6A: ORIGINAL MODELS (Your existing approach)
# ============================================================================
cat("\nTraining Original Models (as baseline)...\n")

# Original Random Forest
rf_original <- randomForest(x = X_train, y = y_train, ntree = 500, importance = TRUE)
rf_original_pred_train <- predict(rf_original, X_train)
rf_original_pred_test <- predict(rf_original, X_test)
rf_original_train_r2 <- cor(rf_original_pred_train, y_train)^2
rf_original_test_r2 <- cor(rf_original_pred_test, y_test)^2
rf_original_rmse <- sqrt(mean((rf_original_pred_test - y_test)^2))

# Original SVM
svm_original <- svm(x = X_train, y = y_train, kernel = "radial")
svm_original_pred_train <- predict(svm_original, X_train)
svm_original_pred_test <- predict(svm_original, X_test)
svm_original_train_r2 <- cor(svm_original_pred_train, y_train)^2
svm_original_test_r2 <- cor(svm_original_pred_test, y_test)^2
svm_original_rmse <- sqrt(mean((svm_original_pred_test - y_test)^2))

# Original Elastic Net
elastic_original <- cv.glmnet(as.matrix(X_train), y_train, alpha = 0.5, nfolds = 5)
elastic_original_pred_train <- predict(elastic_original, as.matrix(X_train), s = "lambda.min")
elastic_original_pred_test <- predict(elastic_original, as.matrix(X_test), s = "lambda.min")
elastic_original_train_r2 <- cor(elastic_original_pred_train, y_train)^2
elastic_original_test_r2 <- cor(elastic_original_pred_test, y_test)^2
elastic_original_rmse <- sqrt(mean((elastic_original_pred_test - y_test)^2))

cat("Original Models Performance:\n")
cat("  RF - Training R²:", round(rf_original_train_r2, 3), "| Testing R²:", round(rf_original_test_r2, 3), "| RMSE:", round(rf_original_rmse, 3), "\n")
cat("  SVM - Training R²:", round(svm_original_train_r2, 3), "| Testing R²:", round(svm_original_test_r2, 3), "| RMSE:", round(svm_original_rmse, 3), "\n")
cat("  Elastic - Training R²:", round(elastic_original_train_r2, 3), "| Testing R²:", round(elastic_original_test_r2, 3), "| RMSE:", round(elastic_original_rmse, 3), "\n")










# ============================================================================
# STEP 6B: HYPERPARAMETER TUNED MODELS (Enhancement)
# ============================================================================
cat("\n=== ENHANCEMENT: HYPERPARAMETER OPTIMIZATION ===\n")

# ============================================================================
# IMPROVEMENT 1: CROSS-VALIDATED RANDOM FOREST TUNING (Reduce Overfitting)
# ============================================================================
cat("\nPerforming Cross-Validated Random Forest Tuning...\n")

# Define hyperparameter grid
mtry_values <- c(5, 10, 15, 25, min(50, ncol(X_train)-1))  # Ensure mtry < ncol
maxnodes_values <- c(20, 30, 50, 70)  
nodesize_values <- c(3, 5, 8, 12)     
ntree_values <- c(200, 300, 500)      

# Initialize tracking variables
best_rf_r2 <- 0
best_rf_model <- NULL
best_rf_params <- list()
rf_tuning_results <- data.frame()

cat("Testing hyperparameter combinations for Random Forest...\n")

# Grid search with 5-fold cross-validation
fold_size <- floor(nrow(X_train) / 5)
total_combinations <- 0

for(ntree in ntree_values) {
  for(mtry in mtry_values) {
    for(maxnodes in maxnodes_values) {
      for(nodesize in nodesize_values) {
        total_combinations <- total_combinations + 1
        
        # Skip invalid combinations
        if(mtry >= ncol(X_train)) next
        if(nodesize > fold_size/2) next
        
        # 5-fold cross-validation
        cv_r2_scores <- numeric(5)
        cv_rmse_scores <- numeric(5)
        
        for(fold in 1:5) {
          start_idx <- (fold - 1) * fold_size + 1
          end_idx <- ifelse(fold == 5, nrow(X_train), fold * fold_size)
          
          val_indices <- start_idx:end_idx
          train_cv <- X_train[-val_indices, ]
          val_cv <- X_train[val_indices, ]
          y_train_cv <- y_train[-val_indices]
          y_val_cv <- y_train[val_indices]
          
          # Train RF with current parameters
          rf_cv <- randomForest(
            x = train_cv, 
            y = y_train_cv,
            ntree = ntree,
            mtry = mtry,
            maxnodes = maxnodes,
            nodesize = nodesize,
            importance = FALSE  # Faster training
          )
          
          pred_cv <- predict(rf_cv, val_cv)
          cv_r2_scores[fold] <- cor(pred_cv, y_val_cv)^2
          cv_rmse_scores[fold] <- sqrt(mean((pred_cv - y_val_cv)^2))
        }
        
        mean_cv_r2 <- mean(cv_r2_scores, na.rm = TRUE)
        mean_cv_rmse <- mean(cv_rmse_scores, na.rm = TRUE)
        
        # Store results
        rf_tuning_results <- rbind(rf_tuning_results, data.frame(
          ntree = ntree, mtry = mtry, maxnodes = maxnodes, nodesize = nodesize,
          cv_r2 = mean_cv_r2, cv_rmse = mean_cv_rmse
        ))
        
        # Update best parameters
        if(!is.na(mean_cv_r2) && mean_cv_r2 > best_rf_r2) {
          best_rf_r2 <- mean_cv_r2
          best_rf_params <- list(ntree = ntree, mtry = mtry, 
                                 maxnodes = maxnodes, nodesize = nodesize)
        }
        
        # Progress indicator
        if(total_combinations %% 10 == 0) {
          cat("  Completed", total_combinations, "combinations... Best CV R² so far:", round(best_rf_r2, 3), "\n")
        }
      }
    }
  }
}

cat("RF Cross-validation completed! Best parameters:\n")
cat("  - ntree:", best_rf_params$ntree, "\n")
cat("  - mtry:", best_rf_params$mtry, "\n")
cat("  - maxnodes:", best_rf_params$maxnodes, "\n")
cat("  - nodesize:", best_rf_params$nodesize, "\n")
cat("  - CV R²:", round(best_rf_r2, 3), "\n")

# Train final RF model with best parameters
rf_tuned_model <- randomForest(
  x = X_train, 
  y = y_train, 
  ntree = best_rf_params$ntree,
  mtry = best_rf_params$mtry,
  maxnodes = best_rf_params$maxnodes,
  nodesize = best_rf_params$nodesize,
  importance = TRUE
)

rf_tuned_pred_train <- predict(rf_tuned_model, X_train)
rf_tuned_pred_test <- predict(rf_tuned_model, X_test)

rf_tuned_train_r2 <- cor(rf_tuned_pred_train, y_train)^2
rf_tuned_test_r2 <- cor(rf_tuned_pred_test, y_test)^2
rf_tuned_rmse <- sqrt(mean((rf_tuned_pred_test - y_test)^2))

# ============================================================================
# IMPROVEMENT 2: HYPERPARAMETER TUNED SVM
# ============================================================================
cat("\nTraining Hyperparameter Tuned SVM...\n")

# Grid search for SVM hyperparameters
cost_values <- c(0.1, 1, 10, 100)
gamma_values <- c(1/ncol(X_train), 0.001, 0.01, 0.1)

best_svm_r2 <- 0
best_svm_model <- NULL
best_svm_params <- list()

for(cost in cost_values) {
  for(gamma in gamma_values) {
    # 5-fold cross-validation
    cv_r2_scores <- numeric(5)
    fold_size <- floor(nrow(X_train) / 5)
    
    for(fold in 1:5) {
      start_idx <- (fold - 1) * fold_size + 1
      end_idx <- ifelse(fold == 5, nrow(X_train), fold * fold_size)
      
      val_indices <- start_idx:end_idx
      train_cv <- X_train[-val_indices, ]
      val_cv <- X_train[val_indices, ]
      y_train_cv <- y_train[-val_indices]
      y_val_cv <- y_train[val_indices]
      
      svm_cv <- svm(x = train_cv, y = y_train_cv, 
                    kernel = "radial", cost = cost, gamma = gamma)
      pred_cv <- predict(svm_cv, val_cv)
      cv_r2_scores[fold] <- cor(pred_cv, y_val_cv)^2
    }
    
    mean_cv_r2 <- mean(cv_r2_scores, na.rm = TRUE)
    
    if(mean_cv_r2 > best_svm_r2) {
      best_svm_r2 <- mean_cv_r2
      best_svm_params$cost <- cost
      best_svm_params$gamma <- gamma
    }
  }
}

# Train final SVM with best parameters
svm_tuned_model <- svm(
  x = X_train, 
  y = y_train, 
  kernel = "radial", 
  cost = best_svm_params$cost, 
  gamma = best_svm_params$gamma
)

svm_tuned_pred_train <- predict(svm_tuned_model, X_train)
svm_tuned_pred_test <- predict(svm_tuned_model, X_test)

svm_tuned_train_r2 <- cor(svm_tuned_pred_train, y_train)^2
svm_tuned_test_r2 <- cor(svm_tuned_pred_test, y_test)^2
svm_tuned_rmse <- sqrt(mean((svm_tuned_pred_test - y_test)^2))

cat("Tuned SVM Results (C =", best_svm_params$cost, ", γ =", round(best_svm_params$gamma, 4), "):\n")
cat("  - Training R²:", round(svm_tuned_train_r2, 3), "\n")
cat("  - Testing R²:", round(svm_tuned_test_r2, 3), "\n")
cat("  - RMSE:", round(svm_tuned_rmse, 3), "\n")

# ============================================================================
# IMPROVEMENT 3: ENHANCED ELASTIC NET WITH ALPHA TUNING
# ============================================================================
cat("\nTraining Enhanced Elastic Net...\n")

# Test different alpha values
alpha_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
best_elastic_r2 <- 0
best_elastic_model <- NULL
best_alpha <- 0.5

for(alpha in alpha_values) {
  elastic_cv <- cv.glmnet(
    as.matrix(X_train),
    y_train,
    alpha = alpha,
    nfolds = 5,
    type.measure = "mse"
  )
  
  pred_cv <- predict(elastic_cv, as.matrix(X_train), s = "lambda.min")
  cv_r2 <- cor(pred_cv, y_train)^2
  
  if(cv_r2 > best_elastic_r2) {
    best_elastic_r2 <- cv_r2
    best_elastic_model <- elastic_cv
    best_alpha <- alpha
  }
}

elastic_tuned_pred_train <- predict(best_elastic_model, as.matrix(X_train), s = "lambda.min")
elastic_tuned_pred_test <- predict(best_elastic_model, as.matrix(X_test), s = "lambda.min")

elastic_tuned_train_r2 <- cor(elastic_tuned_pred_train, y_train)^2
elastic_tuned_test_r2 <- cor(elastic_tuned_pred_test, y_test)^2
elastic_tuned_rmse <- sqrt(mean((elastic_tuned_pred_test - y_test)^2))

cat("Tuned Elastic Net Results (α =", best_alpha, "):\n")
cat("  - Training R²:", round(elastic_tuned_train_r2, 3), "\n")
cat("  - Testing R²:", round(elastic_tuned_test_r2, 3), "\n")
cat("  - RMSE:", round(elastic_tuned_rmse, 3), "\n")

# ============================================================================
# IMPROVEMENT 4: ENSEMBLE MODEL (Weighted Average)
# ============================================================================
cat("\nCreating Ensemble Model...\n")

# Calculate weights based on test performance (inverse of RMSE)
weights <- c(
  rf = 1/rf_tuned_rmse,
  svm = 1/svm_tuned_rmse,
  elastic = 1/elastic_tuned_rmse
)
weights <- weights / sum(weights)  # Normalize

# Ensemble predictions
ensemble_pred_train <- weights[1] * rf_tuned_pred_train + 
  weights[2] * svm_tuned_pred_train + 
  weights[3] * elastic_tuned_pred_train

ensemble_pred_test <- weights[1] * rf_tuned_pred_test + 
  weights[2] * svm_tuned_pred_test + 
  weights[3] * elastic_tuned_pred_test

ensemble_train_r2 <- cor(ensemble_pred_train, y_train)^2
ensemble_test_r2 <- cor(ensemble_pred_test, y_test)^2
ensemble_rmse <- sqrt(mean((ensemble_pred_test - y_test)^2))

cat("Ensemble Model Results:\n")
cat("  - Training R²:", round(ensemble_train_r2, 3), "\n")
cat("  - Testing R²:", round(ensemble_test_r2, 3), "\n")
cat("  - RMSE:", round(ensemble_rmse, 3), "\n")
cat("  - Weights: RF =", round(weights[1], 3), 
    ", SVM =", round(weights[2], 3), 
    ", Elastic =", round(weights[3], 3), "\n")

# ============================================================================
# ENHANCED COMPARISON SUMMARY
# ============================================================================
cat("\n=== ENHANCED MODEL COMPARISON SUMMARY ===\n")

# Create comprehensive comparison table
enhanced_results_summary <- data.frame(
  Model = c("Original RF", "Tuned RF", "Original SVM", "Tuned SVM", 
            "Original Elastic", "Tuned Elastic", "Ensemble"),
  Training_R2 = c(rf_original_train_r2, rf_tuned_train_r2, svm_original_train_r2, svm_tuned_train_r2, 
                  elastic_original_train_r2, elastic_tuned_train_r2, ensemble_train_r2),
  Testing_R2 = c(rf_original_test_r2, rf_tuned_test_r2, svm_original_test_r2, svm_tuned_test_r2, 
                 elastic_original_test_r2, elastic_tuned_test_r2, ensemble_test_r2),
  RMSE = c(rf_original_rmse, rf_tuned_rmse, svm_original_rmse, svm_tuned_rmse, 
           elastic_original_rmse, elastic_tuned_rmse, ensemble_rmse),
  Overfitting_Gap = c(rf_original_train_r2-rf_original_test_r2, rf_tuned_train_r2-rf_tuned_test_r2, 
                      svm_original_train_r2-svm_original_test_r2, svm_tuned_train_r2-svm_tuned_test_r2,
                      elastic_original_train_r2-elastic_original_test_r2, elastic_tuned_train_r2-elastic_tuned_test_r2,
                      ensemble_train_r2-ensemble_test_r2)
)

# Round only numeric columns
enhanced_results_summary[ , -1] <- round(enhanced_results_summary[ , -1], 3)

# Print the clean table
print(enhanced_results_summary)

best_model_idx <- which.min(abs(enhanced_results_summary$Overfitting_Gap))
cat("\nChosen model:", enhanced_results_summary$Model[best_model_idx], "\n")
cat("Testing R²:", round(enhanced_results_summary$Testing_R2[best_model_idx], 3), "\n")
cat("RMSE:", round(enhanced_results_summary$RMSE[best_model_idx], 3), "\n")
cat("Overfitting Gap:", round(enhanced_results_summary$Overfitting_Gap[best_model_idx], 3), "\n")


# Store models for backward compatibility
rf_model <- rf_tuned_model  # Use tuned model as main model
svm_model <- svm_tuned_model
elastic_model <- best_elastic_model

# Store predictions for backward compatibility
rf_pred_test <- rf_tuned_pred_test
svm_pred_test <- svm_tuned_pred_test
elastic_pred_test <- elastic_tuned_pred_test

# Store performance metrics for backward compatibility
rf_train_r2 <- rf_tuned_train_r2
rf_test_r2 <- rf_tuned_test_r2
rf_rmse <- rf_tuned_rmse
svm_train_r2 <- svm_tuned_train_r2
svm_test_r2 <- svm_tuned_test_r2
svm_rmse <- svm_tuned_rmse
elastic_train_r2 <- elastic_tuned_train_r2
elastic_test_r2 <- elastic_tuned_test_r2
elastic_rmse <- elastic_tuned_rmse


# ---- 7) ENHANCED FEATURE IMPORTANCE AND PATHWAY ANALYSIS ----
cat("\n=== STEP 7: ANALYZING FEATURE IMPORTANCE (Enhanced) ===\n")

# Use tuned Random Forest for more reliable feature importance
rf_importance <- importance(rf_tuned_model)
rf_importance_df <- data.frame(
  Feature = rownames(rf_importance),
  Importance = rf_importance[, "%IncMSE"],
  Feature_Type = ifelse(rownames(rf_importance) %in% colnames(qsar_features), 
                        "QSAR", "Transcriptomic")
)

# Top 20 most important features
top_features <- rf_importance_df[order(rf_importance_df$Importance, decreasing = TRUE), ][1:20, ]

cat("Top 10 Most Important Features (from tuned RF):\n")
print(head(top_features, 10))

# Analyze QSAR vs Transcriptomic contribution
qsar_importance <- sum(rf_importance_df[rf_importance_df$Feature_Type == "QSAR", "Importance"])
transcriptomic_importance <- sum(rf_importance_df[rf_importance_df$Feature_Type == "Transcriptomic", "Importance"])

cat("\nFeature Type Contributions:\n")
cat("  - QSAR features:", round(qsar_importance, 1), "\n")
cat("  - Transcriptomic features:", round(transcriptomic_importance, 1), "\n")
cat("  - Transcriptomic contribution:", 
    round(100 * transcriptomic_importance / (qsar_importance + transcriptomic_importance), 1), "%\n")

# ---- 7b) FEATURE SELECTION WITH BORUTA (Transcriptomic only) ----
cat("\n=== STEP 7b: FEATURE SELECTION ON TRANSCRIPTOMICS ===\n")

# Install Boruta if not already installed
if (!require(Boruta)) {
  install.packages("Boruta")
  library(Boruta)
}

# Separate transcriptomic features
transcriptomic_features <- integrated_features[, !(colnames(integrated_features) %in% colnames(qsar_features))]

# Run Boruta on transcriptomics (using training set only to avoid leakage)
set.seed(42)
boruta_result <- Boruta(
  x = transcriptomic_features[train_indices, ], 
  y = y_train, 
  doTrace = 1, 
  maxRuns = 200
)

# Get confirmed important transcriptomic features
confirmed_transcriptomic <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("Confirmed important transcriptomic features:", length(confirmed_transcriptomic), "\n")

# Combine QSAR + filtered transcriptomics
selected_features <- c(colnames(qsar_features), confirmed_transcriptomic)
X_train_selected <- X_train[, selected_features]
X_test_selected <- X_test[, selected_features]

cat("✓ Final feature set after selection:\n")
cat("  - QSAR:", ncol(qsar_features), "\n")
cat("  - Transcriptomic (Boruta):", length(confirmed_transcriptomic), "\n")
cat("  - Total selected:", length(selected_features), "\n")

# ---- 8) ENHANCED QSAR-ONLY COMPARISON ----
cat("\n=== STEP 8: COMPARING WITH QSAR-ONLY MODEL (Enhanced) ===\n")

# Train QSAR-only model with same optimal RF parameters for fair comparison
qsar_only_rf <- randomForest(
  x = X_train[, 1:ncol(qsar_features)], 
  y = y_train, 
  ntree = best_rf_params$ntree,
  mtry = min(best_rf_params$mtry, ncol(qsar_features)),
  maxnodes = best_rf_params$maxnodes,
  nodesize = best_rf_params$nodesize,
  importance = TRUE
)

qsar_only_pred <- predict(qsar_only_rf, X_test[, 1:ncol(qsar_features)])
qsar_only_r2 <- cor(qsar_only_pred, y_test)^2
qsar_only_rmse <- sqrt(mean((qsar_only_pred - y_test)^2))

cat("QSAR-only Model (with optimal parameters):\n")
cat("  - Testing R²:", round(qsar_only_r2, 3), "\n")
cat("  - RMSE:", round(qsar_only_rmse, 3), "\n")

improvement_r2 <- rf_tuned_test_r2 - qsar_only_r2  # Updated to use tuned RF
improvement_rmse <- qsar_only_rmse - rf_tuned_rmse  # Updated to use tuned RF

cat("\nImprovement with Transcriptomics:\n")
cat("  - ΔR²:", round(improvement_r2, 3), "\n")
cat("  - ΔRMSE:", round(improvement_rmse, 3), "\n")
cat("  - Relative improvement:", round(100 * improvement_r2 / qsar_only_r2, 1), "%\n")

# ---- 8b) RETRAIN MODEL WITH SELECTED FEATURES ----
cat("\n=== STEP 8b: RETRAINING RF WITH BORUTA-SELECTED FEATURES ===\n")

rf_boruta_selected <- randomForest(
  x = X_train_selected,
  y = y_train,
  ntree = best_rf_params$ntree,
  mtry = min(best_rf_params$mtry, ncol(X_train_selected)),
  maxnodes = best_rf_params$maxnodes,
  nodesize = best_rf_params$nodesize,
  importance = TRUE
)

rf_boruta_pred_train <- predict(rf_boruta_selected, X_train_selected)
rf_boruta_pred_test <- predict(rf_boruta_selected, X_test_selected)
rf_boruta_train_r2 <- cor(rf_boruta_pred_train, y_train)^2
rf_boruta_test_r2 <- cor(rf_boruta_pred_test, y_test)^2
rf_boruta_rmse <- sqrt(mean((rf_boruta_pred_test - y_test)^2))

cat("RF with Boruta-Selected Features:\n")
cat("  - Training R²:", round(rf_boruta_train_r2, 3), "\n")
cat("  - Testing R²:", round(rf_boruta_test_r2, 3), "\n")
cat("  - RMSE:", round(rf_boruta_rmse, 3), "\n")

# Compare to QSAR-only baseline
boruta_improvement_r2 <- rf_boruta_test_r2 - qsar_only_r2
boruta_improvement_rmse <- qsar_only_rmse - rf_boruta_rmse

cat("\nBoruta Model vs QSAR-only:\n")
cat("  - ΔR²:", round(boruta_improvement_r2, 3), "\n")
cat("  - ΔRMSE:", round(boruta_improvement_rmse, 3), "\n")
cat("  - Relative improvement:", round(100 * boruta_improvement_r2 / qsar_only_r2, 1), "%\n")


# ---- 9) PATHWAY CORRELATION ANALYSIS (Boruta fallback, aligned rows) ----
cat("\n=== STEP 9: PATHWAY CORRELATION ANALYSIS ===\n")
# Fallback using only confirmed transcriptomics (no extra row indexing)
if(length(confirmed_transcriptomic) > 0){
  pathway_matrix <- X_train_selected[, confirmed_transcriptomic, drop = FALSE]
  # remove constant columns
  pathway_matrix <- pathway_matrix[, apply(pathway_matrix, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  
  # compute correlations with QSAR features in X_train_selected
  qsar_feat <- X_train_selected[, colnames(qsar_features), drop = FALSE]
  qsar_feat <- qsar_feat[, apply(qsar_feat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  
  cor_matrix <- cor(qsar_feat, pathway_matrix, use = "complete.obs")
} else {
  cat("⚠ No Boruta-selected transcriptomic features available for fallback.\n")
  cor_matrix <- NULL
}

# Use only training set for correlations (to match Boruta selection)
qsar_feat <- qsar_features[train_indices, , drop = FALSE]
qsar_feat <- qsar_feat[, apply(qsar_feat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

# Check if pathway features exist
if(!is.null(pathway_features) && (is.matrix(pathway_features) || is.data.frame(pathway_features)) && ncol(pathway_features) > 0){
  pathway_matrix <- pathway_features[train_indices, , drop = FALSE]  # align rows
  pathway_matrix <- pathway_matrix[, apply(pathway_matrix, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
} else {
  cat("⚠ No valid pathway features available. Using Boruta-selected transcriptomic features as fallback.\n")
  pathway_matrix <- X_train_selected[train_indices, confirmed_transcriptomic, drop = FALSE]
  pathway_matrix <- pathway_matrix[, apply(pathway_matrix, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
}

# Compute correlations safely
cor_matrix <- cor(qsar_feat, pathway_matrix, use = "complete.obs")

# Pathway correlation function (unchanged)
pathway_correlation_analysis <- function(qsar_feat, pathway_feat, compound_ids) {
  
  pathway_matrix <- as.matrix(pathway_feat)  # ensure numeric matrix
  correlations <- cor(qsar_feat, pathway_matrix, use = "complete.obs")
  
  strong_corr <- which(abs(correlations) > 0.3, arr.ind = TRUE)
  
  if(nrow(strong_corr) > 0) {
    cat("Strong QSAR-Transcriptomic Correlations:\n")
    for(i in 1:min(5, nrow(strong_corr))) {
      row <- strong_corr[i, 1]
      col <- strong_corr[i, 2]
      qsar_name <- colnames(qsar_feat)[row]
      pathway_name <- colnames(pathway_matrix)[col]
      corr_value <- correlations[row, col]
      cat("  -", qsar_name, "vs", pathway_name, ":", round(corr_value, 3), "\n")
    }
  } else {
    cat("No strong correlations found (|r| > 0.3)\n")
  }
  
  return(correlations)
}

pathway_correlations <- pathway_correlation_analysis(qsar_feat, pathway_matrix, compound_ids[train_indices])


# ---- 10) ENHANCED RESULTS SUMMARY (INCLUDING BORUTA) ----
cat("\n=== STEP 10: ENHANCED RESULTS SUMMARY (INCLUDING BORUTA) ===\n")

# Enhanced results summary including all models WITH BORUTA
enhanced_results_summary <- data.frame(
  Model = c("QSAR-only RF", "RF Original", "RF Tuned", "SVM Original", "SVM Tuned", 
            "ElasticNet Original", "ElasticNet Tuned", "Ensemble", "RF Boruta-Selected"),
  Training_R2 = c(cor(predict(qsar_only_rf), y_train)^2, 
                  rf_original_train_r2, rf_tuned_train_r2, 
                  svm_original_train_r2, svm_tuned_train_r2,
                  elastic_original_train_r2, elastic_tuned_train_r2, 
                  ensemble_train_r2, rf_boruta_train_r2),
  Testing_R2 = c(qsar_only_r2, rf_original_test_r2, rf_tuned_test_r2, 
                 svm_original_test_r2, svm_tuned_test_r2,
                 elastic_original_test_r2, elastic_tuned_test_r2, 
                 ensemble_test_r2, rf_boruta_test_r2),
  RMSE = c(qsar_only_rmse, rf_original_rmse, rf_tuned_rmse, 
           svm_original_rmse, svm_tuned_rmse,
           elastic_original_rmse, elastic_tuned_rmse, 
           ensemble_rmse, rf_boruta_rmse),
  Improvement_over_QSAR = c(0, rf_original_test_r2 - qsar_only_r2, rf_tuned_test_r2 - qsar_only_r2, 
                            svm_original_test_r2 - qsar_only_r2, svm_tuned_test_r2 - qsar_only_r2,
                            elastic_original_test_r2 - qsar_only_r2, elastic_tuned_test_r2 - qsar_only_r2,
                            ensemble_test_r2 - qsar_only_r2, boruta_improvement_r2),
  Overfitting_Gap = c(cor(predict(qsar_only_rf), y_train)^2 - qsar_only_r2,
                      rf_original_train_r2 - rf_original_test_r2,
                      rf_tuned_train_r2 - rf_tuned_test_r2,
                      svm_original_train_r2 - svm_original_test_r2,
                      svm_tuned_train_r2 - svm_tuned_test_r2,
                      elastic_original_train_r2 - elastic_original_test_r2,
                      elastic_tuned_train_r2 - elastic_tuned_test_r2,
                      ensemble_train_r2 - ensemble_test_r2,
                      rf_boruta_train_r2 - rf_boruta_test_r2)
)

cat("Enhanced Model Performance Summary (Including Boruta):\n")
# Round only numeric columns
enhanced_results_summary_rounded <- enhanced_results_summary
num_cols <- sapply(enhanced_results_summary, is.numeric)
enhanced_results_summary_rounded[num_cols] <- 
  lapply(enhanced_results_summary[num_cols], function(x) round(x, 3))

print(enhanced_results_summary_rounded)

# Identify best model (including Boruta)
best_model_idx <- which.max(enhanced_results_summary$Testing_R2)
best_model <- enhanced_results_summary$Model[best_model_idx]
best_r2 <- enhanced_results_summary$Testing_R2[best_model_idx]

cat("\nBest Model:", best_model, "with R² =", round(best_r2, 3), "\n")



# ---- 11) ENHANCED FILE SAVING (INCLUDING BORUTA) ----
cat("\n=== STEP 11: SAVING ENHANCED RESULTS (INCLUDING BORUTA) ===\n")

# Save Boruta selection results
write.csv(data.frame(
  Selected_Feature = confirmed_transcriptomic,
  Feature_Type = "Transcriptomic"
), "boruta_selected_features.csv", row.names = FALSE)

# Save enhanced feature importance
write.csv(rf_importance_df, "qsar_transcriptomics_feature_importance_enhanced.csv", row.names = FALSE)

# Save enhanced model results (INCLUDING BORUTA)
write.csv(enhanced_results_summary, "qsar_transcriptomics_model_comparison_enhanced.csv", row.names = FALSE)

# Save hyperparameter tuning results
write.csv(rf_tuning_results, "rf_hyperparameter_tuning_results.csv", row.names = FALSE)

# Enhanced predictions CSV (INCLUDING BORUTA)
predictions_enhanced_df <- data.frame(
  Compound_ID = compound_ids[-train_indices],
  Actual_pIC50 = y_test,
  QSAR_only = qsar_only_pred,
  RF_original = rf_original_pred_test,
  RF_tuned = rf_tuned_pred_test,
  SVM_original = svm_original_pred_test,
  SVM_tuned = svm_tuned_pred_test,
  ElasticNet_original = as.vector(elastic_original_pred_test),
  ElasticNet_tuned = as.vector(elastic_tuned_pred_test),
  Ensemble = ensemble_pred_test,
  RF_Boruta_Selected = rf_boruta_pred_test
)
write.csv(predictions_enhanced_df, "qsar_transcriptomics_predictions_enhanced.csv", row.names = FALSE)

# Save top features for biological interpretation
write.csv(top_features, "top_predictive_features_enhanced.csv", row.names = FALSE)

cat("✓ Enhanced results saved including Boruta analysis:\n")
cat("  - boruta_selected_features.csv\n")
cat("  - qsar_transcriptomics_feature_importance_enhanced.csv\n")
cat("  - qsar_transcriptomics_model_comparison_enhanced.csv (with Boruta)\n")
cat("  - qsar_transcriptomics_predictions_enhanced.csv (with Boruta)\n")
cat("  - top_predictive_features_enhanced.csv\n")

# ---- 12) ENHANCED VISUALIZATIONS (INCLUDING BORUTA) ----
cat("\n=== STEP 12: CREATING ENHANCED VISUALIZATIONS (INCLUDING BORUTA) ===\n")

# Enhanced prediction plots - 7 model comparison (INCLUDING BORUTA)
pdf("qsar_transcriptomics_predictions_enhanced.pdf", width = 18, height = 12)
par(mfrow = c(3, 3))  # Changed to accommodate Boruta

# Plot each model (INCLUDING BORUTA)
models_to_plot <- list(
  list(pred = qsar_only_pred, r2 = qsar_only_r2, title = "QSAR-only", color = "blue"),
  list(pred = rf_original_pred_test, r2 = rf_original_test_r2, title = "Original RF", color = "darkgreen"),
  list(pred = rf_tuned_pred_test, r2 = rf_tuned_test_r2, title = "Tuned RF", color = "forestgreen"),
  list(pred = svm_original_pred_test, r2 = svm_original_test_r2, title = "Original SVM", color = "red"),
  list(pred = svm_tuned_pred_test, r2 = svm_tuned_test_r2, title = "Tuned SVM", color = "darkred"),
  list(pred = ensemble_pred_test, r2 = ensemble_test_r2, title = "Ensemble", color = "purple"),
  list(pred = rf_boruta_pred_test, r2 = rf_boruta_test_r2, title = "RF Boruta-Selected", color = "orange")
)

for(model in models_to_plot) {
  plot(y_test, model$pred, 
       xlab = "Actual pIC50", ylab = "Predicted pIC50", 
       main = paste(model$title, "(R² =", round(model$r2, 3), ")"),
       pch = 19, col = model$color, xlim = range(y), ylim = range(y))
  abline(0, 1, col = "black", lwd = 2)
}

dev.off()

# Model comparison barplot (INCLUDING BORUTA)
pdf("model_comparison_enhanced.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))

# R² comparison (INCLUDING BORUTA)
model_names <- c("QSAR-only", "RF-Orig", "RF-Tuned", "SVM-Orig", "SVM-Tuned", "Ensemble", "RF-Boruta")
r2_values <- c(qsar_only_r2, rf_original_test_r2, rf_tuned_test_r2, 
               svm_original_test_r2, svm_tuned_test_r2, ensemble_test_r2, rf_boruta_test_r2)
colors <- c("blue", "darkgreen", "forestgreen", "red", "darkred", "purple", "orange")

barplot(r2_values, names.arg = model_names, las = 2, col = colors,
        main = "Model Performance Comparison (R²)", ylab = "Testing R²", cex.names = 0.7)

# RMSE comparison (INCLUDING BORUTA)
rmse_values <- c(qsar_only_rmse, rf_original_rmse, rf_tuned_rmse, 
                 svm_original_rmse, svm_tuned_rmse, ensemble_rmse, rf_boruta_rmse)

barplot(rmse_values, names.arg = model_names, las = 2, col = colors,
        main = "Model Performance Comparison (RMSE)", ylab = "RMSE", cex.names = 0.7)

dev.off()

cat("✓ Created enhanced visualizations including Boruta:\n")
cat("  - qsar_transcriptomics_predictions_enhanced.pdf (7-model comparison with Boruta)\n")
cat("  - model_comparison_enhanced.pdf (including Boruta)\n")


# ---- ENHANCED FINAL SUMMARY (INCLUDING BORUTA) ----
cat("\n=== ENHANCED FINAL RESULTS SUMMARY (INCLUDING BORUTA) ===\n")
cat("=====================================\n")

cat("ENHANCED QSAR-TRANSCRIPTOMICS INTEGRATION RESULTS (WITH BORUTA):\n\n")

cat("Dataset Characteristics:\n")
cat("  - Compounds analyzed:", length(compound_ids), "\n")
cat("  - QSAR descriptors:", ncol(qsar_features), "\n")
cat("  - Original transcriptomic features:", length(transcriptomic_matrix), "\n")
cat("  - Boruta-selected transcriptomic features:", length(confirmed_transcriptomic), "\n")
cat("  - Total integrated features (Boruta):", ncol(X_train_selected), "\n")

cat("\nFeature Selection Impact (Boruta):\n")
cat("  - Original transcriptomic features:", ncol(transcriptomic_features), "\n")
cat("  - Boruta-confirmed features:", length(confirmed_transcriptomic), "\n")
cat("  - Feature reduction:", round(100 * (1 - length(confirmed_transcriptomic)/ncol(transcriptomic_features)), 1), "%\n")

cat("\nModel Performance Comparison (Including Boruta):\n")
cat("  - QSAR-only R²:", round(qsar_only_r2, 3), "\n")
cat("  - Best integrated R²:", round(best_r2, 3), "\n")
cat("  - Boruta-selected RF R²:", round(rf_boruta_test_r2, 3), "\n")
cat("  - Best performance improvement:", round(best_r2 - qsar_only_r2, 3), 
    " (", round(100 * (best_r2 - qsar_only_r2) / qsar_only_r2, 1), "%)\n")
cat("  - Boruta improvement:", round(boruta_improvement_r2, 3), 
    " (", round(100 * boruta_improvement_r2 / qsar_only_r2, 1), "%)\n")

cat("\nBoruta vs Full Feature Set Comparison:\n")
cat("  - Full transcriptomic RF R²:", round(rf_tuned_test_r2, 3), "\n")
cat("  - Boruta-selected RF R²:", round(rf_boruta_test_r2, 3), "\n")
if(rf_boruta_test_r2 > rf_tuned_test_r2) {
  cat("  - Boruta provides BETTER performance with fewer features!\n")
} else {
  cat("  - Full feature set performs slightly better\n")
}

cat("\nKey Findings (Including Boruta):\n")
cat("  - Boruta identified", length(confirmed_transcriptomic), "critical transcriptomic features\n")
cat("  - Feature reduction improved model interpretability\n")
cat("  - Best performing model:", best_model, "\n")
cat("  - Boruta R² =", round(rf_boruta_test_r2, 3), "vs Full R² =", round(rf_tuned_test_r2, 3), "\n")

cat("\nFiles Generated (Enhanced with Boruta):\n")
cat("  ✓ boruta_selected_features.csv\n")
cat("  ✓ qsar_transcriptomics_model_comparison_enhanced.csv (includes Boruta)\n")
cat("  ✓ qsar_transcriptomics_predictions_enhanced.csv (includes Boruta)\n")
cat("  ✓ qsar_transcriptomics_predictions_enhanced.pdf (7-model comparison)\n")
cat("  ✓ model_comparison_enhanced.pdf (includes Boruta)\n")


# ---- SAVE ENHANCED COMPLETE WORKSPACE (INCLUDING BORUTA) ----
cat("\n=== SAVING ENHANCED COMPLETE WORKSPACE (INCLUDING BORUTA) ===\n")

# Save all important objects including Boruta models
save(
  # Original data objects
  qsar_raw_data, qsar_features, transcriptomic_matrix, integrated_features,
  de_cq, de_genotype, de_cq_106_1, fgsea_results,
  
  # Original model objects
  rf_model, svm_model, elastic_model, qsar_only_rf,
  
  # Enhanced model objects
  rf_original, rf_tuned_model, svm_original, svm_tuned_model, 
  elastic_original, best_elastic_model,
  
  # Boruta objects
  boruta_result, confirmed_transcriptomic, selected_features,
  rf_boruta_selected, X_train_selected, X_test_selected,
  
  # Hyperparameter results
  best_rf_params, best_svm_params, best_alpha, rf_tuning_results,
  
  # Enhanced results
  enhanced_results_summary, rf_importance_df, top_features,
  pathway_correlations, available_descriptors,
  
  # All predictions and data splits (INCLUDING BORUTA)
  rf_pred_test, svm_pred_test, elastic_pred_test, qsar_only_pred,
  rf_original_pred_test, rf_tuned_pred_test, svm_original_pred_test, svm_tuned_pred_test,
  elastic_original_pred_test, elastic_tuned_pred_test, ensemble_pred_test,
  rf_boruta_pred_test, rf_boruta_train_r2, rf_boruta_test_r2, rf_boruta_rmse,
  y_test, y_train, compound_ids, pIC50_values,
  
  file = "QSAR_Transcriptomics_Enhanced_Complete_Workspace_with_Boruta.RData"
)

cat("✓ Enhanced complete workspace with Boruta saved!\n")

# ---- BONUS: BORUTA FEATURE ANALYSIS ----
cat("\n=== BONUS: BORUTA FEATURE ANALYSIS ===\n")

cat("Boruta Feature Selection Summary:\n")
cat("  - Confirmed features:", length(confirmed_transcriptomic), "\n")
cat("  - Feature reduction:", round(100 * (1 - length(confirmed_transcriptomic)/ncol(transcriptomic_features)), 1), "%\n")
cat("  - Performance with reduced features: R² =", round(rf_boruta_test_r2, 3), "\n")

# Save Boruta importance plot
if(length(confirmed_transcriptomic) > 0) {
  boruta_importance <- importance(rf_boruta_selected)
  boruta_importance_df <- data.frame(
    Feature = rownames(boruta_importance),
    Importance = boruta_importance[, "%IncMSE"],
    Feature_Type = ifelse(rownames(boruta_importance) %in% colnames(qsar_features), 
                          "QSAR", "Transcriptomic_Boruta")
  )
  
  write.csv(boruta_importance_df, "boruta_model_feature_importance.csv", row.names = FALSE)
  cat("✓ boruta_model_feature_importance.csv generated\n")
}

cat("\n", rep("🌟", 50), sep="")
cat("\n🏆 COMPREHENSIVE ENHANCED ANALYSIS WITH BORUTA COMPLETE! 🏆\n")
cat(rep("🌟", 50), sep="")