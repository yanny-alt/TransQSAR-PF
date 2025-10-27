# TransQSAR-pf: Bio-Informed QSAR for Antimalarial Drug Discovery


**Integrating *Plasmodium falciparum* transcriptomic stress signatures with classical QSAR achieves 6.1% improved prediction accuracy with 98.3% feature reduction**

## 🎯 Overview

Traditional QSAR models predict drug activity from molecular structure alone, ignoring the biological state of target organisms. TransQSAR-pf addresses this limitation by incorporating transcriptomic stress response signatures from chloroquine-resistant *Plasmodium falciparum* strains.

**Key Innovation:** Biology-guided feature selection identifies conserved stress pathways as universal predictors of compound efficacy.

## 📊 Performance at a Glance

| Model | Testing R² | RMSE | Features | Improvement |
|-------|------------|------|----------|-------------|
| QSAR-only (baseline) | 0.719 | 0.529 | 15 | - |
| Full Integration | 0.602 | 0.649 | 779 | -16.3% ❌ (overfitting) |
| **TransQSAR-pf (Boruta)** | **0.762** | **0.470** | **28** | **+6.1% ✅** |

## 🔬 Scientific Abstract

Public microarray data (GSE10022) from chloroquine-resistant *P. falciparum* strains were analyzed using limma differential expression and fgsea pathway enrichment. We engineered 764 transcriptomic features (expression signatures, pathway scores, variability metrics) and integrated them with 15 QSAR descriptors from 125 triazolopyrimidine derivatives.

Boruta feature selection reduced the feature space to 13 critical transcriptomic predictors representing:

- Drug response signatures (40%)
- Genotype-specific effects (30%)
- Strain variability (30%)

Machine learning models (Random Forest, SVM, Elastic Net) were optimized via 5-fold cross-validation with hyperparameter tuning. The Boruta-selected Random Forest achieved R² = 0.762 (6.1% improvement over QSAR-only baseline).

**Biological Mapping:** 71.2% of predictive importance derived from conserved unknown-function genes, identifying high-priority targets for mechanistic validation.


## 🧬 Methodology

### 1️⃣ Transcriptomic Analysis
- **Data:** GSE10022 (24,563 probes, 18 RNA samples, 3 *P. falciparum* genotypes)
- **Method:** Limma moderated t-tests with Benjamini-Hochberg FDR correction
- **Output:** Differential expression for CQ treatment, genotype effects

### 2️⃣ Pathway Enrichment (GSEA)
- **Algorithm:** fgsea (fast preranked GSEA)
- **Database:** PlasmoDB GO annotations (18,000 annotations, 2,500 terms)
- **Significant pathways (p < 0.05):**
  - Conserved *Plasmodium* proteins (p = 0.005, NES = 1.68)
  - RNA-binding proteins (p = 0.020, NES = -1.52)
  - PfEMP1 virulence factors (p = 0.036, NES = 1.45)

### 3️⃣ Feature Engineering
- **Transcriptomic Features (764 total):**
  - Differential expression signatures: 600
  - Pathway enrichment scores: 3
  - Expression variability: 100
  - Functional group profiles: 61
- **QSAR Descriptors (15):**
  - Key descriptors from Nnadi et al.: slogP, MW, HBD, HBA, TPSA, nRotB, vsurf_W2, vsurf_CW2, npr1, pmi3

### 4️⃣ Feature Selection
- **Algorithm:** Boruta (Random Forest-based importance with shadow features)
- **Iterations:** 200
- **Result:** 764 → 13 features (98.3% reduction)
- **Validation:** 5-fold cross-validation

### 5️⃣ Machine Learning
- **Models tested:**
  - Random Forest (grid search: ntree, mtry, maxnodes, nodesize)
  - SVM with RBF kernel (grid search: C, γ)
  - Elastic Net (α tuning: 0.1-0.9)
  - Weighted ensemble
- **Best model:** Boruta-selected Random Forest
  - Training R² = 0.899 (minimal overfitting)
  - Testing R² = 0.762
  - RMSE = 0.470

## 📊 Performance at a Glance

| Model | Testing R² | RMSE | Features | Improvement |
|-------|------------|------|----------|-------------|
| QSAR-only (baseline) | 0.719 | 0.529 | 15 | - |
| Full Integration | 0.602 | 0.649 | 779 | -16.3% ❌ (overfitting) |
| **TransQSAR-pf (Boruta)** | **0.762** | **0.470** | **28** | **+6.1% ✅** |

## 📊 Key Results

### Model Performance
[View Model Performance PDF](https://github.com/yanny-alt/TransQSAR-PF/blob/main/figures/model_performance/model_comparison_enhanced.pdf)

### Top 5 Predictive Features (Boruta Selection)
| Rank | Feature | Importance | Biological Function |
|------|---------|------------|---------------------|
| 1 | CQ_106_1_DE_40 | 9.43 | Direct chloroquine response in wild-type strain |
| 2 | CQ_DE_128 | 5.31 | Cross-strain drug response signature |
| 3 | Variability_Pf.12.198 | 4.61 | Conserved unknown protein (71.2% importance) |
| 4 | CQ_DE_169 | 4.54 | Metabolic adaptation to drug pressure |
| 5 | Genotype_DE_88 | 4.17 | Genetic background effect (17.7% importance) |

### Pathway Enrichment
[View GSEA Enrichment Plot](https://github.com/yanny-alt/TransQSAR-PF/blob/main/figures/gsea/gsea_correct_plot.pdf)

**Significant Findings:**
- PfEMP1 virulence factors (p = 0.036): Linked to compound efficacy
- Conserved unknown proteins (p = 0.005): Strongest predictors → drug target candidates
- RNA-binding proteins (p = 0.020): Post-transcriptional regulation under drug stress

## 🎯 Biological Insights

**Novel Drug Targets Identified**
- Conserved unknown-function genes contribute 71.2% of predictive importance
  - Represent unexplored essential pathways
  - Prioritized for functional validation

**PfEMP1 virulence pathway correlates with drug susceptibility**
- Surface antigen variation affects compound uptake/efficacy
- Therapeutic strategy: target antigenic switching machinery

**Genotype-specific responses account for 17.7% of variance**
- Drug resistance mutations alter transcriptional landscape
- Personalized treatment strategies based on strain genotype

## Implications for Drug Discovery
- ✅ **Virtual Screening:** Use TransQSAR-pf scores to prioritize compounds
- ✅ **Mechanism Prediction:** Transcriptomic signatures reveal mode of action
- ✅ **Resistance Profiling:** Genotype features predict resistance likelihood
- ✅ **Target Identification:** High-importance genes = candidate drug targets

## 📈 Comparison with Literature

| Study | Method | Dataset | R² | Improvement |
|-------|---------|---------|----|-------------|
| Nnadi et al. (2025) | Classical QSAR | 125 Tpz | 0.67 | Baseline |
| **TransQSAR-pf (this work)** | **QSAR + Transcriptomics** | **125 Tpz + GSE10022** | **0.762** | **+13.7%** |

**Advantage:** Incorporates biological context without additional in vitro assays
