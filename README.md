# 🎗️ Breast Cancer Survival & Risk Modelling

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=flat-square&logo=python)
![Lifelines](https://img.shields.io/badge/Lifelines-0.27%2B-green?style=flat-square)
![Scikit-Learn](https://img.shields.io/badge/Scikit--Learn-1.0%2B-orange?style=flat-square&logo=scikit-learn)
![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)

A clinical and genomic survival risk stratification model built on the **METABRIC breast cancer dataset** (~2,000 patients). This project combines clinical features with gene expression data to predict overall survival using Kaplan-Meier estimators and Cox Proportional Hazards models.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Project Structure](#project-structure)
- [Setup & Installation](#setup--installation)
- [Pipeline](#pipeline)
- [Results](#results)
- [Plots](#plots)
- [Key Findings](#key-findings)
- [Future Work](#future-work)
- [References](#references)

---

## Overview

Breast cancer is one of the most heterogeneous malignancies, with survival outcomes varying widely across molecular subtypes. This project builds a survival risk model that:

- Stratifies patients by **PAM50 + Claudin-low molecular subtype** using Kaplan-Meier curves
- Quantifies the effect of clinical and treatment variables on survival using a **Cox Proportional Hazards model**
- Integrates **20,000+ gene expression features** via PCA to assess the added predictive value of genomics over clinical data alone
- Evaluates model performance using the **concordance index (C-index)** on a held-out test set

---

## Dataset

**METABRIC (Molecular Taxonomy of Breast Cancer International Consortium)**
- Source: [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric)
- 2,509 primary breast cancer patients (Nature 2012 & Nat Commun 2016)
- After preprocessing: **1,980 patients** with complete clinical + genomic data
- Events (deceased): **1,144 / 1,981** patients (57.7% event rate)

### Data Files (not included in repo — download from cBioPortal)

| File | Description |
|---|---|
| `data/raw/brca_metabric_clinical_data.tsv` | Clinical, treatment & survival data |
| `data/raw/data_mrna_illumina_microarray.txt` | mRNA expression (Illumina microarray, 20,603 genes) |

> **Note:** Raw data files exceed GitHub's 100MB limit and are excluded from this repository. Download directly from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) or via the [METABRIC data bundle](https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz).

---

## Project Structure

```
Breast-Cancer-Survival-Risk-Model/
├── src/
│   ├── preprocessing.py            # Clinical data cleaning & feature engineering
│   ├── genomic_preprocessing.py    # mRNA expression PCA pipeline
│   └── survival_analysis.py        # KM curves, Cox model, C-index evaluation
├── data/
│   ├── raw/                        # Raw input files (not tracked in git)
│   └── processed/                  # Processed CSVs (not tracked in git)
├── outputs/
│   └── plots/                      # Generated figures
│       ├── km_by_subtype.png
│       ├── cox_forest_plot.png
│       └── cindex_comparison.png
├── notebooks/
│   └── analysis.ipynb              # Exploratory analysis notebook
├── requirements.txt
└── README.md
```

---

## Setup & Installation

### 1. Clone the repository
```bash
git clone https://github.com/g-Poulami/Breast-Cancer-Survival-Risk-Model.git
cd Breast-Cancer-Survival-Risk-Model
```

### 2. Install dependencies
```bash
pip install -r requirements.txt
```

### 3. Download the data
Download the METABRIC bundle from cBioPortal:
```bash
mkdir -p data/raw
# Download clinical data and mRNA expression from:
# https://www.cbioportal.org/study/summary?id=brca_metabric
# Place files in data/raw/
```

Or download the full bundle:
```bash
cd data/raw
wget https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz
tar -xzf brca_metabric.tar.gz
cp brca_metabric/data_mrna_illumina_microarray.txt .
```

### 4. Run the pipeline
```bash
# Step 1: Clean clinical data
python src/preprocessing.py

# Step 2: Process gene expression data (PCA)
python src/genomic_preprocessing.py

# Step 3: Run survival analysis
python src/survival_analysis.py
```

---

## Pipeline

### Step 1 — Clinical Preprocessing (`src/preprocessing.py`)
- Loads the raw METABRIC TSV, skipping `#` metadata rows
- Selects 11 clinically relevant features: survival status/time, age, subtype, tumour size/stage, grade, and treatment indicators
- Drops rows missing survival-critical fields (`OS_STATUS`, `OS_MONTHS`) only — preserving maximum sample size
- Converts `OS_STATUS` to binary (1 = Deceased, 0 = Living/Censored)
- Resolves paths relative to `__file__` for reproducibility across environments
- **Output:** `data/processed/metabric_merged.csv` (1,981 patients, 11 columns)

### Step 2 — Genomic Preprocessing (`src/genomic_preprocessing.py`)
- Loads the Illumina microarray expression matrix (20,603 genes × 1,980 patients)
- Transposes to patient × gene format and drops genes with missing values (20,592 genes retained)
- Filters to the **top 5,000 most variable genes** to reduce noise
- Applies `StandardScaler` then **PCA (50 components)**, explaining **64.6% of variance**
- Merges PCA scores with clinical data on `PATIENT_ID`
- **Output:** `data/processed/metabric_genomic_merged.csv` (1,980 patients, 61 columns)

### Step 3 — Survival Analysis (`src/survival_analysis.py`)
- **Kaplan-Meier estimator** stratified by PAM50 + Claudin-low molecular subtype
- **Cox Proportional Hazards model** with `penalizer=0.1` to handle correlated features
- Categorical variables (subtype, treatments) one-hot encoded with `drop_first=True`
- **80/20 train-test split** (random state 42) for unbiased C-index evaluation
- Two models compared: clinical-only vs clinical + 50 genomic PCs

---

## Results

### Model Performance (Held-out Test Set)

| Model | Features | C-index |
|---|---|---|
| Clinical Only | Age, tumour size/stage, grade, subtype, treatments | 0.6536 |
| Clinical + Genomics | Clinical + 50 gene expression PCs | **0.6634** |
| Improvement | +50 genomic PCs (64.6% variance explained) | **+0.0098** |

> A C-index of 0.5 = random prediction; 1.0 = perfect. Values of 0.65–0.75 are considered good for clinical survival models.

### Cox Model — Significant Predictors (p < 0.05)

| Feature | Hazard Ratio | Interpretation |
|---|---|---|
| Age at Diagnosis | 1.03 | 3% increased risk per year of age |
| Tumour Size | 1.01 | Larger tumours associated with worse survival |
| Tumour Stage | 1.38 | Higher stage = 38% increased hazard |
| Grade | 1.14 | Higher grade = 14% increased hazard |
| Chemotherapy | 1.32 | Reflects selection bias (sicker patients treated) |
| Radiotherapy | 0.86 | Associated with 14% reduced hazard |

> **Concordance (training set): 0.67** | **Partial AIC: 10065.42**

---

## Plots

### Kaplan-Meier Survival Curves by Molecular Subtype
![KM by Subtype](outputs/plots/km_by_subtype.png)

Patients stratified by PAM50 + Claudin-low subtype show clearly divergent survival trajectories. LumA (Luminal A) patients have the best prognosis, while Basal-like and claudin-low subtypes show substantially worse survival, consistent with published literature.

---

### Cox Proportional Hazards — Forest Plot
![Cox Forest Plot](outputs/plots/cox_forest_plot.png)

Hazard ratios with 95% confidence intervals for all clinical features. Features with confidence intervals crossing 1.0 (the vertical dashed line) are not statistically significant. Age, tumour stage, and grade are the strongest significant predictors of mortality.

---

### Model Comparison — C-index
![C-index Comparison](outputs/plots/cindex_comparison.png)

Adding 50 gene expression principal components to the clinical model improves the held-out C-index from 0.6536 to 0.6634 (+0.0098). The modest improvement suggests clinical features already capture the majority of survival signal in this dataset, with genomics providing incremental gains.

---

## Key Findings

1. **Molecular subtype is the strongest prognostic stratifier** — LumA patients have substantially better survival than Basal or claudin-low subtypes, confirming the clinical utility of PAM50 classification.

2. **Clinical features dominate survival prediction** — Age at diagnosis, tumour stage, and grade are the most statistically significant predictors in the Cox model (all p < 0.005).

3. **Gene expression adds modest but real value** — Integrating 50 PCs from 20,000+ genes improves the C-index by ~1%, suggesting diminishing returns from raw expression data without targeted feature selection.

4. **Chemotherapy shows a positive hazard ratio** — This reflects indication bias: chemotherapy is preferentially given to higher-risk patients, so its apparent association with worse outcomes is confounded by disease severity, not a true causal effect.

5. **Radiotherapy shows a protective association** — HR = 0.86, consistent with its established role in reducing local recurrence and improving survival in early-stage breast cancer.

---

## Future Work

- [ ] **Random Survival Forest** — capture non-linear interactions between features
- [ ] **Targeted gene selection** — use known prognostic genes (BRCA1, TP53, ESR1, ERBB2) instead of PCA
- [ ] **Add mutation data** — incorporate somatic mutation burden from `data_mutations.txt`
- [ ] **Cross-validation** — replace single train/test split with 5-fold CV for more robust C-index estimation
- [ ] **Deep survival models** — DeepSurv or DRSA for end-to-end genomic survival modelling
- [ ] **Calibration analysis** — assess whether predicted survival probabilities are well-calibrated

---

## References

1. Curtis C. et al. *The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.* Nature, 2012.
2. Rueda O.M. et al. *Dynamics of breast-cancer relapse reveal late-recurring ER-positive genomic subgroups.* Nature, 2019.
3. Cerami E. et al. *The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data.* Cancer Discovery, 2012.
4. Davidson-Pilon C. *lifelines: survival analysis in Python.* Journal of Open Source Software, 2019.
5. Cox D.R. *Regression Models and Life-Tables.* Journal of the Royal Statistical Society, 1972.
