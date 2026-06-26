# Breast cancer survival risk stratification (METABRIC)

Survival modelling of the METABRIC breast cancer cohort (1,980 patients),
combining clinical covariates, molecular subtype, and gene expression to ask one
question: **how much does the molecular profile of a tumour add to survival
prediction beyond clinical variables alone, and which subtypes drive the
difference?**

All numbers below are computed on a held-out test set and reproduced by
`verify_survival_model.py`. Cohort after complete-case filtering on the model
covariates: 1,875 patients (1,312 train / 563 held-out test).

## Headline findings

**1. Gene expression adds real, statistically supported prognostic value beyond
clinical variables.** On the held-out test set, a clinical Cox model achieves a
concordance index of **0.611**; adding 50 genomic principal components raises it
to **0.642**, a gain of **+0.031** (bootstrap 95% CI **[+0.008, +0.054]**). The
confidence interval excludes zero, so the genomic contribution is distinguishable
from noise, though it is modest in absolute terms.

**2. Molecular subtype separates survival strongly.** Kaplan-Meier survival
differs markedly across PAM50 / Claudin-low subtypes (multivariate log-rank
**p = 8.4e-11**). The poorest-prognosis groups are Her2 and Luminal B; Luminal A
and Claudin-low fare best.

**Median overall survival by subtype**

| Subtype | n | Median OS (months) |
| --- | --- | --- |
| claudin-low | 201 | 220.9 |
| LumA | 667 | 188.5 |
| Normal | 135 | 153.9 |
| Basal | 203 | 137.9 |
| LumB | 453 | 123.0 |
| Her2 | 211 | 104.0 |

*(An "NC" / not-classified group of 5 patients is omitted as too small to
interpret.)*

## How to read the C-index result

The concordance index measures how often the model ranks a higher-risk patient
as having the shorter survival; 0.5 is random, 1.0 is perfect. The headline is
not the point estimate but the bootstrap confidence interval on the
clinical-versus-genomic difference: because it excludes zero
([+0.008, +0.054]), the 50 genomic PCs carry prognostic information that clinical
variables do not. The effect is real but modest, which is the honest and
expected result for bulk expression PCs over an already-informative clinical
baseline.

## Methods

- **Cohort:** METABRIC (1,980 patients) with overall survival, clinical
  covariates (age at diagnosis, tumour size, grade, chemotherapy, hormone
  therapy, radiotherapy), PAM50 / Claudin-low subtype, and microarray gene
  expression. Complete-case analysis on the model covariates retains 1,875
  patients. `TUMOR_STAGE` is excluded as a covariate because it is missing for
  26% of patients (514/1,980); dropping those patients would discard a quarter
  of the cohort with likely non-random missingness, so the covariate is dropped
  instead, with `TUMOR_SIZE` and `GRADE` retained.
- **Dimensionality reduction:** PCA over the expression matrix, fit on the
  training split only (no test-set leakage), retaining 50 components.
- **Survival models:** Cox proportional hazards, clinical and clinical +
  genomic, with L2 penalisation.
- **Subtype analysis:** Kaplan-Meier estimators per subtype with a multivariate
  log-rank test.
- **Evaluation:** concordance index on a 30% held-out test split, with a
  1,000-iteration bootstrap confidence interval on the clinical-versus-genomic
  difference.

## Reproduce

```bash
pip install -r requirements.txt
python verify_survival_model.py --data data/processed/metabric_genomic_merged.csv
```

`verify_survival_model.py` recomputes the held-out C-index for both models, the
bootstrap CI on their difference, the proportional-hazards diagnostic, and the
subtype log-rank with medians. It is the single source of truth for every number
quoted here.

## Limitations

- **Modest genomic effect.** The held-out gain from 50 expression PCs is +0.031
  C-index (CI [+0.008, +0.054]): statistically supported but small. Clinical
  variables and subtype carry most of the prognostic signal in this cohort.
- **Proportional-hazards violations.** The PH assumption is violated for
  `AGE_AT_DIAGNOSIS`, `GRADE`, `CHEMOTHERAPY`, and `HORMONE_THERAPY` (Schoenfeld
  test p < 0.005 each); `TUMOR_SIZE` and `RADIO_THERAPY` are consistent with PH.
  With ~1,300 observations the test is highly sensitive, so these flags are
  expected, but a more rigorous model would stratify on grade and the binary
  treatment variables and consider a non-linear or time-interaction term for age.
- **Retrospective cohort.** METABRIC has its own ascertainment and treatment-era
  characteristics; this is a prognostic stratification exercise, not a validated
  clinical tool.
- **PCs are not interpretable as biology.** They capture expression variance, not
  specific pathways.

## References

- Curtis C, et al. The genomic and transcriptomic architecture of 2,000 breast
  tumours reveals novel subgroups. *Nature* 486, 346-352 (2012).
- Pereira B, et al. The somatic mutation profiles of 2,433 breast cancers refine
  their genomic and transcriptomic landscapes. *Nat Commun* 7, 11479 (2016).
- Harrell FE, et al. Evaluating the yield of medical tests. *JAMA* 247,
  2543-2546 (1982). (concordance index)
