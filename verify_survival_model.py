"""
Verify the headline claims of the METABRIC survival risk model.

Purpose
-------
Before featuring this project, confirm (or correct) the three claims a reviewer
will scrutinise:

  1. The C-index improvement from adding genomic PCs is measured on a HELD-OUT
     test set, not in-sample.
  2. That improvement (0.654 -> 0.663, i.e. +0.009) is reported with a
     confidence interval, so it can be judged as real or noise.
  3. The Cox proportional-hazards assumption is checked, and the PAM50 /
     Claudin-low survival separation is supported by a log-rank test.

This script is deliberately defensive: point it at your actual METABRIC data
and it will recompute everything honestly. If you cannot reproduce a claimed
number, change the claim, do not change the data.

Expected input
--------------
A single table (CSV/TSV) with, at minimum:
  - a survival time column (months)
  - an event column (1 = death/event, 0 = censored)
  - clinical covariates (e.g. age, tumour size, grade, stage, nodes)
  - genomic principal components (e.g. PC1..PC50) OR an expression matrix from
    which to compute them
  - a subtype column (e.g. PAM50 / Claudin-low) for the stratified KM

Edit the COLUMN CONFIG block to match your file, then run:
    python verify_survival_model.py --data your_metabric_table.csv
"""

import argparse
import sys
import numpy as np
import pandas as pd

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ----------------------------------------------------------------------------
# COLUMN CONFIG  --  edit these to match your data
# ----------------------------------------------------------------------------
TIME_COL = "OS_MONTHS"
EVENT_COL = "OS_STATUS"
CLINICAL_COLS = ["AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "GRADE", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY"]
SUBTYPE_COL = "SUBTYPE"
# Either provide ready-made PC columns:
PC_COLS = [f"PC{i}" for i in range(1, 51)]
# ...or an expression-matrix prefix to compute PCs from (set PC_COLS = None):
EXPR_PREFIX = None
N_PCS = 50
# ----------------------------------------------------------------------------


def load(path):
    sep = "	" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    yesno = {"YES": 1, "NO": 0, "Yes": 1, "No": 0, "yes": 1, "no": 0}
    for c in ["CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY"]:
        if c in df.columns and df[c].dtype == object:
            df[c] = df[c].map(yesno)
    print(f"Loaded {df.shape[0]} samples x {df.shape[1]} columns")
    return df


def get_pcs(train_df, test_df):
    """Return PC arrays for train/test, fitting PCA on TRAIN ONLY (no leakage)."""
    if PC_COLS and all(c in train_df.columns for c in PC_COLS):
        return train_df[PC_COLS].values, test_df[PC_COLS].values, PC_COLS
    if EXPR_PREFIX is None:
        raise ValueError("Set PC_COLS to existing columns or EXPR_PREFIX to "
                         "compute PCs from an expression matrix.")
    expr_cols = [c for c in train_df.columns if c.startswith(EXPR_PREFIX)]
    scaler = StandardScaler().fit(train_df[expr_cols])
    pca = PCA(n_components=N_PCS, random_state=0).fit(
        scaler.transform(train_df[expr_cols]))
    names = [f"PC{i+1}" for i in range(N_PCS)]
    tr = pca.transform(scaler.transform(train_df[expr_cols]))
    te = pca.transform(scaler.transform(test_df[expr_cols]))
    return tr, te, names


def fit_and_score(train_df, test_df, covariate_cols):
    """Fit Cox on train covariates, return held-out C-index on test."""
    cph = CoxPHFitter(penalizer=0.1)
    cols = covariate_cols + [TIME_COL, EVENT_COL]
    cph.fit(train_df[cols].dropna(), duration_col=TIME_COL, event_col=EVENT_COL)
    # risk = partial hazard; higher = worse. concordance_index expects risk
    # to be NEGATIVELY associated with survival time, so negate.
    risk = cph.predict_partial_hazard(test_df[covariate_cols])
    c = concordance_index(test_df[TIME_COL], -risk, test_df[EVENT_COL])
    return cph, c


def bootstrap_delta(train_df, test_df, clin_cols, pc_names, n_boot=1000, seed=0):
    """Bootstrap the test-set C-index difference (clinical+PC) - (clinical)."""
    rng = np.random.default_rng(seed)
    deltas = []
    test_idx = np.arange(len(test_df))
    for _ in range(n_boot):
        samp = rng.choice(test_idx, size=len(test_idx), replace=True)
        tb = test_df.iloc[samp].reset_index(drop=True)
        _, c_clin = fit_and_score(train_df, tb, clin_cols)
        _, c_full = fit_and_score(train_df, tb, clin_cols + pc_names)
        deltas.append(c_full - c_clin)
    deltas = np.array(deltas)
    return deltas.mean(), np.percentile(deltas, 2.5), np.percentile(deltas, 97.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="METABRIC table (CSV/TSV)")
    ap.add_argument("--test-size", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()

    df = load(args.data)
    needed = [TIME_COL, EVENT_COL]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns {missing}. Edit COLUMN CONFIG.")

    model_cols = [TIME_COL, EVENT_COL] + \
                 [c for c in CLINICAL_COLS if c in df.columns] + \
                 [c for c in PC_COLS if c in df.columns]
    before = len(df)
    df = df.dropna(subset=model_cols)
    print(f'Complete-case filter: kept {len(df)} of {before} samples')
    train_df, test_df = train_test_split(
        df, test_size=args.test_size, random_state=args.seed,
        stratify=df[EVENT_COL])
    print(f"Train {len(train_df)} / Test {len(test_df)} (held out)\n")

    # PCs (fit on train only)
    tr_pc, te_pc, pc_names = get_pcs(train_df, test_df)
    train_df = train_df.copy(); test_df = test_df.copy()
    for i, name in enumerate(pc_names):
        train_df[name] = tr_pc[:, i]
        test_df[name] = te_pc[:, i]

    clin = [c for c in CLINICAL_COLS if c in df.columns]

    # 1 & 2. Held-out C-index, clinical vs clinical+genomic, with bootstrap CI
    _, c_clin = fit_and_score(train_df, test_df, clin)
    _, c_full = fit_and_score(train_df, test_df, clin + pc_names)
    print("HELD-OUT C-index")
    print(f"  clinical only        : {c_clin:.4f}")
    print(f"  clinical + {len(pc_names)} PCs   : {c_full:.4f}")
    print(f"  delta                : {c_full - c_clin:+.4f}")

    mean_d, lo, hi = bootstrap_delta(train_df, test_df, clin, pc_names,
                                     n_boot=args.n_boot, seed=args.seed)
    print(f"  bootstrap delta 95% CI: [{lo:+.4f}, {hi:+.4f}] (mean {mean_d:+.4f})")
    verdict = ("CI excludes 0 -> improvement is supported"
               if lo > 0 else
               "CI includes 0 -> improvement is NOT distinguishable from noise")
    print(f"  verdict              : {verdict}\n")

    # 3a. Proportional hazards assumption check on the full model
    print("PROPORTIONAL HAZARDS CHECK (clinical model)")
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(train_df[clin + [TIME_COL, EVENT_COL]].dropna(),
            duration_col=TIME_COL, event_col=EVENT_COL)
    try:
        cph.check_assumptions(train_df[clin + [TIME_COL, EVENT_COL]].dropna(),
                              show_plots=False)
    except Exception as e:
        print(f"  (assumption check raised: {e})")
    print()

    # 3b. Subtype survival separation with a log-rank test
    if SUBTYPE_COL in df.columns:
        print(f"SUBTYPE SURVIVAL SEPARATION ({SUBTYPE_COL})")
        sub = df.dropna(subset=[SUBTYPE_COL])
        res = multivariate_logrank_test(sub[TIME_COL], sub[SUBTYPE_COL],
                                        sub[EVENT_COL])
        print(f"  multivariate log-rank p = {res.p_value:.3e}")
        kmf = KaplanMeierFitter()
        for name, grp in sub.groupby(SUBTYPE_COL):
            kmf.fit(grp[TIME_COL], grp[EVENT_COL], label=str(name))
            med = kmf.median_survival_time_
            print(f"    {str(name):20s} n={len(grp):4d}  median OS = {med}")
    else:
        print(f"(No {SUBTYPE_COL} column found; skipping subtype log-rank.)")


if __name__ == "__main__":
    main()
