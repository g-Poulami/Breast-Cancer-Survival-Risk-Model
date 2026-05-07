import os
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split


def run_survival_analysis(use_genomic=True):
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if use_genomic:
        data_path = os.path.join(BASE_DIR, "data", "processed", "metabric_genomic_merged.csv")
        print("Using genomic + clinical dataset...")
    else:
        data_path = os.path.join(BASE_DIR, "data", "processed", "metabric_merged.csv")
        print("Using clinical-only dataset...")

    try:
        df = pd.read_csv(data_path)
    except FileNotFoundError:
        print(f"Data file not found at {data_path}.")
        print("Run preprocessing.py and genomic_preprocessing.py first.")
        return

    print(f"Loaded {len(df)} patients, {df.shape[1]} columns.\n")

    # -------------------------------------------------------------------------
    # 1. Kaplan-Meier survival curves by Subtype
    # -------------------------------------------------------------------------
    if "SUBTYPE" in df.columns:
        kmf = KaplanMeierFitter()
        fig, ax = plt.subplots(figsize=(10, 6))
        for subtype in sorted(df["SUBTYPE"].dropna().unique()):
            mask = df["SUBTYPE"] == subtype
            kmf.fit(df.loc[mask, "OS_MONTHS"], df.loc[mask, "OS_STATUS"], label=subtype)
            kmf.plot_survival_function(ax=ax)
        ax.set_title("Overall Survival by PAM50 + Claudin-low Subtype (METABRIC)")
        ax.set_xlabel("Time (Months)")
        ax.set_ylabel("Survival Probability")
        plt.tight_layout()
        plt.show()

    # -------------------------------------------------------------------------
    # 2. Cox PH Model — clinical features only (baseline)
    # -------------------------------------------------------------------------
    print("=" * 60)
    print("BASELINE: Clinical-only Cox Model")
    print("=" * 60)

    clinical_features = ["OS_MONTHS", "OS_STATUS", "AGE_AT_DIAGNOSIS", "TUMOR_SIZE",
                         "TUMOR_STAGE", "GRADE", "SUBTYPE", "CHEMOTHERAPY",
                         "HORMONE_THERAPY", "RADIO_THERAPY"]
    avail_clinical = [c for c in clinical_features if c in df.columns]
    df_clinical = df[avail_clinical].dropna()
    df_clinical = pd.get_dummies(df_clinical,
                                  columns=[c for c in ["SUBTYPE", "CHEMOTHERAPY",
                                                        "HORMONE_THERAPY", "RADIO_THERAPY"]
                                           if c in df_clinical.columns],
                                  drop_first=True)

    train_c, test_c = train_test_split(df_clinical, test_size=0.2, random_state=42)
    cph_c = CoxPHFitter(penalizer=0.1)
    cph_c.fit(train_c, duration_col="OS_MONTHS", event_col="OS_STATUS")
    ci_clinical = concordance_index(test_c["OS_MONTHS"],
                                     -cph_c.predict_partial_hazard(test_c),
                                     test_c["OS_STATUS"])
    print(f"Clinical-only C-index: {ci_clinical:.4f}")

    # -------------------------------------------------------------------------
    # 3. Cox PH Model — clinical + genomic PCs
    # -------------------------------------------------------------------------
    if use_genomic:
        print()
        print("=" * 60)
        print("GENOMIC: Clinical + Gene Expression PCs Cox Model")
        print("=" * 60)

        pc_cols = [c for c in df.columns if c.startswith("PC")]
        genomic_features = avail_clinical + pc_cols
        df_genomic = df[genomic_features].dropna()
        df_genomic = pd.get_dummies(df_genomic,
                                     columns=[c for c in ["SUBTYPE", "CHEMOTHERAPY",
                                                           "HORMONE_THERAPY", "RADIO_THERAPY"]
                                              if c in df_genomic.columns],
                                     drop_first=True)

        train_g, test_g = train_test_split(df_genomic, test_size=0.2, random_state=42)
        cph_g = CoxPHFitter(penalizer=0.1)
        cph_g.fit(train_g, duration_col="OS_MONTHS", event_col="OS_STATUS")
        ci_genomic = concordance_index(test_g["OS_MONTHS"],
                                        -cph_g.predict_partial_hazard(test_g),
                                        test_g["OS_STATUS"])
        print(f"Genomic + Clinical C-index: {ci_genomic:.4f}")

        print()
        print("=" * 60)
        print(f"C-index improvement: {ci_genomic - ci_clinical:+.4f}")
        print("=" * 60)


if __name__ == "__main__":
    run_survival_analysis(use_genomic=True)
