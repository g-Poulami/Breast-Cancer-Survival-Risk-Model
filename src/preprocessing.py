import pandas as pd
import os


def clean_clinical_data(raw_path):
    df = pd.read_csv(raw_path, sep="\t", comment="#")
    print(f"  Raw data shape: {df.shape}")
    cols = [
        "Patient ID", "Overall Survival Status", "Overall Survival (Months)",
        "Age at Diagnosis", "Pam50 + Claudin-low subtype", "Tumor Size",
        "Tumor Stage", "Chemotherapy", "Hormone Therapy", "Radio Therapy",
        "Neoplasm Histologic Grade"
    ]
    rename_map = {
        "Patient ID": "PATIENT_ID", "Overall Survival Status": "OS_STATUS",
        "Overall Survival (Months)": "OS_MONTHS", "Age at Diagnosis": "AGE_AT_DIAGNOSIS",
        "Pam50 + Claudin-low subtype": "SUBTYPE", "Tumor Size": "TUMOR_SIZE",
        "Tumor Stage": "TUMOR_STAGE", "Chemotherapy": "CHEMOTHERAPY",
        "Hormone Therapy": "HORMONE_THERAPY", "Radio Therapy": "RADIO_THERAPY",
        "Neoplasm Histologic Grade": "GRADE"
    }
    existing_cols = [c for c in cols if c in df.columns]
    df = df[existing_cols].rename(columns=rename_map)
    before = len(df)
    df = df.dropna(subset=["OS_STATUS", "OS_MONTHS"])
    print(f"  Dropped {before - len(df)} rows missing OS_STATUS/OS_MONTHS.")
    df["OS_STATUS"] = df["OS_STATUS"].apply(lambda x: 1 if "DECEASED" in str(x).upper() else 0)
    for col in ["OS_MONTHS", "AGE_AT_DIAGNOSIS", "TUMOR_SIZE"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"  Retained {len(df)} patients. Events: {df['OS_STATUS'].sum()}")
    return df


if __name__ == "__main__":
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    clinical_raw = os.path.join(BASE_DIR, "data", "raw", "brca_metabric_clinical_data.tsv")
    output_path = os.path.join(BASE_DIR, "data", "processed", "metabric_merged.csv")
    print("Starting data preprocessing...")
    if os.path.exists(clinical_raw):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        clinical_df = clean_clinical_data(clinical_raw)
        clinical_df.to_csv(output_path, index=False)
        print(f"Success! Saved to {output_path}")
    else:
        print(f"Error: File not found at {clinical_raw}")
