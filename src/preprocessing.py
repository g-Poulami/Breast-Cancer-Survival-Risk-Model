import pandas as pd
from pathlib import Path

RAW_CLINICAL = Path("data/raw/brca_metabric_clinical_data.tsv")
OUT_PATH     = Path("data/processed/metabric_merged.csv")

FEATURES = [
    "PATIENT_ID","OS_STATUS","OS_MONTHS","AGE_AT_DIAGNOSIS",
    "CLAUDIN_SUBTYPE","TUMOR_SIZE","TUMOR_STAGE","GRADE",
    "CHEMOTHERAPY","RADIO_THERAPY","HORMONE_THERAPY",
]

def load_clinical(path):
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    print(f"Loaded {len(df):,} rows from {path.name}")
    return df

def select_and_clean(df):
    available = [c for c in FEATURES if c in df.columns]
    df = df[available].copy()
    before = len(df)
    df.dropna(subset=["OS_STATUS","OS_MONTHS"], inplace=True)
    print(f"Dropped {before - len(df)} rows missing OS fields")
    df["OS_STATUS"] = df["OS_STATUS"].apply(
        lambda x: 1 if str(x).startswith("1") or str(x).upper()=="DECEASED" else 0)
    print(f"Cohort: {len(df):,} patients | Events: {df['OS_STATUS'].sum():,}")
    return df

def main():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df = load_clinical(RAW_CLINICAL)
    df = select_and_clean(df)
    df.to_csv(OUT_PATH, index=False)
    print(f"Saved -> {OUT_PATH}")

if __name__ == "__main__":
    main()
