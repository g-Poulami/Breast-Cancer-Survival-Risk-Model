import pandas as pd
import os

def clean_clinical_data(raw_path):
    """Cleans the raw METABRIC clinical file and skips metadata rows."""
    # skip the comment rows starting with '#'
    df = pd.read_csv(raw_path, sep='\t', comment='#')

    # METABRIC clinical columns often use mixed case or underscores[cite: 4]
    cols = ['Patient ID', 'Overall Survival Status', 'Overall Survival (Months)', 'Age at Diagnosis', 'Claudin Subtype', 'Tumor Size']

    rename_map = {
        'Patient ID': 'PATIENT_ID',
        'Overall Survival Status': 'OS_STATUS',
        'Overall Survival (Months)': 'OS_MONTHS',
        'Age at Diagnosis': 'AGE_AT_DIAGNOSIS',
        'Claudin Subtype': 'CLAUDIN_SUBTYPE',
        'Tumor Size': 'TUMOR_SIZE'
    }

    # Select only if they exist to prevent KeyError, then rename[cite: 4]
    existing_cols = [c for c in cols if c in df.columns]
    df = df[existing_cols].rename(columns=rename_map)

    df = df.dropna()

    # Convert status to binary: 1 for Deceased, 0 for Living[cite: 4]
    if 'OS_STATUS' in df.columns:
        df['OS_STATUS'] = df['OS_STATUS'].apply(lambda x: 1 if 'DECEASED' in str(x).upper() else 0)

    return df

if __name__ == "__main__":
    # Pathing relative to the src directory
    clinical_raw = "../data/raw/brca_metabric_clinical_data.tsv"
    output_path = "../data/processed/metabric_merged.csv"

    print("Starting data preprocessing...")

    if os.path.exists(clinical_raw):
        clinical_df = clean_clinical_data(clinical_raw)
        clinical_df.to_csv(output_path, index=False)
        print(f"Success! Processed data saved to {output_path}[cite: 1, 5]")
    else:
        print(f"File not found: {clinical_raw}. Please ensure data is in data/raw/")
