import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def load_expression_data(mrna_path, n_pcs=50):
    """
    Loads mRNA expression data, filters low-variance genes,
    and reduces to top principal components.
    Returns a DataFrame with PATIENT_ID + PC columns.
    """
    print("  Loading mRNA expression data...")
    df = pd.read_csv(mrna_path, sep="\t", index_col=0)

    # Drop Entrez_Gene_Id column if present
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(columns=["Entrez_Gene_Id"])

    print(f"  Expression matrix: {df.shape[0]} genes x {df.shape[1]} patients")

    # Transpose: patients as rows, genes as columns
    df = df.T
    df.index.name = "PATIENT_ID"

    # Drop genes with any missing values
    df = df.dropna(axis=1)
    print(f"  After dropping NaN genes: {df.shape[1]} genes remaining")

    # Filter to top 5000 most variable genes to reduce noise
    gene_var = df.var(axis=0)
    top_genes = gene_var.nlargest(5000).index
    df = df[top_genes]
    print(f"  Using top 5000 most variable genes")

    # Standardise before PCA
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df)

    # PCA
    n_pcs = min(n_pcs, X_scaled.shape[1], X_scaled.shape[0])
    pca = PCA(n_components=n_pcs, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    explained = pca.explained_variance_ratio_.cumsum()[-1]
    print(f"  PCA: {n_pcs} PCs explain {explained*100:.1f}% of variance")

    pc_cols = [f"PC{i+1}" for i in range(n_pcs)]
    pca_df = pd.DataFrame(X_pca, index=df.index, columns=pc_cols)
    pca_df = pca_df.reset_index()  # PATIENT_ID becomes a column
    return pca_df


def merge_clinical_and_expression(clinical_path, mrna_path, output_path, n_pcs=50):
    """Merges cleaned clinical data with PCA-reduced expression data."""
    clinical_df = pd.read_csv(clinical_path)
    print(f"  Clinical data: {len(clinical_df)} patients")

    pca_df = load_expression_data(mrna_path, n_pcs=n_pcs)
    print(f"  Expression PCA data: {len(pca_df)} patients")

    merged = clinical_df.merge(pca_df, on="PATIENT_ID", how="inner")
    print(f"  Merged dataset: {len(merged)} patients")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    merged.to_csv(output_path, index=False)
    print(f"  Saved to {output_path}")
    return merged


if __name__ == "__main__":
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    clinical_path = os.path.join(BASE_DIR, "data", "processed", "metabric_merged.csv")
    mrna_path = os.path.join(BASE_DIR, "data", "raw", "data_mrna_illumina_microarray.txt")
    output_path = os.path.join(BASE_DIR, "data", "processed", "metabric_genomic_merged.csv")

    print("Starting genomic preprocessing...")
    merged = merge_clinical_and_expression(clinical_path, mrna_path, output_path, n_pcs=50)
    print(f"\nFinal dataset shape: {merged.shape}")
    print(f"Columns: clinical + {sum(1 for c in merged.columns if c.startswith('PC'))} PCs")
