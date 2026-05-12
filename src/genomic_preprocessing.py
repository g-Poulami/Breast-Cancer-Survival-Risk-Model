import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

EXPR_PATH     = Path("data/raw/data_mrna_illumina_microarray.txt")
CLINICAL_PATH = Path("data/processed/metabric_merged.csv")
OUT_PATH      = Path("data/processed/metabric_genomic_merged.csv")
N_TOP_GENES   = 5000
N_PCS         = 50

def main():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    print("Loading expression matrix...")
    expr = pd.read_csv(EXPR_PATH, sep="\t", index_col=0)
    expr = expr.select_dtypes(include=[np.number]).T
    expr.index.name = "PATIENT_ID"
    expr.dropna(axis=1, inplace=True)
    print(f"  {expr.shape[1]:,} genes retained after dropping NaNs")
    top_genes = expr.std(axis=0).nlargest(N_TOP_GENES).index
    expr = expr[top_genes]
    print(f"  Selected top {N_TOP_GENES:,} most-variable genes")
    X = StandardScaler().fit_transform(expr)
    pca = PCA(n_components=N_PCS, random_state=42)
    X_pca = pca.fit_transform(X)
    print(f"  {N_PCS} PCs explain {pca.explained_variance_ratio_.sum()*100:.1f}% variance")
    pc_df = pd.DataFrame(X_pca, index=expr.index,
                         columns=[f"PC{i+1}" for i in range(N_PCS)])
    clinical = pd.read_csv(CLINICAL_PATH)
    merged = clinical.merge(pc_df, on="PATIENT_ID", how="inner")
    print(f"Merged cohort: {len(merged):,} patients, {merged.shape[1]} columns")
    merged.to_csv(OUT_PATH, index=False)
    print(f"Saved -> {OUT_PATH}")

if __name__ == "__main__":
    main()
