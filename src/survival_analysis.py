import warnings, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
from pathlib import Path
from sklearn.model_selection import train_test_split
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
warnings.filterwarnings("ignore")

CLINICAL_PATH = Path("data/processed/metabric_merged.csv")
GENOMIC_PATH  = Path("data/processed/metabric_genomic_merged.csv")
PLOT_DIR      = Path("outputs/plots")
DURATION_COL, EVENT_COL, SUBTYPE_COL = "OS_MONTHS","OS_STATUS","CLAUDIN_SUBTYPE"
CLINICAL_FEATURES = ["AGE_AT_DIAGNOSIS","TUMOR_SIZE","TUMOR_STAGE","GRADE",
                     "CHEMOTHERAPY","RADIO_THERAPY","HORMONE_THERAPY",SUBTYPE_COL]
COLOURS = {"LumA":"#1D9E75","LumB":"#378ADD","Her2":"#BA7517",
           "claudin-low":"#D4537E","Basal":"#E24B4A","Normal":"#888780","NC":"#B4B2A9"}

def load(path, features):
    df = pd.read_csv(path)
    pc = [c for c in df.columns if c.startswith("PC")]
    cols = [EVENT_COL,DURATION_COL]+[f for f in features if f in df.columns]+pc
    return df[cols].dropna(subset=[EVENT_COL,DURATION_COL]).query(f"{DURATION_COL}>0")

def onehot(df, col):
    if col in df.columns and df[col].dtype==object:
        d = pd.get_dummies(df[col], prefix=col, drop_first=True)
        df = pd.concat([df.drop(columns=[col]),d], axis=1)
    return df

def plot_km(df):
    fig,ax = plt.subplots(figsize=(10,6))
    res = multivariate_logrank_test(df[DURATION_COL],df[SUBTYPE_COL],df[EVENT_COL])
    for st in sorted(df[SUBTYPE_COL].dropna().unique()):
        m = df[SUBTYPE_COL]==st
        if m.sum()<10: continue
        kmf = KaplanMeierFitter()
        kmf.fit(df.loc[m,DURATION_COL],df.loc[m,EVENT_COL],label=st)
        kmf.plot_survival_function(ax=ax,ci_show=True,color=COLOURS.get(st,"#888780"),linewidth=2)
    ax.set_title(f"KM Survival by PAM50 Subtype  |  log-rank p={res.p_value:.2e}",fontsize=12)
    ax.set_xlabel("Months"); ax.set_ylabel("Survival probability"); ax.set_ylim(0,1.05)
    ax.spines[["top","right"]].set_visible(False)
    fig.savefig(PLOT_DIR/"km_by_subtype.png",dpi=150,bbox_inches="tight"); plt.close(fig)
    print("Saved -> km_by_subtype.png")

def plot_forest(cph):
    s = cph.summary[cph.summary["p"]<0.05].sort_values("exp(coef)")
    fig,ax = plt.subplots(figsize=(9,max(4,len(s)*0.55+1.5)))
    cols = ["#E24B4A" if h>=1 else "#1D9E75" for h in s["exp(coef)"]]
    for i,(_, row) in enumerate(s.iterrows()):
        ax.plot([row["exp(coef) lower 95%"],row["exp(coef) upper 95%"]],[i,i],color=cols[i],lw=1.5)
        ax.scatter(row["exp(coef)"],i,color=cols[i],s=50,zorder=3)
    ax.axvline(1,color="#444441",lw=1,ls="--",alpha=0.6)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels([x.replace("CLAUDIN_SUBTYPE_","Subtype: ").replace("_"," ").title()
                        for x in s.index],fontsize=9)
    ax.set_xlabel("Hazard ratio (95% CI)"); ax.set_title("Cox PH — Significant Predictors (p<0.05)")
    ax.spines[["top","right"]].set_visible(False)
    fig.savefig(PLOT_DIR/"cox_forest_plot.png",dpi=150,bbox_inches="tight"); plt.close(fig)
    print("Saved -> cox_forest_plot.png")

def plot_cindex(c1,c2):
    fig,ax = plt.subplots(figsize=(7,3.5))
    bars = ax.barh(["Clinical only","Clinical + 50 PCs"],[c1,c2],
                   color=["#378ADD","#1D9E75"],height=0.45,alpha=0.85)
    for b,v in zip(bars,[c1,c2]):
        ax.text(v+0.002,b.get_y()+b.get_height()/2,f"{v:.4f}",va="center",fontweight="bold")
    ax.axvline(0.5,color="#888780",lw=1,ls=":",label="Random (0.5)")
    ax.set_xlim(0.45,max(c1,c2)+0.05)
    ax.set_xlabel("C-index"); ax.set_title(f"Model Comparison  |  genomic gain: +{c2-c1:.4f}")
    ax.spines[["top","right"]].set_visible(False); ax.legend(fontsize=8)
    fig.savefig(PLOT_DIR/"cindex_comparison.png",dpi=150,bbox_inches="tight"); plt.close(fig)
    print("Saved -> cindex_comparison.png")

def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    df0 = pd.read_csv(CLINICAL_PATH)
    df0 = df0[df0[DURATION_COL]>0].dropna(subset=[EVENT_COL,DURATION_COL,SUBTYPE_COL])
    plot_km(df0)
    dc = onehot(load(CLINICAL_PATH,CLINICAL_FEATURES),SUBTYPE_COL)
    feats_c = [c for c in dc.columns if c not in [DURATION_COL,EVENT_COL]]
    tr_c,te_c = train_test_split(dc,test_size=0.2,random_state=42)
    cph_c = CoxPHFitter(penalizer=0.1)
    cph_c.fit(tr_c[[DURATION_COL,EVENT_COL]+feats_c],duration_col=DURATION_COL,event_col=EVENT_COL)
    cph_c.print_summary()
    c1 = cph_c.score(te_c,scoring_method="concordance_index")
    print(f"Clinical C-index (test): {c1:.4f}")
    plot_forest(cph_c)
    dg = onehot(load(GENOMIC_PATH,CLINICAL_FEATURES),SUBTYPE_COL)
    feats_g = [c for c in dg.columns if c not in [DURATION_COL,EVENT_COL]]
    tr_g,te_g = train_test_split(dg,test_size=0.2,random_state=42)
    cph_g = CoxPHFitter(penalizer=0.1)
    cph_g.fit(tr_g[[DURATION_COL,EVENT_COL]+feats_g],duration_col=DURATION_COL,event_col=EVENT_COL)
    c2 = cph_g.score(te_g,scoring_method="concordance_index")
    print(f"Genomic C-index (test):  {c2:.4f}  |  gain: +{c2-c1:.4f}")
    plot_cindex(c1,c2)
    print("Done.")

if __name__ == "__main__":
    main()
