"""
CTGAN and TVAE baseline comparison on K=23 pregnancy coagulation dataset.
Demonstrates that deep generative models fail on very small cohorts.
"""

import pandas as pd
import numpy as np
from sdv.single_table import CTGANSynthesizer, TVAESynthesizer
from sdv.metadata import Metadata
from scipy.stats import ks_2samp, spearmanr
import warnings
import os
import json

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(ROOT, "data")

# ── Load real data ──────────────────────────────────────────────────────────
print("Loading data...")
df_real = pd.read_csv(os.path.join(DATA_PATH, "cleaned_full_data.csv"))
df_sa = pd.read_csv(os.path.join(DATA_PATH, "synthetic_full_longitudinal.csv"))

# Get feature columns (shared between real and SA)
meta_cols = ["SubjectID", "Visit", "DrawTag", "Condition", "PCOS", "DevelopedPE", "DevelopedHTN"]
sa_meta = ["SyntheticID", "Visit"]
feature_cols = [c for c in df_sa.columns if c not in sa_meta]
n_features = len(feature_cols)
print(f"  {len(df_real)} real records, {n_features} features")

# Prepare real data for SDV (features + Visit only)
df_train = df_real[["Visit"] + feature_cols].copy()
df_train = df_train.dropna()
print(f"  {len(df_train)} complete records for training")

# ── Configure metadata ──────────────────────────────────────────────────────
metadata = Metadata.detect_from_dataframe(data=df_train, table_name="patients")

# ── Run CTGAN ───────────────────────────────────────────────────────────────
print("\n=== CTGAN ===")
results = {}

for epochs in [300, 1000, 3000]:
    print(f"\n  Training CTGAN (epochs={epochs})...")
    try:
        ctgan = CTGANSynthesizer(
            metadata,
            epochs=epochs,
            batch_size=min(len(df_train), 500),
            verbose=False,
        )
        ctgan.fit(df_train)
        df_ctgan = ctgan.sample(num_rows=300)  # 100 patients × 3 visits

        # Compute MRE per feature
        mres = []
        for col in feature_cols:
            real_mean = df_real[col].mean()
            synth_mean = df_ctgan[col].mean()
            if abs(real_mean) > 1e-12:
                mres.append(abs(synth_mean - real_mean) / abs(real_mean))

        median_mre = np.median(mres)
        mean_mre = np.mean(mres)
        max_mre = np.max(mres)

        # KS statistics
        ks_stats = []
        for col in feature_cols:
            stat, _ = ks_2samp(df_real[col].dropna(), df_ctgan[col].dropna())
            ks_stats.append(stat)

        median_ks = np.median(ks_stats)

        # Correlation preservation
        real_corr = df_real[feature_cols].corr().values
        synth_corr = df_ctgan[feature_cols].corr().values
        # Replace NaN with 0
        real_corr = np.nan_to_num(real_corr)
        synth_corr = np.nan_to_num(synth_corr)
        mask = np.triu(np.ones_like(real_corr, dtype=bool), k=1)
        corr_mae = np.mean(np.abs(real_corr[mask] - synth_corr[mask]))

        # Check for mode collapse: number of unique values per feature
        n_unique = df_ctgan[feature_cols].nunique().median()

        # Negative values where they shouldn't be (coagulation factors are positive)
        pct_negative = 0
        for col in feature_cols:
            if df_real[col].min() >= 0:
                pct_negative += (df_ctgan[col] < 0).sum()
        pct_negative = pct_negative / (len(df_ctgan) * n_features) * 100

        result = {
            "method": f"CTGAN-{epochs}",
            "median_mre": round(median_mre, 4),
            "mean_mre": round(mean_mre, 4),
            "max_mre": round(max_mre, 4),
            "median_ks": round(median_ks, 4),
            "corr_mae": round(corr_mae, 4),
            "pct_negative": round(pct_negative, 2),
            "median_unique": int(n_unique),
        }
        results[f"CTGAN-{epochs}"] = result
        print(f"    MRE: median={median_mre:.3f}, mean={mean_mre:.3f}, max={max_mre:.3f}")
        print(f"    KS: median={median_ks:.3f}, Corr MAE={corr_mae:.3f}")
        print(f"    Negative values: {pct_negative:.1f}%, Unique values/feature: {int(n_unique)}")

    except Exception as e:
        print(f"    FAILED: {e}")
        results[f"CTGAN-{epochs}"] = {"method": f"CTGAN-{epochs}", "error": str(e)}

# ── Run TVAE ────────────────────────────────────────────────────────────────
print("\n=== TVAE ===")
for epochs in [300, 1000, 3000]:
    print(f"\n  Training TVAE (epochs={epochs})...")
    try:
        tvae = TVAESynthesizer(
            metadata,
            epochs=epochs,
            batch_size=min(len(df_train), 500),
            verbose=False,
        )
        tvae.fit(df_train)
        df_tvae = tvae.sample(num_rows=300)

        mres = []
        for col in feature_cols:
            real_mean = df_real[col].mean()
            synth_mean = df_tvae[col].mean()
            if abs(real_mean) > 1e-12:
                mres.append(abs(synth_mean - real_mean) / abs(real_mean))

        median_mre = np.median(mres)
        mean_mre = np.mean(mres)
        max_mre = np.max(mres)

        ks_stats = []
        for col in feature_cols:
            stat, _ = ks_2samp(df_real[col].dropna(), df_tvae[col].dropna())
            ks_stats.append(stat)
        median_ks = np.median(ks_stats)

        real_corr = df_real[feature_cols].corr().values
        synth_corr = df_tvae[feature_cols].corr().values
        real_corr = np.nan_to_num(real_corr)
        synth_corr = np.nan_to_num(synth_corr)
        mask = np.triu(np.ones_like(real_corr, dtype=bool), k=1)
        corr_mae = np.mean(np.abs(real_corr[mask] - synth_corr[mask]))

        n_unique = df_tvae[feature_cols].nunique().median()

        pct_negative = 0
        for col in feature_cols:
            if df_real[col].min() >= 0:
                pct_negative += (df_tvae[col] < 0).sum()
        pct_negative = pct_negative / (len(df_tvae) * n_features) * 100

        result = {
            "method": f"TVAE-{epochs}",
            "median_mre": round(median_mre, 4),
            "mean_mre": round(mean_mre, 4),
            "max_mre": round(max_mre, 4),
            "median_ks": round(median_ks, 4),
            "corr_mae": round(corr_mae, 4),
            "pct_negative": round(pct_negative, 2),
            "median_unique": int(n_unique),
        }
        results[f"TVAE-{epochs}"] = result
        print(f"    MRE: median={median_mre:.3f}, mean={mean_mre:.3f}, max={max_mre:.3f}")
        print(f"    KS: median={median_ks:.3f}, Corr MAE={corr_mae:.3f}")
        print(f"    Negative values: {pct_negative:.1f}%, Unique values/feature: {int(n_unique)}")

    except Exception as e:
        print(f"    FAILED: {e}")
        results[f"TVAE-{epochs}"] = {"method": f"TVAE-{epochs}", "error": str(e)}

# ── SA comparison values (from our results) ─────────────────────────────────
results["SA"] = {
    "method": "SA",
    "median_mre": 0.019,
    "mean_mre": 0.027,
    "max_mre": 0.157,
    "median_ks": 0.169,
    "corr_mae": 0.200,
    "pct_negative": 0.0,
    "median_unique": 300,
}

results["MVN"] = {
    "method": "MVN",
    "median_mre": 0.024,
    "mean_mre": 0.030,
    "max_mre": 0.145,
    "median_ks": 0.161,
    "corr_mae": 0.090,
    "pct_negative": 0.0,  # approximate
    "median_unique": 300,
}

# ── Save results ────────────────────────────────────────────────────────────
print("\n=== Summary ===")
print(f"{'Method':<15} {'Med MRE':>8} {'Med KS':>8} {'Corr MAE':>9} {'%Neg':>6}")
print("-" * 50)
for key in ["SA", "MVN", "CTGAN-300", "CTGAN-1000", "CTGAN-3000", "TVAE-300", "TVAE-1000", "TVAE-3000"]:
    r = results.get(key, {})
    if "error" in r:
        print(f"{r['method']:<15} FAILED")
    elif "median_mre" in r:
        print(f"{r['method']:<15} {r['median_mre']:>8.3f} {r['median_ks']:>8.3f} {r['corr_mae']:>9.3f} {r['pct_negative']:>6.1f}")

# Save as JSON
with open(os.path.join(DATA_PATH, "ctgan_tvae_comparison.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved results to {os.path.join(DATA_PATH, 'ctgan_tvae_comparison.json')}")
