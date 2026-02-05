#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Empirical p-value & FDR for JSS (1 - JSD) specificity, à la 'fdrtool':
- Read CPM matrix (rows=peaks index, cols=cell types).
- Read LR-significant results (your schema).
- For each (peak, target group): compute JSS vs one-hot.
- Standardize JSS to z-scores (across peaks per target group).
- Fit 2-component GaussianMixture on z to approximate null/alt.
- Compute one-sided p-value under null component and BH q-value.
- Also report local-FDR-like metric: posterior prob(null | z).

Usage example:
python jss_empirical_fdr.py \
  --matrix cpm_matrix.csv \
  --lr lr_results.csv \
  --out jss_empirical.csv \
  --lr-max-padj 0.05 --lr-min-log2fc 1.0 \
  --min-jss 0.90 --min-top-share 0.65 --max-second-share 0.20 --min-ratio 2.5
"""

import argparse
import glob
import sys
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(description="JSS empirical p/q via Gaussian mixture on z-scores")
    p.add_argument("--matrix", required=True, help="CPM matrix CSV (rows=peaks index, cols=cell types)")
    p.add_argument("--lr", required=True, nargs="+", help="LR CSV(s) or glob(s)")
    p.add_argument("--out", required=True, help="Output CSV")

    # LR schema (matches your columns)
    p.add_argument("--lr-peak-col", default="feature name")
    p.add_argument("--lr-log2fc-col", default="log2(fold_change)")
    p.add_argument("--lr-pval-col", default="p-value")              # unused, but kept for completeness
    p.add_argument("--lr-padj-col", default="adjusted p-value")
    p.add_argument("--lr-group-col", default="subclass_corrected")

    # LR thresholds
    p.add_argument("--lr-max-padj", type=float, default=0.05)
    p.add_argument("--lr-min-log2fc", type=float, default=1.0)

    # JSS & exclusivity thresholds (optional)
    p.add_argument("--min-jss", type=float, default=0.90)
    p.add_argument("--min-top-share", type=float, default=0.65)
    p.add_argument("--max-second-share", type=float, default=0.20)
    p.add_argument("--min-ratio", type=float, default=2.5)

    # Mixture options
    p.add_argument("--components", type=int, default=2, help="GMM components (default 2)")
    p.add_argument("--random-state", type=int, default=17)
    p.add_argument("--epsilon", type=float, default=1e-12)

    p.add_argument("--mode", choices=["enrich", "deplete"], default="enrich",
                   help="JSS target: enrich=one-hot at group (default), deplete=zero at group, uniform over others")
    return p.parse_args()


def load_matrix(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if df.shape[0] == 0 or df.shape[1] == 0:
        sys.exit("[Error] Matrix empty after parsing.")
    if (df.values < 0).any():
        sys.exit("[Error] Negative values found; expected non-negative CPMs.")
    return df


def load_lr_tables(patterns, peak_col, group_col, padj_col, log2fc_col, max_padj, min_log2fc):
    paths = []
    for pat in patterns:
        paths.extend(glob.glob(pat))
    if not paths:
        sys.exit(f"[Error] No LR files matched: {patterns}")
    dfs = []
    for p in paths:
        t = pd.read_csv(p)
        need = [peak_col, group_col, padj_col, log2fc_col]
        miss = [c for c in need if c not in t.columns]
        if miss:
            sys.exit(f"[Error] Missing LR columns in {p}: {miss}")
        t = t[[peak_col, group_col, padj_col, log2fc_col]].copy()
        t = t.rename(columns={
            peak_col: "peak_id",
            group_col: "group",
            padj_col: "padj",
            log2fc_col: "log2fc",
        })
        t = t[(t["padj"] <= max_padj) & (np.abs(t["log2fc"]) >= min_log2fc)]
        dfs.append(t)
    if not dfs:
        sys.exit("[Error] No LR rows passed thresholds.")
    lr = pd.concat(dfs, axis=0, ignore_index=True).drop_duplicates()
    return lr


def row_normalize(v, eps=1e-12):
    v = np.asarray(v, dtype=float)
    v = np.maximum(v, 0.0)
    s = v.sum()
    if s <= 0:
        v = v + eps
        s = v.sum()
    return v / s


def compute_jss_target(p_vec, tgt_idx, mode="enrich", base=2):
    """
    mode='enrich'  -> q is one-hot at tgt_idx (what you originally had)
    mode='deplete' -> q[tgt_idx]=0, q[others]=1/(C-1)  (anti-k target)
    """
    C = p_vec.size
    if mode == "enrich":
        q = np.zeros_like(p_vec, dtype=float)
        q[tgt_idx] = 1.0
    elif mode == "deplete":
        q = np.full(C, 1.0 / max(C - 1, 1), dtype=float)
        q[tgt_idx] = 0.0
    else:
        raise ValueError(f"Unknown mode: {mode}")

    jsd_sqrt = jensenshannon(p_vec, q, base=base)
    jsd = jsd_sqrt ** 2
    jss = 1.0 - jsd
    return jss, jsd, jsd_sqrt


def exclusivity_checks(p, min_top, max_second, min_ratio, eps):
    order = np.argsort(p)[::-1]
    p1 = p[order[0]]
    p2 = p[order[1]] if p.size > 1 else 0.0
    ratio = p1 / (p2 + eps)
    return p1, p2, ratio, (p1 >= min_top) and (p2 <= max_second) and (ratio >= min_ratio)


def gmm_empirical_pq(z, n_components=2, random_state=17):
    """
    Fit a Gaussian mixture on z-scores and extract:
      - which component is 'null' (the one with the smaller mean),
      - local FDR proxy = posterior P(null | z),
      - one-sided p-values under null N(mu0, sigma0),
      - BH q-values (FDR) across all rows.

    Returns: dict with arrays for post_null, pval, qval, (mu0, sd0).
    """
    z = np.asarray(z, dtype=float).reshape(-1, 1)
    gmm = GaussianMixture(n_components=n_components, random_state=random_state)
    gmm.fit(z)

    means = gmm.means_.ravel()
    sds = np.sqrt(gmm.covariances_.ravel())
    weights = gmm.weights_.ravel()

    # Define null as the component with the smaller mean (background near 0)
    null_idx = np.argmin(means)
    mu0, sd0, w0 = means[null_idx], sds[null_idx], weights[null_idx]

    # Posterior responsibilities for each component
    resp = gmm.predict_proba(z)  # (n, k)
    post_null = resp[:, null_idx]  # local FDR proxy

    # One-sided p-value for "large positive" z (specificity): p = 1 - CDF_null(z)
    # If you prefer two-sided, use 2*(1 - CDF(|z|)).
    pval = 1.0 - norm.cdf(z.ravel(), loc=mu0, scale=sd0)

    # BH FDR across all p-values
    qval = multipletests(pval, method="fdr_bh")[1]

    return {
        "post_null": post_null,
        "pval": pval,
        "qval": qval,
        "mu0": mu0,
        "sd0": sd0,
        "w0": w0,
        "means": means,
        "sds": sds,
        "weights": weights,
    }


def main():
    args = parse_args()

    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    
    mat = load_matrix(args.matrix)
    lr = load_lr_tables(
        args.lr,
        args.lr_peak_col,
        args.lr_group_col,
        args.lr_padj_col,
        args.lr_log2fc_col,
        args.lr_max_padj,
        args.lr_min_log2fc,
    )

    # filter to overlapping peaks and valid groups
    lr = lr[lr["peak_id"].isin(mat.index)]
    lr = lr[lr["group"].isin(mat.columns)]
    if lr.empty:
        sys.exit("[Error] After filtering to matrix rows/columns, no LR rows remain.")

    celltypes = list(mat.columns)
    ct_to_idx = {ct: i for i, ct in enumerate(celltypes)}
    peak_to_row = {pk: i for i, pk in enumerate(mat.index)}
    M = mat.values

    rows = []
    # Compute JSS & exclusivity per (peak, group)
    for peak, grp, padj, l2fc in lr[["peak_id", "group", "padj", "log2fc"]].itertuples(index=False):
        ridx = peak_to_row[peak]
        vec = M[ridx, :].copy()
        p = row_normalize(vec, eps=args.epsilon)

        jss, jsd, jsd_sqrt = compute_jss_target(p, ct_to_idx[grp], mode=args.mode, base=2)
        p1, p2, ratio, pass_excl = exclusivity_checks(
            p, args.min_top_share, args.max_second_share, args.min_ratio, args.epsilon
        )
        pass_jss = (jss >= args.min_jss)

        rows.append({
            "peak_id": peak,
            "target_group": grp,
            "padj": padj,
            "log2fc": l2fc,
            "JSS": jss,
            "JSD": jsd,
            "JS_distance": jsd_sqrt,
            "top_share": p1,
            "second_share": p2,
            "top_second_ratio": ratio,
            "pass_JSS": pass_jss,
            "pass_exclusivity": pass_excl,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("[Error] No rows computed.")

    # --- Empirical inference per target group (like the R code does by attr) ---
    out_blocks = []
    for grp, sub in df.groupby("target_group", sort=False):
        # z-score JSS within the group
        z = (sub["JSS"] - sub["JSS"].mean()) / (sub["JSS"].std(ddof=0) + 1e-12)
        sub = sub.copy()
        sub["JSS_z"] = z

        # Fit GMM on z and compute p/q + local FDR proxy
        g = gmm_empirical_pq(
            z.values,
            n_components=args.components,
            random_state=args.random_state
        )
        sub["p_empirical"] = g["pval"]
        sub["q_empirical"] = g["qval"]
        sub["lfdr_gmm"] = g["post_null"]  # posterior prob of null (local FDR proxy)
        sub["null_mu"] = g["mu0"]
        sub["null_sd"] = g["sd0"]
        sub["null_weight"] = g["w0"]
        out_blocks.append(sub)

    out = pd.concat(out_blocks, axis=0, ignore_index=True)

    # Final boolean flags
    out["is_specific_score"] = out["pass_JSS"] & out["pass_exclusivity"]
    # Optionally also require empirical significance (e.g., q <= 0.05)
    out["is_specific_empirical"] = out["is_specific_score"] & (out["q_empirical"] <= 0.05)

    out_path = Path(args.out)
    out.to_csv(out_path, index=False)
    print(f"[done] Wrote {out_path} with {out.shape[0]} rows.")
    print("Columns added: JSS_z, p_empirical, q_empirical, lfdr_gmm (posterior null prob), null_mu, null_sd.")
    print("Use is_specific_empirical for a conservative set (score filters + empirical FDR).")


if __name__ == "__main__":
    main()


