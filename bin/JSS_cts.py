#!/usr/bin/env python3

"""
Call cell-type-specific peaks from:
  (1) CPM matrix: rows = peaks (row index), cols = cell types
  (2) LR CSV with columns:
      'feature name','log2(fold_change)','p-value','adjusted p-value',
      'subclass_corrected','n_cells_group1','n_cells_group2'

Workflow:
  - Keep LR-significant peaks (padj/log2FC thresholds)
  - For each (peak, target group) compute:
      * JSS vs one-hot target (1 - JSD)
      * top share p_(1), second share p_(2), ratio p_(1)/p_(2)
  - Apply exclusivity filters (JSS / top / second / ratio)
  - Output a table with metrics + is_specific boolean

Note:
  - Matrix values should be non-negative CPMs.
  - Each peak's CPM vector is normalized to a probability distribution across cell types.
  - No prevalence thresholds are used here.
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


def parse_args():
    p = argparse.ArgumentParser(description="Cell-type-specific peaks via LR + JSS (CPM matrix, LR schema fixed)")
    p.add_argument("--matrix", required=True, help="CPM matrix CSV (rows=peaks as index, cols=cell types)")
    p.add_argument("--lr", required=True, nargs="+", help="LR CSV file(s) or globs")
    p.add_argument("--out", required=True, help="Output CSV")

    # Thresholds for LR filtering
    p.add_argument("--lr-max-padj", type=float, default=0.05, help="Max adjusted p-value (default: 0.05)")
    p.add_argument("--lr-min-log2fc", type=float, default=1.0, help="Min log2 fold-change (default: 1.0)")

    # JSS & exclusivity thresholds
    p.add_argument("--min-jss", type=float, default=0.90, help="Min JSS (1 - JSD) (default: 0.90)")
    p.add_argument("--min-top-share", type=float, default=0.65, help="Min top share p1 (default: 0.65)")
    p.add_argument("--max-second-share", type=float, default=0.20, help="Max second share p2 (default: 0.20)")
    p.add_argument("--min-ratio", type=float, default=2.5, help="Min p1/p2 ratio (default: 2.5)")

    # Column names in the LR CSV (defaults match your schema)
    p.add_argument("--lr-peak-col", default="feature name")
    p.add_argument("--lr-log2fc-col", default="log2(fold_change)")
    p.add_argument("--lr-padj-col", default="adjusted p-value")
    p.add_argument("--lr-group-col", default="subclass_corrected")

    p.add_argument("--epsilon", type=float, default=1e-12, help="Small value to stabilize probabilities")
    return p.parse_args()


def load_matrix(path):
    # Peaks as row index
    df = pd.read_csv(path, index_col=0)
    # ensure numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if df.shape[0] == 0 or df.shape[1] == 0:
        sys.exit("[Error] Matrix appears empty after parsing.")
    if (df.values < 0).any():
        sys.exit("[Error] Matrix contains negative values; expected non-negative CPMs.")
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
        missing = [c for c in [peak_col, group_col, padj_col, log2fc_col] if c not in t.columns]
        if missing:
            sys.exit(f"[Error] Missing LR columns in {p}: {missing}")
        t = t[[peak_col, group_col, padj_col, log2fc_col]].copy()
        t = t.rename(columns={
            peak_col: "peak_id",
            group_col: "group",
            padj_col: "padj",
            log2fc_col: "log2fc",
        })
        t = t[(t["padj"] <= max_padj) & (t["log2fc"] >= min_log2fc)]
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


def compute_js_similarity_to_one_hot(p_vec, target_index, base=2):
    """
    Jensen–Shannon-based similarity:
      JSD = JS divergence, JS_distance = sqrt(JSD)
      JSS := 1 - JSD
    """
    k = p_vec.size
    q = np.zeros(k, dtype=float)
    q[target_index] = 1.0
    js_dist = jensenshannon(p_vec, q, base=base)  # sqrt divergence
    jsd = js_dist ** 2
    jss = 1.0 - jsd
    return jss, jsd, js_dist


def main():
    args = parse_args()

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

    # Keep LR peaks present in matrix and LR groups present as matrix columns
    lr = lr[lr["peak_id"].isin(mat.index)]
    lr = lr[lr["group"].isin(mat.columns)]
    if lr.empty:
        sys.exit("[Error] After filtering to matrix rows/columns, no LR rows remain.")

    celltypes = list(mat.columns)
    ct_to_idx = {ct: i for i, ct in enumerate(celltypes)}
    peak_to_row = {pk: i for i, pk in enumerate(mat.index)}

    mat_np = mat.values  # (n_peaks x n_celltypes)

    out_rows = []
    for peak, grp, padj, l2fc in lr[["peak_id", "group", "padj", "log2fc"]].itertuples(index=False):
        ridx = peak_to_row[peak]
        vec = mat_np[ridx, :].copy()

        # Normalize CPM row to a probability distribution across cell types
        p = row_normalize(vec, eps=args.epsilon)

        # Exclusivity components
        order = np.argsort(p)[::-1]
        p1 = p[order[0]]
        p2 = p[order[1]] if p.size > 1 else 0.0
        ratio = p1 / (p2 + args.epsilon)

        # JSS vs one-hot at target group
        tgt_idx = ct_to_idx[grp]
        jss, jsd, js_dist = compute_js_similarity_to_one_hot(p, tgt_idx, base=2)

        # Decision
        pass_jss = (jss >= args.min_jss)
        pass_top = (p1 >= args.min_top_share)
        pass_leak = (p2 <= args.max_second_share)
        pass_ratio = (ratio >= args.min_ratio)
        is_specific = (pass_jss and pass_top and pass_leak and pass_ratio)

        out_rows.append({
            "peak_id": peak,
            "target_group": grp,
            "padj": padj,
            "log2fc": l2fc,
            "JSS": jss,
            "JSD": jsd,
            "JS_distance": js_dist,
            "top_share": p1,
            "second_share": p2,
            "top_second_ratio": ratio,
            "pass_JSS": pass_jss,
            "pass_top_share": pass_top,
            "pass_second_share": pass_leak,
            "pass_ratio": pass_ratio,
            "is_specific": is_specific,
        })

    out_df = pd.DataFrame(out_rows)
    if out_df.empty:
        sys.exit("[Error] No rows computed. Check inputs/filters.")

    # Optional ranking
    with np.errstate(divide="ignore", invalid="ignore"):
        rank_score = out_df["JSS"] * np.log2(out_df["top_second_ratio"] + 1e-12) * (out_df["top_share"] - out_df["second_share"])
    out_df["specificity_score"] = rank_score.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    out_path = Path(args.out)
    out_df.to_csv(out_path, index=False)
    print(f"[done] Wrote {out_path} with {out_df.shape[0]} rows.")
    print("[hint] Use is_specific == True for your specific peak set; sort by specificity_score to rank.")


if __name__ == "__main__":
    main()

