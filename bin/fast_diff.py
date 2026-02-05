#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import snapatac2 as snap


def parse_args():
    p = argparse.ArgumentParser(description="snapatac2 (v2.7.3) differential peaks with stratified background")
    p.add_argument("--ifile",  dest="infile",     required=True, help="input .h5ad path")
    p.add_argument("--group",  dest="group_col",  required=True, help="column in .obs to group by")
    p.add_argument("--ofile",  dest="outfile",    required=True, help="output prefix (no .csv)")
    p.add_argument("--min_pct", type=float, default=0.01, help="min_pct for diff_test (default 0.01)")
    p.add_argument("--direction", type=str, default="positive",
                   choices=["positive", "negative", "both"], help="direction for diff_test")
    p.add_argument("--workers", type=int, default=4, help="number of worker processes")

    # NEW: stratified sampling controls
    p.add_argument("--target-cells-per-group", type=int, default=1000,
                   help="max cells to retain from the TARGET group (default 1000)")
    p.add_argument("--bg-cells-per-group", type=int, default=1000,
                   help="max cells to sample FROM EACH NON-TARGET group for background (default 1000)")
    p.add_argument("--seed", type=int, default=42, help="random seed base (default 42)")
    return p.parse_args()


def sanitize_for_filename(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9._+-]", "_", s)


def list_groups(infile: Path, group_col: str):
    ad = snap.read(str(infile), backed=None)
    if group_col not in ad.obs.columns:
        raise ValueError(f"Column '{group_col}' not found in .obs. Available: {list(ad.obs.columns)}")
    groups = pd.Index(ad.obs[group_col].unique()).tolist()
    return [g for g in groups if pd.notna(g)]


def _subsample(arr: np.ndarray, n_max: int, rng: np.random.Generator) -> np.ndarray:
    """Return up to n_max unique elements from arr (no replacement)."""
    if arr.size <= n_max:
        return arr
    idx = rng.choice(arr.size, size=n_max, replace=False)
    return arr[np.sort(idx)]


def _make_stratified_background(obs_df: pd.DataFrame,
                                obs_names: np.ndarray,
                                group_col: str,
                                target_group,
                                bg_cells_per_group: int,
                                seed_base: int) -> np.ndarray:
    """
    For each non-target group, sample up to bg_cells_per_group cells;
    then concatenate across all these groups to form background.
    Deterministic via per-(target, other) seeds.
    """
    bg_list = []
    # loop over unique groups except target
    for other in obs_df[group_col].unique():
        if pd.isna(other) or other == target_group:
            continue
        mask_other = (obs_df[group_col] == other).values
        cells_other = obs_names[mask_other]
        if cells_other.size == 0:
            continue
        # per-(target,other) seed → stable & distinct draws
        combo_seed = (hash(f"{target_group}__{other}") + int(seed_base)) % (2**32 - 1)
        rng = np.random.default_rng(combo_seed)
        bg_sample = _subsample(cells_other, bg_cells_per_group, rng)
        if bg_sample.size > 0:
            bg_list.append(bg_sample)

    if len(bg_list) == 0:
        return np.array([], dtype=obs_names.dtype)
    return np.concatenate(bg_list, axis=0)


def run_one_group(infile: str,
                  group_col: str,
                  group_name,
                  min_pct: float,
                  direction: str,
                  out_prefix: str,
                  target_cells_per_group: int,
                  bg_cells_per_group: int,
                  seed: int) -> str:
    """
    Worker: load data, construct TARGET and STRATIFIED BACKGROUND, run diff_test, write CSV.
    """
    ad = snap.read(infile, backed=None)
    obs_df = ad.obs
    obs_names = np.asarray(ad.obs_names)

    # Target cells
    target_mask = (obs_df[group_col] == group_name).values
    target_barcodes = obs_names[target_mask]
    if target_barcodes.size == 0:
        gsafe = sanitize_for_filename(group_name)
        out_path = f"{out_prefix}_{gsafe}.csv"
        pd.DataFrame().to_csv(out_path, index=False)
        return out_path

    # Deterministic target subsample
    tgt_seed = (hash(f"TARGET__{group_name}") + int(seed)) % (2**32 - 1)
    tgt_rng = np.random.default_rng(tgt_seed)
    target_barcodes = _subsample(target_barcodes, target_cells_per_group, tgt_rng)

    # Stratified background: up to bg_cells_per_group from EACH other group
    background_barcodes = _make_stratified_background(
        obs_df=obs_df,
        obs_names=obs_names,
        group_col=group_col,
        target_group=group_name,
        bg_cells_per_group=bg_cells_per_group,
        seed_base=seed,
    )

    if background_barcodes.size == 0:
        gsafe = sanitize_for_filename(group_name)
        out_path = f"{out_prefix}_{gsafe}.csv"
        pd.DataFrame().to_csv(out_path, index=False)
        return out_path

    # snapatac2 diff test
    res = snap.tl.diff_test(
        ad,
        cell_group1=target_barcodes,
        cell_group2=background_barcodes,
        direction=direction,
        min_pct=float(min_pct),
    )

    # Normalize result type to pandas
    try:
        import polars as pl  # noqa: F401
        is_polars = "polars" in type(res).__module__
    except Exception:
        is_polars = False
    if is_polars:
        res = res.to_pandas()
    elif not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)

    # annotate
    res = res.assign(
        **{
            group_col: group_name,
            "n_cells_group1": int(target_barcodes.size),
            "n_cells_group2": int(background_barcodes.size),
        }
    )

    gsafe = sanitize_for_filename(group_name)
    out_path = f"{out_prefix}_{gsafe}.csv"
    res.to_csv(out_path, index=False)
    return out_path


def main():
    args = parse_args()
    in_path = Path(args.infile)
    out_pref = Path(args.outfile)

    if not in_path.exists():
        sys.exit(f"[Error] Input file not found: {in_path}")

    # Keep BLAS threads sane when parallel
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    try:
        groups = list_groups(in_path, args.group_col)
    except Exception as e:
        sys.exit(f"[Error] {e}")
    if len(groups) == 0:
        sys.exit(f"[Error] No groups in '{args.group_col}'.")

    print(f"[info] Groups ({len(groups)}): {groups}")
    print(f"[info] Target cap: {args.target_cells_per_group}  |  Background cap per OTHER group: {args.bg_cells_per_group}")
    print(f"[info] Seed base: {args.seed}")

    per_group_files = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {
            ex.submit(
                run_one_group,
                str(in_path),
                args.group_col,
                g,
                args.min_pct,
                args.direction,
                str(out_pref),
                int(args.target_cells_per_group),
                int(args.bg_cells_per_group),
                int(args.seed),
            ): g for g in groups
        }
        for fut in as_completed(futs):
            g = futs[fut]
            try:
                outp = fut.result()
                per_group_files.append(outp)
                print(f"[done] {g} -> {outp}")
            except Exception as e:
                print(f"[warn] Group {g} failed: {e}")

    # Combine
    dfs = []
    for pth in per_group_files:
        try:
            df = pd.read_csv(pth)
            if df.shape[0] > 0:
                dfs.append(df)
        except Exception as e:
            print(f"[warn] Skipping {pth}: {e}")

    if len(dfs) == 0:
        sys.exit("[Error] No non-empty result files produced.")

    all_results = pd.concat(dfs, axis=0, ignore_index=True)
    combined_csv = out_pref.with_suffix("").parent / f"{out_pref.with_suffix('').name}.csv"
    all_results.to_csv(combined_csv, index=False)
    print(f"[done] Combined results: {combined_csv}")
    print(f"[done] Per-group files prefixed with: {out_pref}_<group>.csv")


if __name__ == "__main__":
    main()
