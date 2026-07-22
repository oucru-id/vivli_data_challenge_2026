#!/usr/bin/env python3

from __future__ import annotations

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fed_traj_lib as ft 

NODES = os.path.join(HERE, "nodes")
OUT = os.path.join(HERE, "output")
TAB = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)
MIN_CELL = 0


def theme():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25,
    })


def int_years(ax):
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(round(x))}"))


def run_method(df, mapping, gt_pct, tag):
    glob = ft.aggregate_global(df, mapping)
    fed_pct = ft.to_matrix(glob, "Year", "global_class", "pct")
    match = ft.hungarian_match_trajectories(fed_pct, gt_pct)
    per, overall = ft.evaluate_vs_gt(fed_pct, gt_pct, match)
    null = ft.permutation_null(fed_pct, gt_pct)
    p = float((null <= overall).mean())
    print(f"\n--- {tag} ---")
    print("mapping:", {k: v for k, v in match.items()})
    print(per.round(3).to_string(index=False))
    print(f"overall MAE = {overall:.3f} pp   mean r = {per['pearson_r'].mean():.3f}   "
          f"perm p = {p:.4f}")
    return glob, fed_pct, match, per, overall, p


def main():
    theme()
    gt = pd.read_csv(os.path.join(TAB, "gt_trajectory.csv"))
    gt_pct = ft.to_matrix(gt, "year", "class", "pct", classes=ft.CLASSES)

    df, weights, supp = ft.collect_trajectories(NODES, min_cell=MIN_CELL)
    print(f"countries {df['Country'].nunique()}  rows {len(df)}  total N {int(df['N'].sum()):,}"
          f"   (GT total {int(gt['n'].sum()):,})")

    mapA = ft.align_from_mic()
    chk = mapA.groupby("Country")["global_class"].agg(lambda s: s.nunique() == len(s))
    print(f"Method A mapping 1-to-1 within country: {int(chk.sum())}/{len(chk)}")
    globA, fedA, matchA, perA, maeA, pA = run_method(df, mapA, gt_pct, "A · MIC-profile alignment")

    mapB, objB = ft.align_from_curves(df, weights)
    chkB = mapB.groupby("Country")["global_class"].agg(lambda s: s.nunique() == len(s))
    print(f"\nMethod B mapping 1-to-1 within country: {int(chkB.sum())}/{len(chkB)}  "
          f"(internal objective {objB:.4f})")
    globB, fedB, matchB, perB, maeB, pB = run_method(df, mapB, gt_pct,
                                                     "B · trajectory-native alignment")

    comp = pd.DataFrame([
        {"method": "A · MIC-profile alignment", "overall_mae_pp": maeA,
         "mean_r": perA["pearson_r"].mean(), "perm_p": pA},
        {"method": "B · trajectory-native alignment", "overall_mae_pp": maeB,
         "mean_r": perB["pearson_r"].mean(), "perm_p": pB},
    ])
    print(comp.round(4).to_string(index=False))
    best = "A" if maeA <= maeB else "B"

    mapA.to_csv(os.path.join(TAB, "alignment_map_A_mic.csv"), index=False)
    mapB.to_csv(os.path.join(TAB, "alignment_map_B_curves.csv"), index=False)
    globA.to_csv(os.path.join(TAB, "global_trajectory_A.csv"), index=False)
    globB.to_csv(os.path.join(TAB, "global_trajectory_B.csv"), index=False)
    perA.to_csv(os.path.join(TAB, "gt_metrics_A.csv"), index=False)
    perB.to_csv(os.path.join(TAB, "gt_metrics_B.csv"), index=False)
    comp.to_csv(os.path.join(TAB, "alignment_head_to_head.csv"), index=False)
    pd.DataFrame([supp]).to_csv(os.path.join(TAB, "suppression_summary.csv"), index=False)

    fed_best, match_best, glob_best = ((fedA, matchA, globA) if best == "A"
                                       else (fedB, matchB, globB))
    inv = {v: k for k, v in match_best.items()}
    pd.DataFrame([{"best_method": best, **{f"gt_{g}": f for g, f in inv.items()}}]).to_csv(
        os.path.join(TAB, "best_method.csv"), index=False)
    (mapA if best == "A" else mapB).to_csv(os.path.join(TAB, "alignment_map_best.csv"), index=False)

if __name__ == "__main__":
    main()
