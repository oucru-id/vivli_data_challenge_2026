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


def main():
    theme()
    gt = pd.read_csv(os.path.join(TAB, "gt_trajectory.csv"))
    gt_pct = ft.to_matrix(gt, "year", "class", "pct", classes=ft.CLASSES)
    gt_n = ft.to_matrix(gt, "year", "class", "n", classes=ft.CLASSES)
    gt_tot = gt_n.sum(axis=1)

    best = pd.read_csv(os.path.join(TAB, "best_method.csv"))["best_method"][0]
    mapping = pd.read_csv(os.path.join(TAB, "alignment_map_best.csv"))
    df, weights, _ = ft.collect_trajectories(NODES, min_cell=0)
    glob = ft.aggregate_global(df, mapping)
    fed_pct = ft.to_matrix(glob, "Year", "global_class", "pct")
    fed_n = ft.to_matrix(glob, "Year", "global_class", "N")
    fed_tot = fed_n.sum(axis=1)
    match = ft.hungarian_match_trajectories(fed_pct, gt_pct)
    print(f"using Method {best}; class map: {match}")

    years = ft.YEARS.astype(int)

    rows = []
    for f_cls, g_cls in match.items():
        a_f, b_f = ft.logistic_trend(years, fed_n[f_cls], fed_tot)
        lo_f, hi_f = ft.boot_trend_ci(years, fed_n[f_cls], fed_tot)
        a_g, b_g = ft.logistic_trend(years, gt_n[g_cls], gt_tot)
        lo_g, hi_g = ft.boot_trend_ci(years, gt_n[g_cls], gt_tot)
        rows.append({"gt_class": g_cls, "fed_class": f_cls,
                     "beta_fed": b_f, "fed_lo": lo_f, "fed_hi": hi_f,
                     "beta_gt": b_g, "gt_lo": lo_g, "gt_hi": hi_g,
                     "abs_diff": abs(b_f - b_g)})
    tr = pd.DataFrame(rows).sort_values("gt_class").reset_index(drop=True)
    tr.to_csv(os.path.join(TAB, "trend_coefficients.csv"), index=False)
    print(tr.round(4).to_string(index=False))
    print(f"mean |beta_fed - beta_gt| = {tr['abs_diff'].mean():.4f}   "
          f"corr(beta) = {np.corrcoef(tr['beta_fed'], tr['beta_gt'])[0,1]:.3f}")

    f1 = [k for k, v in match.items() if v == "Class1"][0]
    e_fed = ft.emergence_year(years, fed_pct[f1].to_numpy())
    e_gt = ft.emergence_year(years, gt_pct["Class1"].to_numpy())
    cp_fed = ft.changepoint(years, fed_pct[f1].to_numpy())
    cp_gt = ft.changepoint(years, gt_pct["Class1"].to_numpy())
    print("\n=== Class1 emergence ===")
    print(f"threshold crossings  federated: {e_fed}")
    print(f"threshold crossings  GT       : {e_gt}")
    print(f"changepoint federated: {cp_fed}   GT: {cp_gt}")
    pd.DataFrame([{"source": "federated", **{f"cross_{k}pct": v for k, v in e_fed.items()},
                   "changepoint": cp_fed},
                  {"source": "GT", **{f"cross_{k}pct": v for k, v in e_gt.items()},
                   "changepoint": cp_gt}]).to_csv(
        os.path.join(TAB, "class1_emergence.csv"), index=False)

if __name__ == "__main__":
    main()
