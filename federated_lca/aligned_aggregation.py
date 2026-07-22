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
import fed_log2_lib as fl 
NODES = os.path.join(HERE, "nodes")
OUT = os.path.join(HERE, "output")
TAB = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)
K = 5


def theme():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
    })


def _rowz(M):
    return (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True, ddof=0)


def evaluate(cent, G, gt_labels):
    col, _ = fl.hungarian_align(cent, G, metric="corr")
    per, overall = fl.evaluate_vs_gt(cent, G, col, None, gt_labels)
    Fz, Gz = _rowz(cent), _rowz(G[col])
    per["rmse_shape"] = np.sqrt(((Fz - Gz) ** 2).mean(1))
    null = fl.permutation_null(cent, G)
    p = float((null <= overall).mean())
    return col, per, overall, p


def main():
    theme()
    gt = pd.read_csv(os.path.join(TAB, "gt_raw_normalized.csv"), index_col=0)
    ABX = list(gt.columns); gt_labels = list(gt.index); G = gt.to_numpy()

    raw, Z, meta, mean, sd = fl.collect_log2(NODES, ABX, mode="pooled")
    w = fl.collect_weights(NODES, meta)
    print(f"profiles {Z.shape[0]}  countries {meta['country'].nunique()}  "
          f"isolate-weight range [{w.min():.0f}, {w.max():.0f}]")

    _, base_cent = fl.central_consensus_cluster(Z, K=K, method="kmeans")
    _, base_per, base_rmse, base_p = evaluate(base_cent, G, gt_labels)
    cent, assign, obj = fl.aligned_aggregate(Z, w, meta, K=K)
    meta["global_class"] = [f"Fed{k+1}" for k in assign]
    col, per, rmse, p = evaluate(cent, G, gt_labels)

    fed = pd.DataFrame(cent, index=[f"Fed{k+1}" for k in range(K)], columns=ABX)
    fed.to_csv(os.path.join(TAB, "federated_global_classes_aligned.csv"))
    per.to_csv(os.path.join(TAB, "gt_match_aligned.csv"), index=False)
    meta.to_csv(os.path.join(TAB, "country_class_assignment_aligned.csv"), index=False)

    print("\nAligned federated and GT (per matched class):")
    print(per.round(3).to_string(index=False))
    print(f"\nMEAN r = {per['pearson_r'].mean():.4f}   RMSE = {rmse:.4f}   perm p = {p:.4f}")

    comp = pd.DataFrame([
        {"method": "k-means, equal weight",
         "mean_r_vs_gt": base_per["pearson_r"].mean(), "rmse": base_rmse, "perm_p": base_p},
        {"method": "aligned and isolate-weighted",
         "mean_r_vs_gt": per["pearson_r"].mean(), "rmse": rmse, "perm_p": p},
    ])
    comp.to_csv(os.path.join(TAB, "aligned_vs_baseline.csv"), index=False)

    matched = gt.iloc[col]


if __name__ == "__main__":
    main()
