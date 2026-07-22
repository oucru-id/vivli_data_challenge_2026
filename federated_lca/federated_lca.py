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
os.makedirs(TAB, exist_ok=True)
K_GLOBAL = 5
SCALED_TAB = os.path.join(HERE, "..", "federated_lca", "tables")


def theme():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
    })


def score(fed, gt, gt_labels, tag):
    col, _ = fl.hungarian_align(fed, gt, metric="corr")
    per, overall = fl.evaluate_vs_gt(fed, gt, col, None, gt_labels)
    null = fl.permutation_null(fed, gt)
    p = float((null <= overall).mean())
    print(f"\n--- {tag} ---")
    print(per.round(3).to_string(index=False))
    print(f"mean r = {per['pearson_r'].mean():.3f}   overall RMSE = {overall:.4f}   "
          f"perm p = {p:.4f}")
    return col, per, overall, p


def main():
    theme()
    gt_raw = pd.read_csv(os.path.join(TAB, "gt_raw_normalized.csv"), index_col=0)
    ABX = list(gt_raw.columns)
    gt_labels = list(gt_raw.index)
    raw, Z, meta, mean, sd = fl.collect_log2(NODES, ABX, mode="pooled")
    print(f"countries {meta['country'].nunique()}  profiles {Z.shape[0]}  features {Z.shape[1]}")
    print(f"\nFederated POOLED standardization from node sufficient statistics: "
          f"mean={mean[0]:.3f}  sd={sd[0]:.3f}")

    _, Z_abx, _, _, _ = fl.collect_log2(NODES, ABX, mode="per_abx")

    lab, cent = fl.central_consensus_cluster(Z, K=K_GLOBAL, method="kmeans")
    meta["global_class"] = [f"Fed{k+1}" for k in lab]
    fed = pd.DataFrame(cent, index=[f"Fed{k+1}" for k in range(K_GLOBAL)], columns=ABX)
    fed.to_csv(os.path.join(TAB, "federated_global_classes.csv"))
    meta.to_csv(os.path.join(TAB, "country_class_assignment.csv"), index=False)
    print("\nFederated global class profiles:")
    print(fed.round(3).to_string())
    print("\nprofiles per global class:")
    print(meta["global_class"].value_counts().sort_index().to_string())

if __name__ == "__main__":
    main()
