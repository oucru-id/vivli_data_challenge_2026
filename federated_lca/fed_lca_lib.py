#!/usr/bin/env python3

from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
_SIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "federated_lca")
sys.path.insert(0, os.path.abspath(_SIB))
from fed_lca_lib import (
    central_consensus_cluster,
    evaluate_vs_gt,
    gt_zscore,
    hungarian_align,
    permutation_null,
)

# LOCAL NODE SIDE                                                                    
def local_profiles_log2(node_dir, antibiotics):
    df = pd.read_csv(os.path.join(node_dir, "profiles.csv"))
    wide = df.pivot_table(index="Class", columns="Antibiotic", values="log2_MIC")
    wide = wide.reindex(columns=antibiotics)
    return list(wide.index), wide.to_numpy(dtype=float)


def local_suffstats(node_dir, antibiotics):
    _, M = local_profiles_log2(node_dir, antibiotics)
    n = np.sum(~np.isnan(M), axis=0).astype(float)
    s = np.nansum(M, axis=0)
    ss = np.nansum(M ** 2, axis=0)
    return n, s, ss


def local_weight(node_dir):
    p = os.path.join(node_dir, "weight.txt")
    if os.path.exists(p):
        with open(p) as fh:
            return float(fh.read().strip())
    return 1.0

# GLOBAL AGGREGATION SIDE                                                                  
def central_global_standardize(suffstats, mode="pooled"):

    N = np.sum([n for n, _, _ in suffstats], axis=0)
    S = np.sum([s for _, s, _ in suffstats], axis=0)
    SS = np.sum([ss for _, _, ss in suffstats], axis=0)
    if mode == "per_abx":
        mean = S / N
        var = np.clip(SS / N - mean ** 2, 1e-12, None)
        return mean, np.sqrt(var)
    if mode == "pooled":
        n, s, ss = N.sum(), S.sum(), SS.sum()
        mean = s / n
        var = max(ss / n - mean ** 2, 1e-12)
        k = len(N)
        return np.full(k, mean), np.full(k, np.sqrt(var))
    raise ValueError(mode)


def collect_log2(nodes_root, antibiotics, mode="pooled"):
    suff, rows, meta = [], [], []
    for country in sorted(os.listdir(nodes_root)):
        node = os.path.join(nodes_root, country)
        if not os.path.isdir(node):
            continue
        suff.append(local_suffstats(node, antibiotics))
        labels, M = local_profiles_log2(node, antibiotics)
        for lab, vec in zip(labels, M):
            rows.append(vec)
            meta.append({"country": country, "local_class": lab})

    raw = np.vstack(rows)
    mean, sd = central_global_standardize(suff, mode=mode)
    Z = (raw - mean) / sd
    return raw, Z, pd.DataFrame(meta), mean, sd


def collect_weights(nodes_root, meta):
    wc = {}
    for country in sorted(os.listdir(nodes_root)):
        node = os.path.join(nodes_root, country)
        if os.path.isdir(node):
            wc[country] = local_weight(node)
    med = float(np.median(list(wc.values())))
    return np.array([wc.get(c, med) for c in meta["country"]], dtype=float)

def aligned_aggregate(Z, w, meta, K=5, n_seeds=20, iters=200):
    groups = [np.where((meta["country"] == c).values)[0]
              for c in meta["country"].unique()]

    def one_run(C):
        for _ in range(iters):
            assign = np.full(len(Z), -1)
            for idx in groups:
                cost = ((Z[idx][:, None, :] - C[None, :, :]) ** 2).sum(2)
                r, cc = linear_sum_assignment(cost)   
                assign[idx[r]] = cc
            newC = np.vstack([np.average(Z[assign == k], axis=0, weights=w[assign == k])
                              if (assign == k).any() else C[k] for k in range(K)])
            if np.allclose(newC, C, atol=1e-9):
                C = newC
                break
            C = newC
        obj = 0.0
        assign = np.full(len(Z), -1)
        for idx in groups:
            cost = ((Z[idx][:, None, :] - C[None, :, :]) ** 2).sum(2)
            r, cc = linear_sum_assignment(cost)
            assign[idx[r]] = cc
            obj += (w[idx[r]] * cost[r, cc]).sum()
        return C, assign, float(obj)

    best = None
    for seed in range(n_seeds):
        C0 = KMeans(K, n_init=25, random_state=seed).fit(Z, sample_weight=w).cluster_centers_
        C, assign, obj = one_run(C0)
        if best is None or obj < best[2]:
            best = (C, assign, obj)
    return best
