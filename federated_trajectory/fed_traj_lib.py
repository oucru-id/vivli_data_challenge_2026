#!/usr/bin/env python3

from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr

_LOG2 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", "federated_lca_log2"))
sys.path.insert(0, _LOG2)
from fed_log2_lib import aligned_aggregate  

YEARS = np.arange(2004, 2025)
CLASSES = [f"Class{i}" for i in range(1, 6)]


# LOCAL NODE SIDE                                                                    
def local_trajectory(node_dir):
    return pd.read_csv(os.path.join(node_dir, "trajectory.csv"))

def local_weight(node_dir):
    p = os.path.join(node_dir, "weight.txt")
    if os.path.exists(p):
        with open(p) as fh:
            return float(fh.read().strip())
    return 1.0


# GLOBAL AGGREGATION SIDE                                                  
def collect_trajectories(nodes_root, min_cell=0):
    frames, weights = [], {}
    for country in sorted(os.listdir(nodes_root)):
        node = os.path.join(nodes_root, country)
        if not os.path.isdir(node):
            continue
        frames.append(local_trajectory(node))
        weights[country] = local_weight(node)
    df = pd.concat(frames, ignore_index=True)

    supp = {"cells_total": len(df), "cells_suppressed": 0, "N_suppressed": 0}
    if min_cell and min_cell > 0:
        m = (df["N"] > 0) & (df["N"] < min_cell)
        supp["cells_suppressed"] = int(m.sum())
        supp["N_suppressed"] = int(df.loc[m, "N"].sum())
        df.loc[m, "N"] = 0
    return df, weights, supp

def align_from_mic(log2_tables=None):
    log2_tables = log2_tables or os.path.join(_LOG2, "tables")
    a = pd.read_csv(os.path.join(log2_tables, "country_class_assignment_aligned.csv"))
    m = pd.read_csv(os.path.join(log2_tables, "gt_match_aligned.csv"))
    f2g = dict(zip(m["fed_class"], m["gt_class"]))
    a = a.rename(columns={"country": "Country", "local_class": "Class"})
    a["global_class"] = a["global_class"].map(f2g)
    return a[["Country", "Class", "global_class"]]


def curve_matrix(df):
    tot = df.groupby(["Country", "Year"])["N"].transform("sum")
    d = df.assign(share=np.where(tot > 0, df["N"] / tot, 0.0))
    wide = d.pivot_table(index=["Country", "Class"], columns="Year",
                         values="share", aggfunc="sum").reindex(columns=YEARS)
    wide = wide.fillna(0.0)
    meta = wide.index.to_frame(index=False)
    meta.columns = ["country", "local_class"]
    return wide.to_numpy(dtype=float), meta


def align_from_curves(df, weights, K=5, n_seeds=20):
    X, meta = curve_matrix(df)
    w = np.array([weights.get(c, 1.0) for c in meta["country"]], dtype=float)
    _, assign, obj = aligned_aggregate(X, w, meta, K=K, n_seeds=n_seeds)
    out = meta.copy()
    out["global_class"] = [f"G{k+1}" for k in assign]
    out = out.rename(columns={"country": "Country", "local_class": "Class"})
    return out[["Country", "Class", "global_class"]], obj

def aggregate_global(df, mapping):
    j = df.merge(mapping, on=["Country", "Class"], how="left")
    if j["global_class"].isna().any():
        raise ValueError(f"{int(j['global_class'].isna().sum())} unmapped rows")
    g = j.groupby(["Year", "global_class"])["N"].sum().reset_index()
    tot = g.groupby("Year")["N"].transform("sum")
    g["pct"] = np.where(tot > 0, 100 * g["N"] / tot, 0.0)
    return g


def to_matrix(long_df, year_col, class_col, val_col, classes=None):
    piv = long_df.pivot_table(index=year_col, columns=class_col, values=val_col,
                              aggfunc="sum").reindex(index=YEARS)
    if classes is not None:
        piv = piv.reindex(columns=classes)
    return piv.fillna(0.0)


def hungarian_match_trajectories(fed_mat, gt_mat):
    F, G = fed_mat.to_numpy().T, gt_mat.to_numpy().T  
    C = np.zeros((F.shape[0], G.shape[0]))
    for i in range(F.shape[0]):
        for j in range(G.shape[0]):
            if np.std(F[i]) < 1e-12 or np.std(G[j]) < 1e-12:
                C[i, j] = np.mean(np.abs(F[i] - G[j]))
            else:
                r = pearsonr(F[i], G[j])[0]
                C[i, j] = -(r if np.isfinite(r) else 0.0)
    row, col = linear_sum_assignment(C)
    return {fed_mat.columns[i]: gt_mat.columns[j] for i, j in zip(row, col)}


def evaluate_vs_gt(fed_mat, gt_mat, mapping):
    rows = []
    for f_cls, g_cls in mapping.items():
        f, g = fed_mat[f_cls].to_numpy(), gt_mat[g_cls].to_numpy()
        r = pearsonr(f, g)[0] if np.std(f) > 1e-12 and np.std(g) > 1e-12 else np.nan
        rows.append({"fed_class": f_cls, "gt_class": g_cls,
                     "mae_pp": float(np.mean(np.abs(f - g))),
                     "rmse_pp": float(np.sqrt(np.mean((f - g) ** 2))),
                     "pearson_r": float(r)})
    per = pd.DataFrame(rows).sort_values("gt_class").reset_index(drop=True)
    F = np.column_stack([fed_mat[k].to_numpy() for k in mapping])
    G = np.column_stack([gt_mat[v].to_numpy() for v in mapping.values()])
    return per, float(np.mean(np.abs(F - G)))


def permutation_null(fed_mat, gt_mat, n=2000, seed=7):
    rng = np.random.default_rng(seed)
    F = fed_mat.to_numpy()
    G = gt_mat.to_numpy()
    k = G.shape[1]
    out = []
    for _ in range(n):
        out.append(np.mean(np.abs(F - G[:, rng.permutation(k)])))
    return np.array(out)

def logistic_trend(years, n_class, n_total):
    from scipy.optimize import minimize
    yc = np.asarray(years, float) - 2004.0
    k = np.asarray(n_class, float)
    N = np.asarray(n_total, float)
    ok = N > 0
    yc, k, N = yc[ok], k[ok], N[ok]

    def nll(p):
        z = np.clip(p[0] + p[1] * yc, -30, 30)
        pr = np.clip(1 / (1 + np.exp(-z)), 1e-9, 1 - 1e-9)
        return -np.sum(k * np.log(pr) + (N - k) * np.log(1 - pr))

    r = minimize(nll, [0.0, 0.0], method="Nelder-Mead", options={"maxiter": 5000})
    return float(r.x[0]), float(r.x[1])


def boot_trend_ci(years, n_class, n_total, B=400, seed=11):
    rng = np.random.default_rng(seed)
    N = np.asarray(n_total, int)
    p = np.clip(np.asarray(n_class, float) / np.maximum(N, 1), 0, 1)
    betas = []
    for _ in range(B):
        ks = rng.binomial(N, p)
        try:
            betas.append(logistic_trend(years, ks, N)[1])
        except Exception:
            pass
    return float(np.percentile(betas, 2.5)), float(np.percentile(betas, 97.5))


def emergence_year(years, pct, thresholds=(1.0, 5.0, 10.0)):
    years = np.asarray(years); pct = np.asarray(pct, float)
    return {t: (int(years[pct >= t][0]) if (pct >= t).any() else None) for t in thresholds}


def changepoint(years, pct):
    years = np.asarray(years); y = np.asarray(pct, float)
    best, best_sse = None, np.inf
    for i in range(2, len(y) - 1):
        sse = ((y[:i] - y[:i].mean()) ** 2).sum() + ((y[i:] - y[i:].mean()) ** 2).sum()
        if sse < best_sse:
            best_sse, best = sse, int(years[i])
    return best
