#!/usr/bin/env python3

from __future__ import annotations
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASEDIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(BASEDIR, "mic_heatmap_data_all_countries.csv")
GT = os.path.join(BASEDIR, "mic_heatmap_global_GT.csv")
ATLAS = os.path.join(BASEDIR, "atlas_vivli_2004_2024.csv")
NODES = os.path.join(HERE, "nodes")
TAB = os.path.join(HERE, "tables")
os.makedirs(NODES, exist_ok=True)
os.makedirs(TAB, exist_ok=True)


def isolate_counts():
    a = pd.read_csv(ATLAS, usecols=["Species", "Source", "Country"], dtype=str)
    for c in a.columns:
        a[c] = a[c].astype(str).str.strip()
    kb = a[(a["Species"] == "Klebsiella pneumoniae") & (a["Source"] == "Blood")]
    return kb.groupby("Country").size()


def norm_abx(name: str) -> str:
    return str(name).replace(".", " ").strip()


def main():
    gt = pd.read_csv(GT, index_col=0)
    gt.columns = [norm_abx(c) for c in gt.columns]
    gt_abx = list(gt.columns)

    d = pd.read_csv(SRC)
    d["Antibiotic"] = d["Antibiotic"].map(norm_abx)
    d = d[d["Antibiotic"].isin(gt_abx)].copy()

    counts = isolate_counts()
    med = float(counts.median())

    reg = []
    for country, sub in d.groupby("Country"):
        node = os.path.join(NODES, country.replace("/", "_"))
        os.makedirs(node, exist_ok=True)
        sub[["Country", "Class", "Antibiotic", "log2_MIC"]].to_csv(
            os.path.join(node, "profiles.csv"), index=False)
        n_iso = int(counts.get(country, med))
        with open(os.path.join(node, "weight.txt"), "w") as fh:
            fh.write(str(n_iso))
        reg.append({"country": country, "K": sub["Class"].nunique(),
                    "n_antibiotics": sub["Antibiotic"].nunique(), "n_isolates": n_iso})

    reg = pd.DataFrame(reg).sort_values("country")
    reg.to_csv(os.path.join(TAB, "node_registry.csv"), index=False)
    gt.to_csv(os.path.join(TAB, "gt_raw_normalized.csv"))

    print(f"nodes written  : {len(reg)}")
    print(f"profile vectors: {int(reg['K'].sum())}")
    print(f"K distribution : {reg['K'].value_counts().sort_index().to_dict()}")
    print(f"payload        : log2_MIC (absolute, cross-country comparable)")
    print(f"Nodes -> {NODES}")


if __name__ == "__main__":
    main()
