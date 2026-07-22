#!/usr/bin/env python3

from __future__ import annotations
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASEDIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(BASEDIR, "all_countries_trajectory_data.csv")
GT = os.path.join(BASEDIR, "GT_global_trajectory_LCA.csv")
NODES = os.path.join(HERE, "nodes")
TAB = os.path.join(HERE, "tables")
os.makedirs(NODES, exist_ok=True)
os.makedirs(TAB, exist_ok=True)


def main():
    t = pd.read_csv(SRC)
    if "Uncertainty" in t.columns:
        n_unc = t["Uncertainty"].notna().sum()
        print(f"dropping 'Uncertainty' column (non-null values: {n_unc})")
        t = t.drop(columns=["Uncertainty"])

    gt = pd.read_csv(GT)
    gt.to_csv(os.path.join(TAB, "gt_trajectory.csv"), index=False)

    reg = []
    for country, sub in t.groupby("Country"):
        node = os.path.join(NODES, str(country).replace("/", "_"))
        os.makedirs(node, exist_ok=True)
        sub[["Country", "Year", "Class", "N", "Total_N"]].to_csv(
            os.path.join(node, "trajectory.csv"), index=False)
        total = int(sub["N"].sum())
        with open(os.path.join(node, "weight.txt"), "w") as fh:
            fh.write(str(total))
        reg.append({"country": country, "K": sub["Class"].nunique(),
                    "year_min": int(sub["Year"].min()), "year_max": int(sub["Year"].max()),
                    "n_years": sub["Year"].nunique(), "total_N": total})

    reg = pd.DataFrame(reg).sort_values("country")
    reg.to_csv(os.path.join(TAB, "node_registry.csv"), index=False)

if __name__ == "__main__":
    main()
