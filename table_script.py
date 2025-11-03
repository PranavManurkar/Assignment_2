# plot_sfe_heatmaps.py
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ------------ edit if needed ------------
CSV_PATH = "sfe_data.csv"
OUT_FIG = "SFE_heatmaps.png"
# ----------------------------------------

# Read data
df = pd.read_csv(CSV_PATH)

# Ensure numeric columns
df["Al_at_pct"] = pd.to_numeric(df["Al_at_pct"], errors="coerce")
df["Co_at_pct"] = pd.to_numeric(df["Co_at_pct"], errors="coerce")
df["Cr_at_pct"] = pd.to_numeric(df["Cr_at_pct"], errors="coerce")
df["Temperature (K)"] = pd.to_numeric(df["Temperature (K)"], errors="coerce")

# Create a friendly composition label and ordering key
df["Composition"] = df["Alloy"]  # already set in CSV as Al-xx-Co-yy-Cr-zz
# order compositions by Al_at_pct for plotting left->right
comp_order = (
    df[["Composition", "Al_at_pct"]]
    .drop_duplicates()
    .sort_values("Al_at_pct")["Composition"]
    .tolist()
)

# Structures to plot (keep the order consistent)
structures = sorted(df["Structure"].unique())

# Temperatures (sorted)
temp_order = sorted(df["Temperature (K)"].dropna().unique())

# Prepare the figure: 1 row, N columns (one per structure)
n_struct = len(structures)
fig, axes = plt.subplots(1, n_struct, figsize=(4.5 * n_struct + 2, 6), constrained_layout=True)

if n_struct == 1:
    axes = [axes]

# Plot each structure as a heatmap
for ax, struct in zip(axes, structures):
    sub = df[df["Structure"] == struct].copy()
    # pivot: rows=temp, cols=composition
    pivot = sub.pivot(index="Temperature (K)", columns="Composition", values="SFE (mJ/m^2)")
    # ensure temperature order
    pivot = pivot.reindex(index=temp_order)
    # ensure column order by composition Al%
    cols_present = [c for c in comp_order if c in pivot.columns]
    pivot = pivot[cols_present]

    # mask for NaNs so they draw as blank
    mask = pivot.isnull()

    sns.heatmap(
        pivot,
        ax=ax,
        mask=mask,
        annot=False,         # set to True if you filled SFE and want annotations
        fmt=".2f",
        cmap="viridis",
        cbar=True if ax is axes[-1] else False,
        linewidths=0.3,
        linecolor="gray"
    )

    ax.set_title(f"{struct.upper()} structure")
    ax.set_xlabel("Composition (Al-Co-Cr)")
    ax.set_ylabel("Temperature (K)")
    # improve x-label readability
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)

# Save figure
plt.suptitle("Stacking Fault Energy (SFE) — Al-Co-Cr (empty cells = missing SFE)", y=1.02, fontsize=14)
plt.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
print(f"Saved heatmap figure to: {OUT_FIG}")
plt.show()
