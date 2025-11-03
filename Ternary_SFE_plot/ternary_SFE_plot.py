# ternary_SFE_350K.py
#
# Ternary stacking fault energy diagrams at 350 K for
#   γ_ESF, γ_ISF and γ_Twin
# using SFE_results_AlCrCo.csv (Al–Co–Cr ternary).

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

CSV_FILE = "SFE_results_AlCrCo.csv"
TEMP_K = 350

# --------------------------------------------------
# 1. Load data and keep only 350 K
# --------------------------------------------------
df = pd.read_csv(CSV_FILE)

if TEMP_K not in df["temperature"].unique():
    raise ValueError(f"{TEMP_K} K not found in file.")

df = df[df["temperature"] == TEMP_K].copy()

# quick check for missing SFE values
if df[["gamma_ESF", "gamma_ISF", "gamma_Twin"]].isna().any().any():
    print("WARNING: some SFE values are missing at 350 K.")

# --------------------------------------------------
# 2. Parse composition: CompXX_AlYY_CoZZ_CrWW -> fractions
# --------------------------------------------------
def parse_comp(s):
    # drop 'Comp##'
    parts = s.split("_")[1:]
    comp = {}
    for p in parts:
        i = 0
        while i < len(p) and not p[i].isdigit():
            i += 1
        elem = p[:i]          # Al / Co / Cr
        pct = float(p[i:])    # e.g. 25
        comp[elem] = pct / 100.0

    # ensure all keys
    for k in ("Al", "Co", "Cr"):
        comp.setdefault(k, 0.0)

    # normalise (just in case of rounding)
    ssum = comp["Al"] + comp["Co"] + comp["Cr"]
    return comp["Al"]/ssum, comp["Co"]/ssum, comp["Cr"]/ssum

fracs = df["composition"].apply(parse_comp)
df["Al_frac"] = fracs.apply(lambda t: t[0])
df["Co_frac"] = fracs.apply(lambda t: t[1])
df["Cr_frac"] = fracs.apply(lambda t: t[2])

# --------------------------------------------------
# 3. Map barycentric (Al, Co, Cr) -> Cartesian
#    vertices: Cr bottom-left, Co bottom-right, Al top
# --------------------------------------------------
SQ3_2 = math.sqrt(3) / 2.0
A = np.array([0.5, SQ3_2])  # Al
B = np.array([1.0, 0.0])    # Co
C = np.array([0.0, 0.0])    # Cr

def bary_to_xy(a, b, c):
    return a*A + b*B + c*C

xy = np.vstack([
    bary_to_xy(a, b, c)
    for a, b, c in zip(df["Al_frac"], df["Co_frac"], df["Cr_frac"])
])
df["x"] = xy[:, 0]
df["y"] = xy[:, 1]

# --------------------------------------------------
# 4. Helper: draw triangle & axis labels
# --------------------------------------------------
def draw_triangle(ax):
    tri_x = [C[0], B[0], A[0], C[0]]
    tri_y = [C[1], B[1], A[1], C[1]]
    ax.plot(tri_x, tri_y, "k-", linewidth=1.0)
    ax.set_aspect("equal")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, SQ3_2 + 0.03)
    ax.set_xticks([])
    ax.set_yticks([])

    # axis labels (similar to your figure)
    ax.text(0.5, -0.08, "Cr-fraction", ha="center", va="top", fontsize=10)
    ax.text(0.06, 0.5, "Al-fraction", rotation=60,
            ha="center", va="center", fontsize=10)
    ax.text(0.94, 0.5, "Co-fraction", rotation=-60,
            ha="center", va="center", fontsize=10)

# --------------------------------------------------
# 5. Generic ternary contour plotting function
# --------------------------------------------------
def plot_ternary(ax, values, title, cbar_label):
    # convert J/m² → mJ/m² (dataset values ~1e4)
    z = values / 1e5

    tri = Triangulation(df["x"], df["y"])
    cont = ax.tricontourf(tri, z, levels=14, cmap="plasma")

    draw_triangle(ax)
    ax.set_title(title, fontsize=11)

    cbar = plt.colorbar(cont, ax=ax, pad=0.02, fraction=0.05)
    cbar.set_label(cbar_label, fontsize=10)

# --------------------------------------------------
# 6. Make three diagrams: ESF, ISF, Twin
# --------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)

plot_ternary(
    axes[0],
    df["gamma_ESF"],
    title=f"γ_ESF (mJ/m²) at {TEMP_K} K",
    cbar_label="γ_ESF (mJ/m²)",
)
plot_ternary(
    axes[1],
    df["gamma_ISF"],
    title=f"γ_ISF (mJ/m²) at {TEMP_K} K",
    cbar_label="γ_ISF (mJ/m²)",
)
plot_ternary(
    axes[2],
    df["gamma_Twin"],
    title=f"γ_Twin (mJ/m²) at {TEMP_K} K",
    cbar_label="γ_Twin (mJ/m²)",
)

plt.suptitle("Ternary stacking fault energies at 350 K (Al–Co–Cr)", y=1.03)
plt.savefig("ternary_SFE_350K.png", dpi=300, bbox_inches="tight")
plt.show()
