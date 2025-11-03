# phase_stability_AlCoCr_350K.py
#
# Phase Stability Map for Al-Co-Cr:
#   ΔE = E_HCP - E_FCC at 350 K
# using SFE_results_AlCrCo.csv

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

CSV_FILE = "SFE_results_AlCrCo.csv"
TEMP_K = 350   # temperature for the map

# ---------------- 1. Load data ----------------
df = pd.read_csv(CSV_FILE)

if TEMP_K not in df["temperature"].unique():
    raise ValueError(f"{TEMP_K} K not found in {CSV_FILE}.")

dfT = df[df["temperature"] == TEMP_K].copy()

if dfT["delta_E_hcp_fcc"].isna().any():
    print("WARNING: some delta_E_hcp_fcc values are missing at", TEMP_K, "K")

# ---------------- 2. Parse compositions ----------------
# Example label: "Comp04_Al25_Co75_Cr00" -> (Al, Co, Cr) fractions
def parse_comp(s):
    parts = s.split("_")[1:]  # drop 'Comp##'
    comp = {}
    for p in parts:
        i = 0
        while i < len(p) and not p[i].isdigit():
            i += 1
        elem = p[:i]          # 'Al', 'Co', 'Cr'
        pct = float(p[i:])    # e.g. 25
        comp[elem] = pct / 100.0

    for k in ("Al", "Co", "Cr"):
        comp.setdefault(k, 0.0)

    ssum = comp["Al"] + comp["Co"] + comp["Cr"]
    return comp["Al"]/ssum, comp["Co"]/ssum, comp["Cr"]/ssum

fracs = dfT["composition"].apply(parse_comp)
dfT["Al_frac"] = fracs.apply(lambda t: t[0])
dfT["Co_frac"] = fracs.apply(lambda t: t[1])
dfT["Cr_frac"] = fracs.apply(lambda t: t[2])

# ---------------- 3. Barycentric -> Cartesian ----------------
# Triangle vertices: Cr (bottom-left), Co (bottom-right), Al (top)
SQ3_2 = math.sqrt(3) / 2.0
A = np.array([0.5, SQ3_2])  # Al
B = np.array([1.0, 0.0])    # Co
C = np.array([0.0, 0.0])    # Cr

def bary_to_xy(a, b, c):
    return a*A + b*B + c*C

xy = np.vstack([
    bary_to_xy(a, b, c)
    for a, b, c in zip(dfT["Al_frac"], dfT["Co_frac"], dfT["Cr_frac"])
])
dfT["x"] = xy[:, 0]
dfT["y"] = xy[:, 1]

# ---------------- 4. Draw main triangle ----------------
def draw_triangle(ax):
    tri_x = [C[0], B[0], A[0], C[0]]
    tri_y = [C[1], B[1], A[1], C[1]]
    ax.plot(tri_x, tri_y, "k-", linewidth=2.0)

    ax.set_aspect("equal")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, SQ3_2 + 0.03)
    ax.set_xticks([])
    ax.set_yticks([])

    # edge axis labels
    ax.text(0.5, -0.08, "Cr fraction", ha="center", va="top", fontsize=10)
    ax.text(0.06, 0.5, "Al fraction", rotation=60,
            ha="center", va="center", fontsize=10)
    ax.text(0.94, 0.5, "Co fraction", rotation=-60,
            ha="center", va="center", fontsize=10)

# ---------------- 5. Make contour plot ----------------
# Convert assumed meV/atom -> eV/atom
deltaE_eV = dfT["delta_E_hcp_fcc"] / 1000.0

fig, ax = plt.subplots(figsize=(6.0, 5.3))

tri = Triangulation(dfT["x"], dfT["y"])
contour = ax.tricontourf(
    tri,
    deltaE_eV,
    levels=14,
    cmap="coolwarm"    # similar to reference figure
)

draw_triangle(ax)

title = f"Phase Stability Map (ΔE = E_HCP - E_FCC at {TEMP_K} K)"
ax.set_title(title)

cbar = plt.colorbar(contour, ax=ax, pad=0.02)
cbar.set_label("ΔE (eV/atom)")

plt.tight_layout()
plt.savefig("phase_stability_AlCoCr_350K.png", dpi=300, bbox_inches="tight")
plt.show()
