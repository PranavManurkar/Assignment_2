# all_sfe_plots.py
#
# Reads SFE_results_AlCrCo.csv and produces one figure per
# energy type (ISF, ESF, Twin, ΔE_dhcp-fcc, ΔE_hcp-fcc).

import pandas as pd
import matplotlib.pyplot as plt

# -------- 1. Read data --------
filename = "SFE_results_AlCrCo.csv"   # change path if needed
df = pd.read_csv(filename)

# -------- 2. Build composition labels (same style as example figure) --------
def make_comp_label(comp):
    """
    Convert 'Comp01_Al100_Co00_Cr00'
    --> 'Al1.00_Co0.00_Cr0.00'
    """
    parts = comp.split("_")[1:]  # ignore 'Comp##'
    new_parts = []
    for p in parts:
        i = 0
        while i < len(p) and not p[i].isdigit():
            i += 1
        elem = p[:i]              # element symbol, e.g. 'Al'
        perc = float(p[i:])       # e.g. 100
        frac = perc / 100.0       # 1.00
        new_parts.append(f"{elem}{frac:.2f}")
    return "_".join(new_parts)

df["comp_label"] = df["composition"].apply(make_comp_label)

# -------- 3. Define which energies to plot & scale factors --------
# NOTE: scale_factor is chosen so the y-values are of order 0.1,
# similar to the example plot. Adjust if your units are different.
ENERGIES = [
    ("gamma_ISF",        "Intrinsic Stacking Fault Energy (ISF)",   1e5),
    ("gamma_ESF",        "Extrinsic Stacking Fault Energy (ESF)",   1e5),
    ("gamma_Twin",       "Twin Fault Energy",                       1e5),
    ("delta_E_dhcp_fcc", "ΔE_dhcp-fcc",                             1e3),
    ("delta_E_hcp_fcc",  "ΔE_hcp-fcc",                              1e3),
]

# -------- 4. Loop over each energy type and create plots --------
temps = sorted(df["temperature"].unique())

for col, nice_name, scale_factor in ENERGIES:
    df[col + "_scaled"] = df[col] / scale_factor

    fig, ax = plt.subplots(figsize=(10, 7))

    for comp, group in df.groupby("comp_label"):
        group = group.sort_values("temperature")
        ax.plot(
            group["temperature"],
            group[col + "_scaled"],
            marker="o",
            linewidth=1.5,
            markersize=5,
            label=comp,
        )

    ax.set_title(f"Variation of {nice_name} with Temperature")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(f"{nice_name} (mJ/m²)")  # label like the example
    ax.set_xticks(temps)

    # grid & legend style similar to your attached figure
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(
        title="Composition",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        fontsize="small",
    )

    plt.tight_layout()

    outname = f"{col}_vs_Temperature.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close(fig)

print("Finished. Saved one PNG per energy type in the current folder.")
