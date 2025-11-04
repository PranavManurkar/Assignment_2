#!/usr/bin/env python3
"""
compare_DMLF_Ising_ESF.py

Compare DMLF extrinsic stacking fault energy (gamma_ESF)
with a simple axial Ising (AIM) estimate for the Al–Co–Cr system.

Inputs (same folder):
  - SFE_results_AlCrCo.csv

Outputs:
  - SFE_DMLF_Ising_ESF_comparison.csv
  - DMLF_vs_Ising_ESF_100K.png / _350K.png / _550K.png
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

E_CHARGE = 1.602176634e-19  # J per eV

# -------------------- I/O and data prep --------------------

def load_sfe_data(csv_path: str = "SFE_results_AlCrCo.csv") -> pd.DataFrame:
    """
    Normalize column names to internal schema.

    File columns (from your data):
      composition, temperature,
      E_fcc, E_hcp, E_dhcp,
      A_fcc (Å^2), natoms,
      gamma_ISF, gamma_ESF, gamma_Twin,
      delta_E_dhcp_fcc, delta_E_hcp_fcc
    """
    df = pd.read_csv(csv_path)

    # Make composition name consistent
    if "CompName" not in df.columns and "composition" in df.columns:
        df.rename(columns={"composition": "CompName"}, inplace=True)

    rename_map = {
        "temperature": "Temp_K",
        "E_fcc": "E_FCC_eV",
        "E_hcp": "E_HCP_eV",
        "E_dhcp": "E_DHCP_eV",
        "A_fcc": "A_FCC_A2",  # Å^2
        "gamma_ISF": "gamma_ISF_DMLF_mJm2",
        "gamma_ESF": "gamma_ESF_DMLF_mJm2",
        "gamma_Twin": "gamma_Twin_DMLF_mJm2",
        "delta_E_hcp_fcc": "delta_E_hcp_fcc_eV",
        "delta_E_dhcp_fcc": "delta_E_dhcp_fcc_eV",
    }
    df.rename(columns=rename_map, inplace=True)

    required = [
        "CompName", "Temp_K", "E_FCC_eV", "E_HCP_eV",
        "A_FCC_A2", "natoms", "gamma_ISF_DMLF_mJm2", "gamma_ESF_DMLF_mJm2"
    ]
    for c in required:
        if c not in df.columns:
            raise ValueError(
                f"Required column '{c}' not found in {csv_path}. "
                f"Available: {list(df.columns)}"
            )

    df["Temp_K"] = df["Temp_K"].astype(int)
    df["natoms"] = pd.to_numeric(df["natoms"], errors="coerce").astype(int)

    num_cols = [
        "E_FCC_eV", "E_HCP_eV", "E_DHCP_eV",
        "A_FCC_A2",
        "gamma_ISF_DMLF_mJm2", "gamma_ESF_DMLF_mJm2", "gamma_Twin_DMLF_mJm2",
        "delta_E_hcp_fcc_eV", "delta_E_dhcp_fcc_eV"
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # If delta_E_hcp_fcc not present, compute it
    if "delta_E_hcp_fcc_eV" not in df.columns or df["delta_E_hcp_fcc_eV"].isna().all():
        df["delta_E_hcp_fcc_eV"] = df["E_HCP_eV"] - df["E_FCC_eV"]

    return df


def compute_ising_gammas(df: pd.DataFrame) -> pd.DataFrame:
    """
    First compute an axial Ising estimate of intrinsic SFE:

        gamma_ISF^Ising = 2 * (E_hcp - E_fcc) * natoms / A_fcc

    (E in eV/atom, natoms dimensionless, A_fcc in Å^2, with conversions.)

    Then assume extrinsic SF is roughly two intrinsic SFs:

        gamma_ESF^Ising = 2 * gamma_ISF^Ising

    Adds:
      - gamma_ISF_Ising_mJm2
      - gamma_ESF_Ising_mJm2
    """
    df = df.copy()

    # convert Å^2 -> m^2
    A_m2 = df["A_FCC_A2"] * 1e-20
    dE_cell_eV = (df["E_HCP_eV"] - df["E_FCC_eV"]) * df["natoms"]

    gamma_ISF_Jm2 = (2.0 * dE_cell_eV * E_CHARGE) / A_m2
    df["gamma_ISF_Ising_mJm2"] = gamma_ISF_Jm2 * 1e3

    df["gamma_ESF_Ising_mJm2"] = 2.0 * df["gamma_ISF_Ising_mJm2"]

    return df


# -------------------- colors / labels --------------------

_comp_pat = re.compile(r"Al(\d+)_Co(\d+)_Cr(\d+)", re.IGNORECASE)

def pretty_label(compname: str) -> str:
    m = _comp_pat.search(compname)
    if not m:
        return compname
    return f"Al{m.group(1)} Co{m.group(2)} Cr{m.group(3)}"

def build_color_map(df_all: pd.DataFrame):
    comps = sorted(df_all["CompName"].unique())
    n = len(comps)
    cmap = plt.cm.get_cmap("tab20", n if n > 1 else 2)
    colors = [cmap(i) for i in range(n)]
    return {comp: colors[i] for i, comp in enumerate(comps)}


# -------------------- plotting --------------------

def plot_dmlf_vs_ising(df: pd.DataFrame,
                       color_map: dict,
                       out_prefix: str = "DMLF_vs_Ising_ESF"):
    temps = sorted(df["Temp_K"].dropna().unique())

    for T in temps:
        sub = df[(df["Temp_K"] == T) &
                 df["gamma_ESF_DMLF_mJm2"].notna() &
                 df["gamma_ESF_Ising_mJm2"].notna()].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(8.5, 5.6))

        handles, labels = [], []
        for comp, grp in sub.groupby("CompName"):
            x = grp["gamma_ESF_DMLF_mJm2"].to_numpy()
            y = grp["gamma_ESF_Ising_mJm2"].to_numpy()
            color = color_map.get(comp, "C0")
            sc = ax.scatter(x, y, s=36, color=color, edgecolors="k",
                            linewidths=0.4, alpha=0.95)
            handles.append(sc)
            labels.append(pretty_label(comp))

        allx = sub["gamma_ESF_DMLF_mJm2"].to_numpy()
        ally = sub["gamma_ESF_Ising_mJm2"].to_numpy()
        lo = float(min(allx.min(), ally.min()))
        hi = float(max(allx.max(), ally.max()))
        pad = 0.06 * (hi - lo) if hi > lo else 1.0
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                "k--", lw=1.0, label="y = x")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)

        ax.set_xlabel(r"DMLF $\gamma_{\mathrm{ESF}}$ [mJ/m$^2$]")
        ax.set_ylabel(r"Ising (AIM) $\gamma_{\mathrm{ESF}}$ [mJ/m$^2$]")
        ax.set_title(f"Extrinsic SFE: DMLF vs Ising (AIM) at {T} K")
        ax.grid(True, linestyle="--", alpha=0.35)

        lgd = ax.legend(handles, labels, title="Composition (at.%)",
                        loc="center left", bbox_to_anchor=(1.02, 0.5),
                        borderaxespad=0.0, frameon=True, fontsize="small",
                        title_fontsize="small", ncol=1)

        fig.tight_layout()
        fname = f"{out_prefix}_{T}K.png"
        fig.savefig(fname, dpi=300, bbox_extra_artists=(lgd,),
                    bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {fname}")


# -------------------- main --------------------

def main():
    sfe_file = "SFE_results_AlCrCo.csv"
    if not Path(sfe_file).exists():
        raise FileNotFoundError(
            f"{sfe_file} not found. Run from the Assignment_2 folder."
        )

    df = load_sfe_data(sfe_file)
    df = compute_ising_gammas(df)

    color_map = build_color_map(df)

    out_csv = "SFE_DMLF_Ising_ESF_comparison.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    plot_dmlf_vs_ising(df, color_map)


if __name__ == "__main__":
    main()
