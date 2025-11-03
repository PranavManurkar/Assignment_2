# lattice_all.py
#
# 1) Collect lattice parameters for EVERY composition and EVERY temperature
#    from all CompXX_Al*_Co*_Cr* folders and write:
#       lattice_all_compositions_all_temps.csv
#
# 2) For each Structure (FCC, HCP, DHCP) and Temp (100, 350, 550 K),
#    generate a ternary contour map of a_final_A over Al–Co–Cr composition.

from pathlib import Path
import math
import re
import csv

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# Temperatures and structures used in the simulations
TEMPERATURES = [100, 350, 550]
STRUCTURES = ["FCC", "HCP", "DHCP"]

# Reference pure FCC lattice parameters (Å)
# HCP/DHCP initial a is taken as a_fcc / sqrt(2)
PURE_FCC_A = {"Al": 4.05, "Co": 3.544, "Cr": 2.885}

# ---------- helpers for reading data ----------

# Matches Al/Co/Cr percentages in folder names like Comp01_Al100_Co00_Cr00
name_pat = re.compile(r"Al(?P<Al>\d+)_Co(?P<Co>\d+)_Cr(?P<Cr>\d+)")


def parse_comp_dir_name(name: str):
    """Return (Al%, Co%, Cr%) from folder name."""
    m = name_pat.search(name)
    if not m:
        raise ValueError(f"Folder name does not contain Al/Co/Cr %: {name}")
    al = int(m.group("Al"))
    co = int(m.group("Co"))
    cr = int(m.group("Cr"))
    if al + co + cr != 100:
        print(f"WARNING: {name}: Al+Co+Cr = {al+co+cr} (expected 100)")
    return al, co, cr


def initial_lattice_params(al_pct, co_pct, cr_pct):
    """Vegard's law for FCC, then derive HCP/DHCP a from FCC."""
    alf, cof, crf = al_pct / 100.0, co_pct / 100.0, cr_pct / 100.0
    a_fcc0 = (
        alf * PURE_FCC_A["Al"]
        + cof * PURE_FCC_A["Co"]
        + crf * PURE_FCC_A["Cr"]
    )
    a_hcp0 = a_fcc0 / math.sqrt(2.0)
    a_dhcp0 = a_hcp0
    return a_fcc0, a_hcp0, a_dhcp0


def read_results_summary(comp_dir: Path):
    """
    Reads comp_dir/results_summary.txt.
    Each valid line is:
      STRUCT TEMP natoms pe_per_atom volume lx ly lz area
    Returns dict[(STRUCT, TEMP)] -> record.
    """
    p = comp_dir / "results_summary.txt"
    data = {}
    if not p.exists():
        print(f"WARNING: {p} not found")
        return data

    with p.open() as f:
        for ln in f:
            parts = ln.split()
            if len(parts) < 9:
                continue
            struct = parts[0].upper()
            if struct not in STRUCTURES:
                continue
            try:
                temp = int(parts[1])
                natoms = int(parts[2])
                volume = float(parts[4])
                lx = float(parts[5])
                ly = float(parts[6])
                lz = float(parts[7])
            except Exception:
                continue
            data[(struct, temp)] = {
                "natoms": natoms,
                "volume": volume,
                "lx": lx,
                "ly": ly,
                "lz": lz,
            }
    return data


def infer_nx_from_lx(lx, a0):
    """Estimate integer replication along x from lx and initial a0."""
    if a0 <= 0:
        return None
    est = max(1, round(lx / a0))
    return est


def final_a(struct, rec, a0_guess):
    """
    Compute final lattice parameter a (Å) from relaxed cell.
    FCC:  a from volume per atom
    HCP/DHCP:  lx ≈ nx * a
    """
    if rec is None:
        return None
    if struct == "FCC":
        # V = (n/4) * a^3  → a = (4V/n)^(1/3)
        return (4.0 * rec["volume"] / rec["natoms"]) ** (1.0 / 3.0)
    else:
        nx = infer_nx_from_lx(rec["lx"], a0_guess)
        return rec["lx"] / nx if nx else None


# ---------- collect rows & write CSV ----------

def collect_rows():
    root = Path(".")
    comp_dirs = sorted(
        d for d in root.iterdir() if d.is_dir() and d.name.startswith("Comp")
    )
    if not comp_dirs:
        raise RuntimeError("No Comp* directories found. Run from Assignment_2 folder.")

    rows = []

    for comp_dir in comp_dirs:
        name = comp_dir.name
        al_pct, co_pct, cr_pct = parse_comp_dir_name(name)
        a0_fcc, a0_hcp, a0_dhcp = initial_lattice_params(al_pct, co_pct, cr_pct)
        summary = read_results_summary(comp_dir)

        for struct in STRUCTURES:
            if struct == "FCC":
                a0 = a0_fcc
            elif struct == "HCP":
                a0 = a0_hcp
            else:
                a0 = a0_dhcp

            for T in TEMPERATURES:
                rec = summary.get((struct, T))
                a_final = final_a(struct, rec, a0) if rec is not None else None
                if rec is None:
                    print(f"NOTE: missing data for {name} {struct} at {T} K")

                rows.append(
                    {
                        "CompName": name,
                        "Al_at_pct": al_pct,
                        "Co_at_pct": co_pct,
                        "Cr_at_pct": cr_pct,
                        "Temp_K": T,
                        "Structure": struct,
                        "a_init_A": round(a0, 5),
                        "a_final_A": "" if a_final is None else round(a_final, 5),
                    }
                )

    return rows


def write_csv(rows, out_path: Path):
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "CompName",
                "Al_at_pct",
                "Co_at_pct",
                "Cr_at_pct",
                "Temp_K",
                "Structure",
                "a_init_A",
                "a_final_A",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"✓ Wrote {out_path}")


# ---------- ternary plotting ----------

def bary_to_xy(a, b, c):
    """
    Convert barycentric (Al, Co, Cr) to 2D coordinates.
    Vertices:
      Cr: (0, 0)      bottom-left
      Co: (1, 0)      bottom-right
      Al: (0.5, √3/2) top
    """
    SQ3_2 = math.sqrt(3) / 2.0
    A = np.array([0.5, SQ3_2])  # Al
    B = np.array([1.0, 0.0])    # Co
    C = np.array([0.0, 0.0])    # Cr
    return a * A + b * B + c * C


def draw_ternary_frame(ax):
    SQ3_2 = math.sqrt(3) / 2.0
    # Triangle
    tri_x = [0.0, 1.0, 0.5, 0.0]
    tri_y = [0.0, 0.0, SQ3_2, 0.0]
    ax.plot(tri_x, tri_y, "k-", linewidth=1.8)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, SQ3_2 + 0.03)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])

    # Edge labels
    ax.text(0.5, -0.08, "Cr fraction", ha="center", va="top", fontsize=10)
    ax.text(0.06, 0.5, "Al fraction", rotation=60,
            ha="center", va="center", fontsize=10)
    ax.text(0.94, 0.5, "Co fraction", rotation=-60,
            ha="center", va="center", fontsize=10)


def make_ternary_maps(df: pd.DataFrame):
    """Create ternary contour maps of a_final_A for each (Structure, Temp)."""

    df = df.copy()
    df["a_final_A"] = pd.to_numeric(df["a_final_A"], errors="coerce")
    df = df[df["a_final_A"].notna()]
    if df.empty:
        print("No valid a_final_A values for plotting.")
        return

    # Composition fractions and XY coordinates
    df["Al_frac"] = df["Al_at_pct"] / 100.0
    df["Co_frac"] = df["Co_at_pct"] / 100.0
    df["Cr_frac"] = df["Cr_at_pct"] / 100.0

    xy = np.vstack(
        [bary_to_xy(a, b, c) for a, b, c in
         zip(df["Al_frac"], df["Co_frac"], df["Cr_frac"])]
    )
    df["x"] = xy[:, 0]
    df["y"] = xy[:, 1]

    for struct in STRUCTURES:
        for T in TEMPERATURES:
            sub = df[(df["Structure"] == struct) & (df["Temp_K"] == T)]
            if len(sub) < 3:   # need at least a triangle
                print(f"Skip ternary map for {struct} at {T} K (not enough points).")
                continue

            x = sub["x"].to_numpy()
            y = sub["y"].to_numpy()
            z = sub["a_final_A"].to_numpy()

            tri = Triangulation(x, y)

            fig, ax = plt.subplots(figsize=(6, 5.5))
            contour = ax.tricontourf(tri, z, levels=12, cmap="plasma")

            draw_ternary_frame(ax)
            ax.set_title(f"{struct} lattice parameter a (Å) at {T} K")

            cbar = plt.colorbar(contour, ax=ax, pad=0.02)
            cbar.set_label("a (Å)")

            plt.tight_layout()
            fname = f"ternary_a_{struct}_{T}K.png"
            fig.savefig(fname, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"✓ Saved {fname}")


# ---------- main ----------

def main():
    rows = collect_rows()
    out_path = Path("lattice_all_compositions_all_temps.csv")
    write_csv(rows, out_path)

    df = pd.DataFrame(rows)
    make_ternary_maps(df)


if __name__ == "__main__":
    main()
