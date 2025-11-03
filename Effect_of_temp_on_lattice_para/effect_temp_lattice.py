# lattice_all.py
#
# Collect lattice parameters for EVERY composition and EVERY temperature
# from all CompXX_Al*_Co*_Cr* folders and write a single CSV:
#
#   lattice_all_compositions_all_temps.csv
#
# Then generate plots to show ONLY:
#   Effect of temperature on lattice parameter for each structure,
#   with each composition explicitly mentioned (legend).
#
# Run from the Assignment_2 directory:
#   python lattice_all.py

from pathlib import Path
import math
import re
import csv

import pandas as pd
import matplotlib.pyplot as plt

# Temperatures and structures used in the simulations
TEMPERATURES = [100, 350, 550]
STRUCTURES = ["FCC", "HCP", "DHCP"]

# Reference pure FCC lattice parameters (Å)
# HCP/DHCP initial a is taken as a_fcc / sqrt(2)
PURE_FCC_A = {"Al": 4.05, "Co": 3.544, "Cr": 2.885}

# ---------- helpers ----------

# Matches the Al/Co/Cr percentages in folder names like Comp01_Al100_Co00_Cr00
name_pat = re.compile(r"Al(?P<Al>\d+)_Co(?P<Co>\d+)_Cr(?P<Cr>\d+)")


def parse_comp_dir_name(name: str):
    """
    Works for:
      Comp01_Al100_Co00_Cr00
      Comp13_Al33_Co34_Cr33
      Comp21_Al60_Co20_Cr20
    Returns (Al%, Co%, Cr%).
    """
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
    """
    Vegard's law for FCC, then derive HCP/DHCP a from FCC.
    Returns (a_fcc0, a_hcp0, a_dhcp0) in Å.
    """
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


# ---------- main collection & CSV write ----------

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
            # choose correct initial a for this structure
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
                    # helpful but not fatal
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


# ---------- plotting functions ----------

def make_plots(df: pd.DataFrame):
    """
    Use the dataframe to generate:
      - For each structure: a_final vs Temp_K (one line per composition),
        with a legend mentioning each composition.
    """

    df = df.copy()
    df["a_final_A"] = pd.to_numeric(df["a_final_A"], errors="coerce")

    # Effect of temperature for each structure, with each composition labeled
    for struct in STRUCTURES:
        sub = df[(df["Structure"] == struct) & df["a_final_A"].notna()]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(9, 7))
        for comp_name, g in sub.groupby("CompName"):
            g_sorted = g.sort_values("Temp_K")
            # Build a concise composition label, e.g. "Al33Co34Cr33"
            comp_info = g_sorted.iloc[0]
            label = (
                f"Al{int(comp_info['Al_at_pct'])}"
                f"Co{int(comp_info['Co_at_pct'])}"
                f"Cr{int(comp_info['Cr_at_pct'])}"
            )
            ax.plot(
                g_sorted["Temp_K"],
                g_sorted["a_final_A"],
                marker="o",
                linewidth=1.2,
                alpha=0.8,
                label=label,
            )

        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Lattice parameter a (Å)")
        ax.set_title(
            f"{struct}: Effect of temperature on lattice parameter\n"
            "for all Al–Co–Cr compositions"
        )
        ax.grid(True, linestyle="--", alpha=0.5)

        # Legend listing each composition
        ax.legend(
            title="Composition (at.%)",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
            fontsize="small",
        )

        fig.tight_layout()
        fig.savefig(f"{struct}_a_vs_T_all_compositions.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"✓ Saved {struct}_a_vs_T_all_compositions.png")


# ---------- main ----------

def main():
    rows = collect_rows()
    out_path = Path("lattice_all_compositions_all_temps.csv")
    write_csv(rows, out_path)

    # Build DataFrame for plotting
    df = pd.DataFrame(rows)
    make_plots(df)


if __name__ == "__main__":
    main()
