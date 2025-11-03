#!/usr/bin/env python3
"""
mc_ternary_al_cocr.py

2D Monte Carlo (Metropolis) phase separation for a ternary alloy (Al - Co - Cr).
Produces lattice snapshots, energy vs step, interface-length vs step,
and optionally an animated GIF/MP4 of microstructure evolution.

Usage:
    python mc_ternary_al_cocr.py

Requirements:
    numpy, matplotlib
    optional: imageio (for GIF), tqdm (for progress bar)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import random
from math import exp
import os
try:
    from tqdm import trange
except Exception:
    def trange(x): return range(x)

# Optional: for saving GIF (pip install imageio)
try:
    import imageio
    HAVE_IMAGEIO = True
except Exception:
    HAVE_IMAGEIO = False

# ----------------------------
# Parameters (edit these)
# ----------------------------
L = 200                 # lattice linear dimension (200x200 -> 40k sites; reduce if slow)
fraction = [0.33, 0.33, 0.34]   # fractions for [Al, Co, Cr], must sum to 1.0
T = 1.0                 # temperature in energy units (k_B = 1)
n_steps = 300           # number of Monte Carlo steps (each step ~ L*L attempted swaps)
snapshot_interval = 30  # save figure every snapshot_interval steps
seed = 42               # RNG seed for reproducibility

# Pairwise bond energies (units arbitrary). Negative = favorable like-bond.
# Species indexing: 0 -> Al, 1 -> Co, 2 -> Cr
E = np.zeros((3, 3))
E[0, 0] = -1.0  # E_AlAl
E[1, 1] = -1.0  # E_CoCo
E[2, 2] = -1.0  # E_CrCr

# Cross interactions (positive means unlike bonds are unfavorable -> drives phase separation)
E[0, 1] = 0.2  # Al-Co
E[1, 0] = E[0, 1]
E[0, 2] = 0.45 # Al-Cr (more repulsive)
E[2, 0] = E[0, 2]
E[1, 2] = 0.05 # Co-Cr (slightly repulsive)
E[2, 1] = E[1, 2]

kB = 1.0  # Boltzmann constant in chosen energy units

# Output folder
outdir = "mc_outputs"
os.makedirs(outdir, exist_ok=True)
# ----------------------------

np.random.seed(seed)
random.seed(seed)

# initialize lattice with species labels 0,1,2 according to fractions
def init_lattice(L, fractions):
    assert len(fractions) == 3
    freqs = np.array(fractions) / np.sum(fractions)
    N = L * L
    counts = (freqs * N).astype(int)
    # adjust rounding: fill remaining randomly
    while counts.sum() < N:
        counts[np.random.randint(3)] += 1
    flat = np.zeros(N, dtype=np.int8)
    idx = 0
    for s in range(3):
        flat[idx:idx+counts[s]] = s
        idx += counts[s]
    np.random.shuffle(flat)
    return flat.reshape((L, L))

# neighbor offsets (4 nearest neighbors)
neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# local energy at site (sum of bond energies with neighbors)
def local_energy(lat, i, j):
    Ls = lat.shape[0]
    s = lat[i, j]
    e = 0.0
    for di, dj in neighbors:
        ni = (i + di) % Ls
        nj = (j + dj) % Ls
        e += E[s, lat[ni, nj]]
    return e

# total energy (each bond counted twice if we sum all local energies, so divide by 2)
def total_energy(lat):
    Ls = lat.shape[0]
    e = 0.0
    for i in range(Ls):
        for j in range(Ls):
            e += local_energy(lat, i, j)
    return 0.5 * e

# interface length (measure of unlike neighbor pairs) - count unique bonds only once
def interface_length(lat):
    Ls = lat.shape[0]
    count = 0
    for i in range(Ls):
        for j in range(Ls):
            s = lat[i, j]
            # consider only right and down neighbor to avoid double counting
            ni = i
            nj = (j + 1) % Ls
            if lat[ni, nj] != s:
                count += 1
            ni = (i + 1) % Ls
            nj = j
            if lat[ni, nj] != s:
                count += 1
    return count  # proportional to domain boundary length

# plotting utility
cmap = ListedColormap(['#ff7f0e', '#1f77b4', '#2ca02c'])  # Al: orange, Co: blue, Cr: green
species_labels = ['Al', 'Co', 'Cr']

def plot_lattice(lat, step, outdir=outdir):
    plt.figure(figsize=(5,5))
    plt.imshow(lat, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=2.5)
    plt.title(f"Microstructure (step {step})")
    plt.axis('off')
    # create legend
    patches = [plt.plot([],[], marker="s", ms=8, ls="", color=cmap(i))[0] for i in range(3)]
    plt.legend(patches, species_labels, loc='lower right', framealpha=0.9)
    fname = os.path.join(outdir, f"micro_step_{step:05d}.png")
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    return fname

# Monte Carlo simulation (Kawasaki swaps)
def run_mc(lat, T, n_steps, snapshot_interval=50):
    Ls = lat.shape[0]
    energy_trace = []
    iface_trace = []
    frame_files = []

    E_tot = total_energy(lat)
    iface = interface_length(lat)
    energy_trace.append(E_tot)
    iface_trace.append(iface)
    print(f"Initial total energy = {E_tot:.3f}, interface = {iface}")

    for step in trange(n_steps):
        # perform L*L attempted swaps per MC step
        n_trials = Ls * Ls
        for trial in range(n_trials):
            i1 = np.random.randint(Ls)
            j1 = np.random.randint(Ls)
            i2 = np.random.randint(Ls)
            j2 = np.random.randint(Ls)

            if lat[i1, j1] == lat[i2, j2]:
                continue

            # compute local energy before swap (only around the two sites)
            E_before = local_energy(lat, i1, j1) + local_energy(lat, i2, j2)

            # do swap
            lat[i1, j1], lat[i2, j2] = lat[i2, j2], lat[i1, j1]

            E_after = local_energy(lat, i1, j1) + local_energy(lat, i2, j2)
            dE = E_after - E_before

            if dE <= 0:
                # accept
                pass
            else:
                if random.random() > exp(-dE / (kB * T)):
                    # reject -> swap back
                    lat[i1, j1], lat[i2, j2] = lat[i2, j2], lat[i1, j1]
                else:
                    # accept
                    pass

        # diagnostics after the MC step
        if (step + 1) % snapshot_interval == 0 or step == n_steps - 1 or step == 0:
            E_tot = total_energy(lat)
            iface = interface_length(lat)
            energy_trace.append(E_tot)
            iface_trace.append(iface)
            print(f"Step {step+1:4d}: E = {E_tot:.3f}, interface = {iface}")
            fname = plot_lattice(lat, step+1)
            frame_files.append(fname)

    return lat, np.array(energy_trace), np.array(iface_trace), frame_files

# ----------------------------
# Main run
# ----------------------------
def main():
    lat = init_lattice(L, fraction)
    # sanity check
    print("Initial composition counts:", np.bincount(lat.flatten(), minlength=3))
    final_lat, energy_trace, iface_trace, frames = run_mc(lat, T, n_steps, snapshot_interval)

    # Save energy and interface traces
    steps_recorded = np.linspace(0, n_steps, len(energy_trace), dtype=int)
    np.savetxt(os.path.join(outdir, "energy_trace.txt"), np.column_stack([steps_recorded, energy_trace]),
               header="step total_energy")
    np.savetxt(os.path.join(outdir, "interface_trace.txt"), np.column_stack([steps_recorded, iface_trace]),
               header="step interface_length")

    # Plot energy and interface vs step
    plt.figure(figsize=(6,3))
    plt.plot(steps_recorded, energy_trace, marker='o')
    plt.xlabel("MC step")
    plt.ylabel("Total Energy (arb. units)")
    plt.title("Energy vs MC step")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "energy_vs_step.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(6,3))
    plt.plot(steps_recorded, iface_trace, marker='o', color='tab:red')
    plt.xlabel("MC step")
    plt.ylabel("Interface length (arb. units)")
    plt.title("Interface length vs MC step")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "interface_vs_step.png"), dpi=150)
    plt.close()

    # Optionally build GIF (requires imageio)
    if HAVE_IMAGEIO and len(frames) > 0:
        gif_path = os.path.join(outdir, "micro_evolution.gif")
        images = []
        for f in frames:
            images.append(imageio.imread(f))
        imageio.mimsave(gif_path, images, duration=0.6)
        print("Saved GIF:", gif_path)
    else:
        print("ImageIO not available or no frames; saved PNG frames in", outdir)

    print("Done. Example frame files:", frames[:5])

if __name__ == "__main__":
    main()
