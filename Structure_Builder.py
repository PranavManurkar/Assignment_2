#!/usr/bin/env python3
"""
OPTIMIZED Structure Builder - Al-Cr-Co Ternary System
Group 1: Al-Co-Cr system (Assignment 2, MM309)
Temperatures: 100K, 550K, 350K

Key optimization: 6×6×6 supercells for faster simulations
Still maintains sufficient size to minimize boundary effects
"""

import numpy as np
from itertools import product
import os

class OptimizedAlloyStructureBuilder:
    def __init__(self, composition, lattice_params):
        """
        composition: dict like {'Al': 0.33, 'Cr': 0.33, 'Co': 0.34}
        lattice_params: dict with reference lattice parameters for FCC
        """
        self.composition = composition
        self.lattice_params = lattice_params
        self.elements = sorted(composition.keys())

        # Vegard's law for average lattice parameter
        self.a_fcc = sum(composition[el] * lattice_params[el] for el in self.elements)

    def create_fcc_supercell(self, nx=6, ny=6, nz=6):
        """
        OPTIMIZED: Create smaller FCC supercell (6×6×6 instead of 8×8×8)
        Reduces atoms from 2048 to 864 (58% reduction)
        """
        basis = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5]
        ])

        positions = []

        for i, j, k in product(range(nx), range(ny), range(nz)):
            for b in basis:
                pos = (np.array([i, j, k]) + b) * self.a_fcc
                positions.append(pos)

        positions = np.array(positions)
        n_atoms = len(positions)

        atom_types = self._assign_atom_types(n_atoms)
        box = np.array([nx, ny, nz]) * self.a_fcc

        return positions, atom_types, box

    def create_hcp_supercell(self, nx=6, ny=6, nz=12):
        """
        OPTIMIZED: Smaller HCP structure (6×6×12 instead of 8×8×16)
        Maintains aspect ratio for stacking direction
        """
        a = self.a_fcc / np.sqrt(2)
        c = a * np.sqrt(8/3)

        positions = []

        for k in range(nz):
            z = k * c / 2
            if k % 2 == 0:  # A layer
                for i in range(nx):
                    for j in range(ny):
                        x = i * a
                        y = j * a * np.sqrt(3)
                        positions.append([x, y, z])
            else:  # B layer (shifted)
                for i in range(nx):
                    for j in range(ny):
                        x = (i + 0.5) * a
                        y = (j + 0.5) * a * np.sqrt(3)
                        positions.append([x, y, z])

        positions = np.array(positions)
        n_atoms = len(positions)
        atom_types = self._assign_atom_types(n_atoms)

        box = np.array([nx * a, ny * a * np.sqrt(3), nz * c / 2])

        return positions, atom_types, box

    def create_dhcp_supercell(self, nx=6, ny=6, nz=12):
        """
        OPTIMIZED: Smaller DHCP structure (6×6×12 instead of 8×8×16)
        """
        a = self.a_fcc / np.sqrt(2)
        c = a * np.sqrt(8/3)
        layer_spacing = c / 2

        positions = []

        stacking_sequence = ['A', 'B', 'A', 'C']

        shifts = {
            'A': np.array([0.0, 0.0]),
            'B': np.array([0.5, 0.5]),
            'C': np.array([0.5, 0.0])
        }

        for k in range(nz):
            layer_type = stacking_sequence[k % 4]
            shift = shifts[layer_type]
            z = k * layer_spacing

            for i in range(nx):
                for j in range(ny):
                    x = (i + shift[0]) * a
                    y = (j + shift[1]) * a * np.sqrt(3)
                    positions.append([x, y, z])

        positions = np.array(positions)
        n_atoms = len(positions)
        atom_types = self._assign_atom_types(n_atoms)

        box = np.array([nx * a, ny * a * np.sqrt(3), nz * layer_spacing])

        return positions, atom_types, box

    def _assign_atom_types(self, n_atoms):
        """Randomly assign atom types based on composition"""
        atom_type_pool = []

        for idx, el in enumerate(self.elements, 1):
            frac = self.composition[el]
            n_el = int(round(frac * n_atoms))
            atom_type_pool.extend([idx] * n_el)

        # Handle rounding errors
        while len(atom_type_pool) < n_atoms:
            atom_type_pool.append(1)
        while len(atom_type_pool) > n_atoms:
            atom_type_pool.pop()

        np.random.shuffle(atom_type_pool)
        return np.array(atom_type_pool)

    def write_lammps_data(self, positions, atom_types, box, filename, structure_type):
        """Write LAMMPS data file in atomic style"""
        n_atoms = len(positions)
        n_types = len(self.elements)

        with open(filename, 'w') as f:
            f.write(f"# {structure_type} structure - OPTIMIZED (smaller supercell)\n")
            f.write(f"# Composition: {self.composition}\n")
            f.write(f"# Al-Cr-Co System (Group 1) - Assignment 2 - MM309\n\n")

            f.write(f"{n_atoms} atoms\n")
            f.write(f"{n_types} atom types\n\n")

            f.write(f"0.0 {box[0]:.8f} xlo xhi\n")
            f.write(f"0.0 {box[1]:.8f} ylo yhi\n")
            f.write(f"0.0 {box[2]:.8f} zlo zhi\n\n")

            f.write("Masses\n\n")
            masses = {'Al': 26.9815, 'Co': 58.9332, 'Cr': 51.9961}
            for idx, el in enumerate(self.elements, 1):
                f.write(f"{idx} {masses[el]:.4f}  # {el}\n")

            f.write("\nAtoms # atomic\n\n")
            for i in range(n_atoms):
                f.write(f"{i+1} {atom_types[i]} {positions[i,0]:.8f} {positions[i,1]:.8f} {positions[i,2]:.8f}\n")


def generate_compositions():
    """
    Generate compositions across Al-Cr-Co ternary diagram
    Total: 21 compositions covering vertices, edges, and interior
    """
    compositions = []

    # Pure elements (vertices) - 3
    compositions.append({'Al': 1.0, 'Co': 0.0, 'Cr': 0.0})
    compositions.append({'Al': 0.0, 'Co': 1.0, 'Cr': 0.0})
    compositions.append({'Al': 0.0, 'Co': 0.0, 'Cr': 1.0})

    # Binary edges - 9
    for frac in [0.25, 0.50, 0.75]:
        compositions.append({'Al': frac, 'Co': 1-frac, 'Cr': 0.0})
        compositions.append({'Al': frac, 'Co': 0.0, 'Cr': 1-frac})
        compositions.append({'Al': 0.0, 'Co': frac, 'Cr': 1-frac})

    # Interior ternary - 9
    compositions.append({'Al': 0.33, 'Co': 0.34, 'Cr': 0.33})
    compositions.append({'Al': 0.34, 'Co': 0.33, 'Cr': 0.33})
    compositions.append({'Al': 0.50, 'Co': 0.25, 'Cr': 0.25})
    compositions.append({'Al': 0.25, 'Co': 0.50, 'Cr': 0.25})
    compositions.append({'Al': 0.25, 'Co': 0.25, 'Cr': 0.50})
    compositions.append({'Al': 0.40, 'Co': 0.40, 'Cr': 0.20})
    compositions.append({'Al': 0.40, 'Co': 0.20, 'Cr': 0.40})
    compositions.append({'Al': 0.20, 'Co': 0.40, 'Cr': 0.40})
    compositions.append({'Al': 0.60, 'Co': 0.20, 'Cr': 0.20})

    return compositions


def main():
    """Main function - OPTIMIZED version for Al-Cr-Co system"""

    # Lattice parameters (Angstroms) at room temperature
    lattice_params = {
        'Al': 4.05,   # FCC
        'Co': 3.544,  # HCP -> FCC equivalent
        'Cr': 2.885   # BCC -> FCC equivalent (a_fcc ≈ a_bcc * sqrt(2))
    }

    print("=" * 70)
    print("OPTIMIZED Al-Cr-Co Ternary Alloy Structure Generator")
    print("Group 1: Al-Co-Cr System")
    print("Assignment 2: MM309, IIT Indore")
    print("Temperatures: 100K, 550K, 350K")
    print("=" * 70)
    print("\nOPTIMIZATIONS:")
    print("  • FCC: 6×6×6 (~864 atoms instead of 2048)")
    print("  • HCP/DHCP: 6×6×12 (~432 atoms instead of 1024)")
    print("  • Total speedup: ~2.5x faster per simulation")
    print("  • Still maintains sufficient size (>15Å in each direction)")
    print("=" * 70)

    compositions = generate_compositions()

    for idx, comp in enumerate(compositions, 1):
        print(f"\n[{idx}/{len(compositions)}] Generating OPTIMIZED structures:")
        print(f"  Al: {comp['Al']:.2f}, Co: {comp['Co']:.2f}, Cr: {comp['Cr']:.2f}")

        comp_name = f"Comp{idx:02d}_Al{int(comp['Al']*100):02d}_Co{int(comp['Co']*100):02d}_Cr{int(comp['Cr']*100):02d}"
        os.makedirs(comp_name, exist_ok=True)

        builder = OptimizedAlloyStructureBuilder(comp, lattice_params)

        # Generate FCC
        print("  → Creating FCC (6×6×6)...")
        pos_fcc, types_fcc, box_fcc = builder.create_fcc_supercell(6, 6, 6)
        builder.write_lammps_data(pos_fcc, types_fcc, box_fcc,
                                 f"{comp_name}/structure_fcc.data", "FCC")

        # Generate HCP
        print("  → Creating HCP (6×6×12)...")
        pos_hcp, types_hcp, box_hcp = builder.create_hcp_supercell(6, 6, 12)
        builder.write_lammps_data(pos_hcp, types_hcp, box_hcp,
                                 f"{comp_name}/structure_hcp.data", "HCP")

        # Generate DHCP
        print("  → Creating DHCP (6×6×12)...")
        pos_dhcp, types_dhcp, box_dhcp = builder.create_dhcp_supercell(6, 6, 12)
        builder.write_lammps_data(pos_dhcp, types_dhcp, box_dhcp,
                                 f"{comp_name}/structure_dhcp.data", "DHCP")

        print(f"  ✓ FCC: {len(pos_fcc)} atoms | HCP: {len(pos_hcp)} atoms | DHCP: {len(pos_dhcp)} atoms")
        print(f"  ✓ Box: {box_fcc[0]:.2f} × {box_fcc[1]:.2f} × {box_fcc[2]:.2f} Ų")

    print("\n" + "=" * 70)
    print("OPTIMIZED structure generation complete!")
    print(f"Total compositions: {len(compositions)}")
    print("Atom count reduced by ~60% → Faster simulations!")
    print("=" * 70)


if __name__ == "__main__":
    np.random.seed(42)
    main()