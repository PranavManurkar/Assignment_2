#!/usr/bin/env python3
"""
LAMMPS Input File Generator - Al-Cr-Co System
Generates input files INSIDE each composition folder
Group 1: Temperatures 100K, 550K, 350K
"""

import os
from pathlib import Path

def generate_lammps_input(structure, temp, output_dir):
    """
    Generate LAMMPS input file for given structure and temperature
    
    Parameters:
    -----------
    structure : str
        Structure type ('FCC', 'HCP', or 'DHCP')
    temp : int
        Temperature in Kelvin
    output_dir : str
        Directory to save the input file
    """
    
    structure_lower = structure.lower()
    
    content = f"""# LAMMPS input script - Al-Cr-Co System
# Structure: {structure}, Temperature: {temp}K
# Group 1: Al-Co-Cr system - Assignment 2, MM309

# ============================================================
# VARIABLE DEFINITIONS
# ============================================================
variable        structure string {structure}
variable        temp equal {temp}
variable        structure_file string structure_{structure_lower}.data

# ============================================================
# INITIALIZATION
# ============================================================
units           metal
atom_style      atomic
boundary        p p p

# Read structure
read_data       ${{structure_file}}

# ============================================================
# INTERATOMIC POTENTIAL (MEAM)
# ============================================================
# Element order: Al Co Cr (alphabetical, matches structure file)
# Make sure library.meam and CrCoAl.meam are in the same directory
pair_style      meam
pair_coeff      * * library.meam Al Co Cr CrCoAl.meam Al Co Cr

# ============================================================
# NEIGHBOR SETTINGS
# ============================================================
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ============================================================
# COMPUTE DEFINITIONS
# ============================================================
variable        natoms equal atoms
variable        area_xy equal lx*ly

# ============================================================
# OUTPUT SETTINGS
# ============================================================
thermo          100
thermo_style    custom step temp pe ke etotal press vol lx ly lz pxx pyy pzz

log             log.${{structure}}.${{temp}}K.lammps

# ============================================================
# STAGE 1: ENERGY MINIMIZATION
# ============================================================
print           "============================================"
print           "STAGE 1: Energy Minimization"
print           "============================================"

min_style       cg
minimize        1.0e-8 1.0e-10 10000 100000

variable        pe_min equal pe
variable        pe_min_atom equal pe/v_natoms

print           "Minimized PE: ${{pe_min}} eV"
print           "Minimized PE/atom: ${{pe_min_atom}} eV/atom"

# ============================================================
# STAGE 2: NPT EQUILIBRATION
# ============================================================
print           "============================================"
print           "STAGE 2: NPT Equilibration at ${{temp}} K"
print           "============================================"

reset_timestep  0
velocity        all create ${{temp}} 87654 dist gaussian
timestep        0.001

fix             npt1 all npt temp ${{temp}} ${{temp}} 0.1 iso 0.0 0.0 1.0
run             30000
unfix           npt1

# ============================================================
# STAGE 3: NPT PRODUCTION RUN
# ============================================================
print           "============================================"
print           "STAGE 3: NPT Production Run at ${{temp}} K"
print           "============================================"

reset_timestep  0
fix             npt2 all npt temp ${{temp}} ${{temp}} 0.1 iso 0.0 0.0 1.0

variable        pe_step equal pe
variable        pe_atom_step equal pe/v_natoms
variable        vol_step equal vol
variable        lx_step equal lx
variable        ly_step equal ly
variable        lz_step equal lz

fix             ave_pe all ave/time 10 100 1000 v_pe_atom_step &
                file pe_vs_time.${{structure}}.${{temp}}K.dat

fix             ave_vol all ave/time 10 100 1000 v_vol_step v_lx_step v_ly_step v_lz_step &
                file vol_vs_time.${{structure}}.${{temp}}K.dat

run             50000

# ============================================================
# EXTRACT RESULTS
# ============================================================
variable        pe_avg equal f_ave_pe
variable        lx_final equal lx
variable        ly_final equal ly
variable        lz_final equal lz
variable        vol_final equal vol
variable        area_final equal v_lx_final*v_ly_final

print           ""
print           "============================================"
print           "FINAL RESULTS: ${{structure}} at ${{temp}} K"
print           "============================================"
print           "Number of atoms:      ${{natoms}}"
print           "Avg PE/atom:          ${{pe_avg}} eV/atom"
print           "Final volume:         ${{vol_final}} A^3"
print           "Final box:            ${{lx_final}} x ${{ly_final}} x ${{lz_final}} A"
print           "Interface area (xy):  ${{area_final}} A^2"
print           "============================================"

print           "${{structure}} ${{temp}} ${{natoms}} ${{pe_avg}} ${{vol_final}} ${{lx_final}} ${{ly_final}} ${{lz_final}} ${{area_final}}" &
                append results_summary.txt

unfix           npt2
unfix           ave_pe
unfix           ave_vol

print           "Simulation complete!"
"""
    
    filename = os.path.join(output_dir, f"in.{structure}_{temp}K.lammps")
    
    # Use UTF-8 encoding explicitly for Windows
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
    
    return filename


def main():
    """Generate all LAMMPS input files for Al-Cr-Co system"""
    
    # Group 1 parameters
    structures = ['FCC', 'HCP', 'DHCP']
    temperatures = [100, 550, 350]  # Group 1 temperatures
    
    print("=" * 70)
    print("LAMMPS Input File Generator - Al-Cr-Co System")
    print("Group 1: Al-Co-Cr")
    print("=" * 70)
    
    # Find all composition directories
    comp_dirs = sorted([d for d in Path('.').iterdir() 
                       if d.is_dir() and d.name.startswith('Comp')])
    
    if not comp_dirs:
        print("\nERROR: No composition directories found!")
        print("Please run Structure_Builder_AlCrCo.py first")
        return
    
    print(f"\nFound {len(comp_dirs)} composition directories")
    print(f"Structures: {', '.join(structures)}")
    print(f"Temperatures: {', '.join(map(str, temperatures))} K")
    print(f"Total files to generate: {len(comp_dirs)} x {len(structures)} x {len(temperatures)} = {len(comp_dirs) * len(structures) * len(temperatures)}")
    print("=" * 70)
    
    total_generated = 0
    
    for comp_dir in comp_dirs:
        print(f"\nProcessing: {comp_dir.name}")
        
        # Check if structure files exist
        missing_files = []
        for struct in structures:
            struct_file = comp_dir / f"structure_{struct.lower()}.data"
            if not struct_file.exists():
                missing_files.append(struct_file.name)
        
        if missing_files:
            print(f"   WARNING: Missing files: {', '.join(missing_files)}")
            print(f"   Skipping {comp_dir.name}")
            continue
        
        # Generate input files for this composition
        for structure in structures:
            for temp in temperatures:
                filename = generate_lammps_input(structure, temp, str(comp_dir))
                rel_path = os.path.relpath(filename)
                print(f"   Created: {rel_path}")
                total_generated += 1
    
    print("\n" + "=" * 70)
    print(f"Generation complete!")
    print(f"Total files created: {total_generated}")
    print("=" * 70)
    
    print("\nNext steps:")
    print("1. Copy library.meam and CrCoAl.meam to each Comp* folder")
    print("2. Run simulations in each folder")
    print("3. Analyze results using calculate_SFE.py")


if __name__ == "__main__":
    main()