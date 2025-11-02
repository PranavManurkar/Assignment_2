#!/usr/bin/env python3
"""
Parallel LAMMPS Simulation Runner - Al-Cr-Co System
Runs 9 simulations simultaneously in each composition folder
"""

import subprocess
import os
from pathlib import Path
import time

LAMMPS_EXE = "lmp"  # Change to your LAMMPS executable
structures = ['FCC', 'HCP', 'DHCP']
temperatures = [100, 350, 550]

def run_simulations_parallel(comp_dir):
    """Run all 9 simulations in parallel for a composition"""
    
    print(f"\n{'='*70}")
    print(f"Processing: {comp_dir.name}")
    print(f"{'='*70}")
    
    os.chdir(comp_dir)
    
    # Check for required files
    if not Path("library.meam").exists():
        print(f"  ERROR: Missing library.meam")
        os.chdir("..")
        return 0
    if not Path("CrCoAl.meam").exists():
        print(f"  ERROR: Missing CrCoAl.meam")
        os.chdir("..")
        return 0
    
    print("Starting 9 parallel simulations...\n")
    
    processes = []
    sim_num = 0
    
    # Start all 9 simulations
    for structure in structures:
        for temp in temperatures:
            sim_num += 1
            input_file = f"in.{structure}_{temp}K.lammps"
            log_file = f"log.{structure}.{temp}K.lammps"
            
            if not Path(input_file).exists():
                print(f"  [{sim_num}/9] ERROR: Missing {input_file}")
                continue
            
            # Start LAMMPS process
            process = subprocess.Popen(
                [LAMMPS_EXE, "-in", input_file, "-log", log_file],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            processes.append(process)
            print(f"  [{sim_num}/9] Started: {structure}_{temp}K (PID: {process.pid})")
    
    print(f"\nWaiting for all {len(processes)} simulations to complete...")
    
    # Wait for all processes to finish
    for i, process in enumerate(processes, 1):
        process.wait()
        print(f"  Completed: {i}/{len(processes)}")
    
    print(f"All simulations in {comp_dir.name} completed!")
    
    os.chdir("..")
    return len(processes)


def main():
    """Main function"""
    
    print("="*70)
    print("Al-Cr-Co Parallel Simulation Runner")
    print("Running 9 simulations simultaneously per composition")
    print("="*70)
    
    # Find all composition directories
    comp_dirs = sorted([d for d in Path('.').iterdir() 
                       if d.is_dir() and d.name.startswith('Comp')])
    
    if not comp_dirs:
        print("\nERROR: No composition directories found!")
        return
    
    print(f"\nFound {len(comp_dirs)} composition directories")
    print(f"Total simulations: {len(comp_dirs) * 9}")
    
    start_time = time.time()
    total_sims = 0
    
    for i, comp_dir in enumerate(comp_dirs, 1):
        print(f"\n[{i}/{len(comp_dirs)}]")
        sims_run = run_simulations_parallel(comp_dir)
        total_sims += sims_run
    
    elapsed_time = time.time() - start_time
    
    print("\n" + "="*70)
    print("All Compositions Completed!")
    print("="*70)
    print(f"Total compositions processed: {len(comp_dirs)}")
    print(f"Total simulations run: {total_sims}")
    print(f"Total time: {elapsed_time/3600:.2f} hours")
    print("="*70)


if __name__ == "__main__":
    main()