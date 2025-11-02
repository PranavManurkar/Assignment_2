#!/usr/bin/env python3
"""
Stacking Fault Energy (SFE) Calculator - Al-Cr-Co System
Calculates SFEs using the DMLF (Diffuse Multi-layer Fault) model
Group 1: Al-Co-Cr system - Assignment 2, MM309

Reference: M.A. Charpagne et al., Acta Materialia 194 (2020) 224-235
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class SFECalculator:
    """Calculate stacking fault energies using DMLF model"""
    
    def __init__(self, composition_dir):
        """
        Initialize with a composition directory
        
        Parameters:
        -----------
        composition_dir : str
            Path to composition directory (e.g., 'Comp01_Al100_Co00_Cr00')
        """
        self.comp_dir = Path(composition_dir)
        self.comp_name = self.comp_dir.name
        self.results = {}
        
    def extract_results(self, structure, temp):
        """
        Extract energy and area from LAMMPS results
        
        Parameters:
        -----------
        structure : str
            Structure type ('FCC', 'HCP', 'DHCP')
        temp : int
            Temperature in Kelvin
            
        Returns:
        --------
        dict : Contains pe_per_atom and area
        """
        # Read from results_summary.txt
        summary_file = self.comp_dir / "results_summary.txt"
        
        if not summary_file.exists():
            raise FileNotFoundError(f"Results file not found: {summary_file}")
        
        with open(summary_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 9:
                    struct = parts[0]
                    temp_val = int(parts[1])
                    
                    if struct == structure and temp_val == temp:
                        return {
                            'natoms': int(parts[2]),
                            'pe_per_atom': float(parts[3]),
                            'volume': float(parts[4]),
                            'lx': float(parts[5]),
                            'ly': float(parts[6]),
                            'lz': float(parts[7]),
                            'area': float(parts[8])
                        }
        
        raise ValueError(f"Results not found for {structure} at {temp}K")
    
    def calculate_sfe(self, temp):
        """
        Calculate SFEs using DMLF model for a given temperature
        
        DMLF Equations:
        γ_ISF = 4(E_dhcp - E_fcc) / A_fcc
        γ_ESF = (E_hcp + 2E_dhcp - 3E_fcc) / A_fcc
        γ_Twin = 2(E_dhcp - E_fcc) / A_fcc
        
        Parameters:
        -----------
        temp : int
            Temperature in Kelvin
            
        Returns:
        --------
        dict : Contains all SFE values in mJ/m²
        """
        # Extract results for all structures
        fcc_data = self.extract_results('FCC', temp)
        hcp_data = self.extract_results('HCP', temp)
        dhcp_data = self.extract_results('DHCP', temp)
        
        # Extract energies (eV/atom)
        E_fcc = fcc_data['pe_per_atom']
        E_hcp = hcp_data['pe_per_atom']
        E_dhcp = dhcp_data['pe_per_atom']
        
        # Extract areas (Ų)
        A_fcc = fcc_data['area']
        
        # Calculate energy differences (eV)
        natoms_fcc = fcc_data['natoms']
        
        delta_E_dhcp_fcc = (E_dhcp - E_fcc) * natoms_fcc
        delta_E_hcp_fcc = (E_hcp - E_fcc) * natoms_fcc
        
        # Calculate SFEs using DMLF model
        # Convert from eV/Ų to mJ/m² (1 eV/Ų = 16021.77 mJ/m²)
        conversion = 16021.77
        
        gamma_ISF = 4 * delta_E_dhcp_fcc / A_fcc * conversion
        gamma_ESF = (delta_E_hcp_fcc + 2 * delta_E_dhcp_fcc) / A_fcc * conversion
        gamma_Twin = 2 * delta_E_dhcp_fcc / A_fcc * conversion
        
        results = {
            'temperature': temp,
            'E_fcc': E_fcc,
            'E_hcp': E_hcp,
            'E_dhcp': E_dhcp,
            'A_fcc': A_fcc,
            'natoms': natoms_fcc,
            'gamma_ISF': gamma_ISF,
            'gamma_ESF': gamma_ESF,
            'gamma_Twin': gamma_Twin,
            'delta_E_dhcp_fcc': delta_E_dhcp_fcc,
            'delta_E_hcp_fcc': delta_E_hcp_fcc
        }
        
        return results


def process_all_compositions(temperatures=[100, 350, 550]):
    """
    Process all composition directories and calculate SFEs
    
    Parameters:
    -----------
    temperatures : list
        List of temperatures to process
    """
    
    # Find all composition directories
    comp_dirs = sorted([d for d in Path('.').iterdir() 
                       if d.is_dir() and d.name.startswith('Comp')])
    
    if not comp_dirs:
        print("No composition directories found!")
        return None
    
    print("=" * 70)
    print(f"SFE Analysis - Al-Cr-Co System")
    print(f"Found {len(comp_dirs)} compositions")
    print(f"Temperatures: {temperatures} K")
    print("=" * 70)
    
    all_results = []
    
    for comp_dir in comp_dirs:
        print(f"\nProcessing: {comp_dir.name}")
        
        calculator = SFECalculator(comp_dir)
        
        for temp in temperatures:
            try:
                results = calculator.calculate_sfe(temp)
                results['composition'] = comp_dir.name
                all_results.append(results)
                
                print(f"  {temp}K - ISF: {results['gamma_ISF']:.2f}, "
                      f"ESF: {results['gamma_ESF']:.2f}, "
                      f"Twin: {results['gamma_Twin']:.2f} mJ/m²")
                
            except (FileNotFoundError, ValueError) as e:
                print(f"  {temp}K - Error: {e}")
    
    # Create DataFrame
    df = pd.DataFrame(all_results)
    
    # Save to CSV
    output_file = 'SFE_results_AlCrCo.csv'
    df.to_csv(output_file, index=False)
    print(f"\n✓ Results saved to: {output_file}")
    
    return df


def plot_sfe_results(df, output_dir='plots'):
    """
    Create plots of SFE vs temperature
    
    Parameters:
    -----------
    df : DataFrame
        Results DataFrame from process_all_compositions
    output_dir : str
        Directory to save plots
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot for each composition
    compositions = df['composition'].unique()
    
    for comp in compositions:
        comp_data = df[df['composition'] == comp].sort_values('temperature')
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(comp_data['temperature'], comp_data['gamma_ISF'], 
                'o-', label='γ_ISF', linewidth=2, markersize=8)
        ax.plot(comp_data['temperature'], comp_data['gamma_ESF'], 
                's-', label='γ_ESF', linewidth=2, markersize=8)
        ax.plot(comp_data['temperature'], comp_data['gamma_Twin'], 
                '^-', label='γ_Twin', linewidth=2, markersize=8)
        
        ax.set_xlabel('Temperature (K)', fontsize=12)
        ax.set_ylabel('Stacking Fault Energy (mJ/m²)', fontsize=12)
        ax.set_title(f'SFE vs Temperature - {comp}', fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/SFE_{comp}.png', dpi=300)
        plt.close()
    
    print(f"\n✓ Plots saved to: {output_dir}/")


def main():
    """Main function"""
    
    # Group 1 temperatures
    temperatures = [100, 350, 550]
    
    # Process all compositions
    df = process_all_compositions(temperatures)
    
    if df is not None:
        # Create plots
        plot_sfe_results(df)
        
        # Print summary statistics
        print("\n" + "=" * 70)
        print("Summary Statistics (all compositions, all temperatures)")
        print("=" * 70)
        print("\nIntrinsic SFE (γ_ISF):")
        print(df['gamma_ISF'].describe())
        print("\nExtrinsic SFE (γ_ESF):")
        print(df['gamma_ESF'].describe())
        print("\nTwin SFE (γ_Twin):")
        print(df['gamma_Twin'].describe())
        
        print("\n" + "=" * 70)
        print("Analysis complete!")
        print("=" * 70)


if __name__ == "__main__":
    main()