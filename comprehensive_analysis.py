#!/usr/bin/env python3
"""
Comprehensive Analysis Script for Al-Cr-Co Assignment
Generates all required plots and analysis for the report
Group 1: Al-Co-Cr System, Temperatures: 100K, 350K, 550K
"""

import os
from pathlib import Path
import numpy as np

class ComprehensiveAnalyzer:
    """Complete analysis including SFE, lattice parameters, and temperature effects"""
    
    def __init__(self):
        self.temperatures = [100, 350, 550]
        self.structures = ['FCC', 'HCP', 'DHCP']
        self.all_results = []
        
    def extract_composition_from_name(self, comp_name):
        """Extract Al, Co, Cr fractions from directory name"""
        parts = comp_name.split('_')
        al_frac = float(parts[1].replace('Al', '')) / 100.0
        co_frac = float(parts[2].replace('Co', '')) / 100.0
        cr_frac = float(parts[3].replace('Cr', '')) / 100.0
        return al_frac, co_frac, cr_frac
    
    def extract_results(self, comp_dir, structure, temp):
        """Extract all data from results_summary.txt"""
        summary_file = comp_dir / "results_summary.txt"
        
        if not summary_file.exists():
            return None
        
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
        return None
    
    def calculate_lattice_parameter(self, volume, natoms):
        """Calculate average lattice parameter from volume"""
        # For FCC: V = (n/4) * a^3, so a = (4V/n)^(1/3)
        a = (4 * volume / natoms) ** (1/3)
        return a
    
    def calculate_sfe(self, comp_dir, temp):
        """Calculate SFEs using DMLF model"""
        
        fcc_data = self.extract_results(comp_dir, 'FCC', temp)
        hcp_data = self.extract_results(comp_dir, 'HCP', temp)
        dhcp_data = self.extract_results(comp_dir, 'DHCP', temp)
        
        if not all([fcc_data, hcp_data, dhcp_data]):
            return None
        
        # Extract energies
        E_fcc = fcc_data['pe_per_atom']
        E_hcp = hcp_data['pe_per_atom']
        E_dhcp = dhcp_data['pe_per_atom']
        
        # Extract areas and atoms
        A_fcc = fcc_data['area']
        natoms_fcc = fcc_data['natoms']
        
        # Calculate lattice parameter
        lattice_param = self.calculate_lattice_parameter(fcc_data['volume'], natoms_fcc)
        
        # Energy differences
        delta_E_dhcp_fcc = (E_dhcp - E_fcc) * natoms_fcc
        delta_E_hcp_fcc = (E_hcp - E_fcc) * natoms_fcc
        
        # SFEs using DMLF (convert eV/A^2 to mJ/m^2)
        conversion = 16021.77
        
        gamma_ISF = 4 * delta_E_dhcp_fcc / A_fcc * conversion
        gamma_ESF = (delta_E_hcp_fcc + 2 * delta_E_dhcp_fcc) / A_fcc * conversion
        gamma_Twin = 2 * delta_E_dhcp_fcc / A_fcc * conversion
        
        return {
            'E_fcc': E_fcc,
            'E_hcp': E_hcp,
            'E_dhcp': E_dhcp,
            'lattice_param': lattice_param,
            'volume': fcc_data['volume'],
            'natoms': natoms_fcc,
            'gamma_ISF': gamma_ISF,
            'gamma_ESF': gamma_ESF,
            'gamma_Twin': gamma_Twin
        }
    
    def analyze_all_compositions(self):
        """Analyze all composition directories"""
        
        comp_dirs = sorted([d for d in Path('.').iterdir() 
                           if d.is_dir() and d.name.startswith('Comp')])
        
        if not comp_dirs:
            print("ERROR: No composition directories found!")
            return False
        
        print("="*70)
        print("Comprehensive Analysis - Al-Cr-Co System")
        print(f"Compositions: {len(comp_dirs)}")
        print(f"Temperatures: {self.temperatures}")
        print("="*70)
        
        for comp_dir in comp_dirs:
            al, co, cr = self.extract_composition_from_name(comp_dir.name)
            
            print(f"\nProcessing: {comp_dir.name}")
            print(f"  Composition: Al={al:.2f}, Co={co:.2f}, Cr={cr:.2f}")
            
            for temp in self.temperatures:
                results = self.calculate_sfe(comp_dir, temp)
                
                if results:
                    results['composition'] = comp_dir.name
                    results['temperature'] = temp
                    results['Al_frac'] = al
                    results['Co_frac'] = co
                    results['Cr_frac'] = cr
                    
                    self.all_results.append(results)
                    
                    print(f"  {temp}K: a={results['lattice_param']:.4f} A, "
                          f"ISF={results['gamma_ISF']:.2f} mJ/m^2")
                else:
                    print(f"  {temp}K: Missing data")
        
        return True
    
    def save_comprehensive_csv(self):
        """Save all results to comprehensive CSV"""
        
        if not self.all_results:
            print("No results to save!")
            return
        
        output_file = 'comprehensive_results.csv'
        
        with open(output_file, 'w') as f:
            # Header
            f.write("composition,temperature,Al_frac,Co_frac,Cr_frac,")
            f.write("E_fcc,E_hcp,E_dhcp,lattice_param,volume,natoms,")
            f.write("gamma_ISF,gamma_ESF,gamma_Twin\n")
            
            # Data
            for r in self.all_results:
                f.write(f"{r['composition']},{r['temperature']},")
                f.write(f"{r['Al_frac']:.4f},{r['Co_frac']:.4f},{r['Cr_frac']:.4f},")
                f.write(f"{r['E_fcc']:.8f},{r['E_hcp']:.8f},{r['E_dhcp']:.8f},")
                f.write(f"{r['lattice_param']:.8f},{r['volume']:.8f},{r['natoms']},")
                f.write(f"{r['gamma_ISF']:.4f},{r['gamma_ESF']:.4f},{r['gamma_Twin']:.4f}\n")
        
        print(f"\n✓ Comprehensive results saved to: {output_file}")
    
    def generate_analysis_summary(self):
        """Generate text summary for report"""
        
        output_file = 'analysis_summary.txt'
        
        with open(output_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("COMPREHENSIVE ANALYSIS SUMMARY\n")
            f.write("Al-Cr-Co Ternary System (Group 1)\n")
            f.write("="*70 + "\n\n")
            
            # Temperature effect on lattice parameters
            f.write("1. LATTICE PARAMETER ANALYSIS\n")
            f.write("-"*70 + "\n")
            
            for temp in self.temperatures:
                temp_data = [r for r in self.all_results if r['temperature'] == temp]
                lattice_params = [r['lattice_param'] for r in temp_data]
                
                f.write(f"\nAt {temp}K:\n")
                f.write(f"  Average lattice parameter: {np.mean(lattice_params):.4f} A\n")
                f.write(f"  Range: {np.min(lattice_params):.4f} - {np.max(lattice_params):.4f} A\n")
                f.write(f"  Std deviation: {np.std(lattice_params):.4f} A\n")
            
            # SFE analysis
            f.write("\n\n2. STACKING FAULT ENERGY ANALYSIS\n")
            f.write("-"*70 + "\n")
            
            for sfe_type in ['gamma_ISF', 'gamma_ESF', 'gamma_Twin']:
                f.write(f"\n{sfe_type}:\n")
                
                for temp in self.temperatures:
                    temp_data = [r for r in self.all_results if r['temperature'] == temp]
                    sfe_values = [r[sfe_type] for r in temp_data]
                    
                    f.write(f"  At {temp}K:\n")
                    f.write(f"    Average: {np.mean(sfe_values):.2f} mJ/m^2\n")
                    f.write(f"    Range: {np.min(sfe_values):.2f} - {np.max(sfe_values):.2f} mJ/m^2\n")
            
            # Composition effects
            f.write("\n\n3. COMPOSITION EFFECTS (at 350K)\n")
            f.write("-"*70 + "\n")
            
            temp_350_data = [r for r in self.all_results if r['temperature'] == 350]
            
            # Find extremes
            max_isf = max(temp_350_data, key=lambda x: x['gamma_ISF'])
            min_isf = min(temp_350_data, key=lambda x: x['gamma_ISF'])
            
            f.write(f"\nHighest ISF: {max_isf['composition']}\n")
            f.write(f"  Al={max_isf['Al_frac']:.2f}, Co={max_isf['Co_frac']:.2f}, Cr={max_isf['Cr_frac']:.2f}\n")
            f.write(f"  ISF = {max_isf['gamma_ISF']:.2f} mJ/m^2\n")
            
            f.write(f"\nLowest ISF: {min_isf['composition']}\n")
            f.write(f"  Al={min_isf['Al_frac']:.2f}, Co={min_isf['Co_frac']:.2f}, Cr={min_isf['Cr_frac']:.2f}\n")
            f.write(f"  ISF = {min_isf['gamma_ISF']:.2f} mJ/m^2\n")
            
            f.write("\n" + "="*70 + "\n")
        
        print(f"✓ Analysis summary saved to: {output_file}")
    
    def print_table_for_latex(self):
        """Generate LaTeX table code"""
        
        output_file = 'latex_tables.tex'
        
        with open(output_file, 'w') as f:
            f.write("% Table 1: Pure element properties at 350K\n")
            f.write("\\begin{table}[h]\n")
            f.write("\\centering\n")
            f.write("\\caption{Calculated properties of pure elements at 350K}\n")
            f.write("\\begin{tabular}{lccc}\n")
            f.write("\\hline\n")
            f.write("Element & Lattice Parameter (\\AA) & Energy (eV/atom) & SFE$_{ISF}$ (mJ/m$^2$) \\\\\n")
            f.write("\\hline\n")
            
            # Pure elements
            pure_comps = ['Comp01', 'Comp02', 'Comp03']
            elements = ['Al', 'Co', 'Cr']
            
            for comp, elem in zip(pure_comps, elements):
                data = [r for r in self.all_results 
                       if r['composition'].startswith(comp) and r['temperature'] == 350]
                if data:
                    d = data[0]
                    f.write(f"{elem} & {d['lattice_param']:.4f} & {d['E_fcc']:.4f} & {d['gamma_ISF']:.2f} \\\\\n")
            
            f.write("\\hline\n")
            f.write("\\end{tabular}\n")
            f.write("\\label{tab:pure_elements}\n")
            f.write("\\end{table}\n\n")
            
            # Temperature effect table
            f.write("% Table 2: Temperature effect on equiatomic composition\n")
            f.write("\\begin{table}[h]\n")
            f.write("\\centering\n")
            f.write("\\caption{Temperature effect on Al$_{33}$Co$_{34}$Cr$_{33}$ alloy}\n")
            f.write("\\begin{tabular}{lcccc}\n")
            f.write("\\hline\n")
            f.write("T (K) & $a$ (\\AA) & $\\gamma_{ISF}$ & $\\gamma_{ESF}$ & $\\gamma_{Twin}$ \\\\\n")
            f.write("      &           & (mJ/m$^2$) & (mJ/m$^2$) & (mJ/m$^2$) \\\\\n")
            f.write("\\hline\n")
            
            # Equiatomic composition (Comp10)
            for temp in self.temperatures:
                data = [r for r in self.all_results 
                       if 'Al33_Co34_Cr33' in r['composition'] and r['temperature'] == temp]
                if data:
                    d = data[0]
                    f.write(f"{temp} & {d['lattice_param']:.4f} & {d['gamma_ISF']:.2f} & ")
                    f.write(f"{d['gamma_ESF']:.2f} & {d['gamma_Twin']:.2f} \\\\\n")
            
            f.write("\\hline\n")
            f.write("\\end{tabular}\n")
            f.write("\\label{tab:temp_effect}\n")
            f.write("\\end{table}\n")
        
        print(f"✓ LaTeX tables saved to: {output_file}")


def main():
    """Main analysis function"""
    
    analyzer = ComprehensiveAnalyzer()
    
    # Run comprehensive analysis
    if analyzer.analyze_all_compositions():
        # Save results
        analyzer.save_comprehensive_csv()
        analyzer.generate_analysis_summary()
        analyzer.print_table_for_latex()
        
        print("\n" + "="*70)
        print("Comprehensive analysis complete!")
        print("Generated files:")
        print("  - comprehensive_results.csv  (All numerical data)")
        print("  - analysis_summary.txt       (Text summary)")
        print("  - latex_tables.tex           (LaTeX table code)")
        print("\nNext: Install matplotlib/pandas and run plotting script")
        print("="*70)
    else:
        print("Analysis failed! Check if simulations are complete.")


if __name__ == "__main__":
    main()