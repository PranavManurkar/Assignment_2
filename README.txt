Assignment 2 – Al–Co–Cr Stacking Fault Energies and Lattice Parameters
=====================================================================

Group  : 1
Course : MM 309/309N – Computational Methods for Materials
System : Al–Co–Cr ternary alloys
Temperatures : 100 K, 350 K, 550 K
Structures : FCC, HCP, DHCP


0. PROJECT OVERVIEW
-------------------

This mini–project studies how composition and temperature affect:

  - Lattice parameters of FCC, HCP and DHCP Al–Co–Cr alloys
  - Intrinsic, extrinsic and twin stacking fault energies
  - Phase stability between FCC and HCP
  - Comparison of DMLF-based fault energies with a simple axial Ising model

The workflow starts from LAMMPS simulations for 21 compositions, collects all
outputs into CSV files, post-processes them with Python to generate:
  - Tables of lattice parameters
  - Line plots (lattice vs. temperature, SFE vs. temperature)
  - Ternary contour maps (SFE, phase stability)
  - DMLF vs. Ising comparison graphs
and then uses these figures in three separate LaTeX reports.



1. REPOSITORY STRUCTURE
-----------------------

High-level structure (names may vary slightly):

  Assignment_2/
    ├── Comp01_Al100_Co00_Cr00/
    ├── Comp02_Al75_Co25_Cr00/
    ├── ...
    ├── Comp21_Al00_Co00_Cr100/
    │     ├── log_* (LAMMPS logs)
    │     ├── results_summary.txt
    │     └── other intermediate files
    │
    ├── lattice_all.py
    ├── SFE_processing_scripts/ (optional, e.g. compute_SFE_AlCoCr.py)
    ├── compare_DMLF_Ising_ISF.py
    ├── compare_DMLF_Ising_ESF.py
    ├── compare_DMLF_Ising_Twin.py
    ├── lattice_all_compositions_all_temps.csv
    ├── SFE_results_AlCrCo.csv
    ├── figures/ (all PNG plots used in reports)
    │     ├── FCC_a_vs_T_all_compositions.png
    │     ├── HCP_a_vs_T_all_compositions.png
    │     ├── DHCP_a_vs_T_all_compositions.png
    │     ├── gamma_ISF_vs_Temperature.png
    │     ├── gamma_ESF_vs_Temperature.png
    │     ├── gamma_Twin_vs_Temperature.png
    │     ├── ternary_SFE_350K.png
    │     ├── phase_stability_AlCoCr_350K.png
    │     ├── DMLF_vs_Ising_ISF_100K.png / 350K / 550K
    │     ├── DMLF_vs_Ising_ESF_100K.png / 350K / 550K
    │     ├── DMLF_vs_Ising_Twin_100K.png / 350K / 550K
    │     └── (any extra ternary or contour plots)
    │
    ├── Report          
    └── README.txt      


2. WORKFLOW STEP-BY-STEP
------------------------


2.1. LAMMPS SIMULATIONS AND RAW DATA
------------------------------------

(1) Composition and structure setup
    - For each of the 21 compositions, a folder named
      CompXX_AlYY_CoZZ_CrWW was created.
      Example: Comp01_Al100_Co00_Cr00, Comp13_Al33_Co34_Cr33, etc.
    - Within each folder, three supercells were generated:
      FCC, HCP and DHCP.
    - Atomic species were distributed randomly according to Al/Co/Cr
      atomic percentages.

(2) Initial lattice parameter (Vegard’s law)
    - FCC initial lattice parameter estimated via:
      a_FCC^0 = x_Al a_Al + x_Co a_Co + x_Cr a_Cr
    - HCP and DHCP in-plane parameter derived from:
      a_HCP^0 = a_DHCP^0 = a_FCC^0 / sqrt(2)

(3) LAMMPS simulation protocol
    - Potential: Sharifi–Wick MEAM potential (CrCoAl.meam, library.meam)
      from NIST.
    - Ensemble: NPT (periodic boundary conditions).
    - Timestepping:
        * Time step Δt = 1 fs
        * Equilibration: ~20,000 MD steps
        * Production: additional 20,000–30,000 steps
        * Averages taken over the last 10,000–15,000 steps
    - Temperatures: 100 K, 350 K, 550 K
    - Structures: FCC, HCP, DHCP

(4) Summary extraction
    - For each CompXX folder, a script or post-processing routine collected:
        * natoms
        * final volume and box lengths (lx, ly, lz)
        * per-atom energies for FCC, HCP, DHCP
      into a file: results_summary.txt

This gives the raw data needed to compute:
  - Equilibrium lattice parameter a for each structure and T
  - DMLF stacking fault energies and twin boundary energy


2.2. LATTICE PARAMETERS: DATA COLLECTION AND PLOTS
--------------------------------------------------

(5) Script: lattice_all.py
    - Purpose:
        * Scan all Comp* folders.
        * Parse composition from folder names (Al% / Co% / Cr%).
        * Read results_summary.txt in each folder.
    - Composition parsing:
        * Uses a regex on directory names (AlXX_CoYY_CrZZ).
        * Checks that Al+Co+Cr = 100 (issues a warning if not).
    - Initial lattice parameters:
        * Computes a_FCC^0 using Vegard’s law.
        * Computes a_HCP^0 = a_DHCP^0 = a_FCC^0 / sqrt(2).

(6) Computing final lattice parameter a
    - For FCC:
        * From relaxed volume:
          V = (natoms / 4) * a^3 → a = (4V / natoms)^(1/3)
    - For HCP/DHCP:
        * Uses lx ~ n_x * a
        * Estimates n_x from lx / a_0 and then computes a = lx / n_x

(7) Writing master lattice CSV
    - Output file: lattice_all_compositions_all_temps.csv
    - Columns:
        * CompName
        * Al_at_pct, Co_at_pct, Cr_at_pct
        * Temp_K
        * Structure (FCC, HCP, DHCP)
        * a_init_A (Vegard guess)
        * a_final_A (relaxed lattice parameter)

(8) Lattice plots
    - Using pandas + matplotlib:
        * FCC_a_vs_T_all_compositions.png
        * HCP_a_vs_T_all_compositions.png
        * DHCP_a_vs_T_all_compositions.png
    - Each plot:
        * x-axis: Temperature (K)
        * y-axis: a_final_A (Å)
        * One line per composition (labelled by composition).
    - These figures are used to show:
        * Thermal expansion behaviour.
        * Relative spread of a across compositions.


2.3. DMLF STACKING FAULT AND TWIN ENERGIES
------------------------------------------

(9) Generating SFE_results_AlCrCo.csv
    - For each composition, temperature and structure:
      * Read energies E_FCC, E_HCP, E_DHCP and A_FCC from processed data.
    - Use the DMLF formulas:
        - γ_ISF  = 4 (E_DHCP – E_FCC) / A_FCC
        - γ_ESF  = (E_HCP + 2 E_DHCP – 3 E_FCC) / A_FCC
        - γ_Twin = 2 (E_DHCP – E_FCC) / A_FCC
    - Store results in SFE_results_AlCrCo.csv with columns like:
        * Temp_K
        * composition identifier
        * E_FCC_eV, E_HCP_eV, E_DHCP_eV
        * A_FCC_m2, natoms
        * gamma_ISF_mJm2, gamma_ESF_mJm2, gamma_Twin_mJm2
        * delta_E_hcp_fcc, delta_E_dhcp_fcc, etc.

(10) SFE vs. temperature
    - From SFE_results_AlCrCo.csv, build:
        * gamma_ISF_vs_Temperature.png
        * gamma_ESF_vs_Temperature.png
        * gamma_Twin_vs_Temperature.png
    - Each graph:
        * x-axis: Temperature (K)
        * y-axis: corresponding SFE (mJ/m²)
        * Multiple lines, one per composition.

(11) Ternary maps at 350 K
    - Using python + matplotlib.tri:
        * Read γ_ISF, γ_ESF, γ_Twin at 350 K.
        * Convert compositions (Al, Co, Cr) into ternary coordinates.
        * Plot interpolated contour maps:
          ternary_SFE_350K.png (combined panel, or separate maps).
    - These maps highlight:
        * Cr-rich region → low SFE, low twin energy.
        * Al/Co-rich regions → higher SFE and twin energies.


2.4. PHASE STABILITY MAP
------------------------

(12) ΔE = E_HCP – E_FCC
    - From SFE_results_AlCrCo.csv, compute:
        ΔE = E_HCP – E_FCC
    - Use ternary plotting to map ΔE across composition space:
        * phase_stability_AlCoCr_350K.png

(13) Interpretation
    - ΔE < 0 → FCC more stable.
    - ΔE > 0 → HCP more stable.
    - The phase-stability map is compared against SFE and twin-energy
      maps to link low SFE / low γ_Twin regions to proximity of HCP
      stability.


2.5. AXIAL ISING (AIM) COMPARISON WITH DMLF
-------------------------------------------

(14) Motivation
    - Adapt the idea of a 2D Ising model (spin variables, J coupling)
      to a 1D axial model representing stacking along ⟨111⟩.
    - Model perfect FCC and HCP as two “spin” states.
    - Use E_HCP – E_FCC and area A_FCC to define an effective coupling.

(15) For ISF (compare_DMLF_Ising_ISF.py)
    - Input: SFE_results_AlCrCo.csv
    - Steps:
        * For each composition and T:
          - Extract ΔE = E_HCP – E_FCC and A_FCC.
          - Build a simple Ising-like estimate of γ_ISF.
        * Compare DMLF γ_ISF vs. Ising γ_ISF on scatter plots:
          - DMLF_vs_Ising_ISF_100K.png
          - DMLF_vs_Ising_ISF_350K.png
          - DMLF_vs_Ising_ISF_550K.png
    - Interpretation:
        * Points near y = x → good agreement.
        * Systematic offsets show limits of the simple Ising model.

(16) For ESF (compare_DMLF_Ising_ESF.py)
    - Similar logic, but targeting γ_ESF:
        * Uses an AIM-based relation between ISF and ESF.
        * Produces plots:
          - DMLF_vs_Ising_ESF_100K.png
          - DMLF_vs_Ising_ESF_350K.png
          - DMLF_vs_Ising_ESF_550K.png
    - Used to check whether the simplified model captures
      the main ESF trends.

(17) For Twin Boundary Energy (compare_DMLF_Ising_Twin.py)
    - Treat a coherent twin as two intrinsic faults back-to-back.
    - Use effective Ising parameters to approximate γ_Twin.
    - Scatter plots:
        * DMLF_vs_Ising_Twin_100K.png
        * DMLF_vs_Ising_Twin_350K.png
        * DMLF_vs_Ising_Twin_550K.png
    - Compare AIM vs DMLF twin energies and analyse deviations.

Each comparison script:
  - Loads the SFE CSV.
  - Computes AIM predictions.
  - Generates scatter plots with a y = x reference line.
  - Uses different colours / labels per composition to reduce clutter.


2.6. REPORT GENERATION 
-----------------------------------------------

(18) Report  – overview 
    - LaTeX file: Report.tex (title: “Stacking Fault Energies and Lattice
      Parameters in Al–Co–Cr Ternary Alloys”).
    - Focus:
        * Overall lattice parameter trends.
        * ISF as primary SFE.
        * Ternary SFE maps at 350 K.
        * Phase stability map.
        * Detailed comparison: DMLF vs. AIM for ISF.

Report :
  - Use the same base methodology section (potential, LAMMPS, DMLF).
  - Share core figures (lattice vs. T, ternary maps, phase stability).
  - Emphasise different parts of the data to make the documents
    distinct and non-overlapping in focus.


3. HOW TO RE-RUN THE WORKFLOW
------------------------------

Prerequisites:
  - LAMMPS installed (for running MD, if raw simulations are repeated).
  - Python 3.x with:
      * numpy
      * pandas
      * matplotlib
      * matplotlib.tri
  - MEAM potential files from NIST (CrCoAl.meam, library.meam).

Typical workflow from Assignment_2 directory:

  1) Run or reuse LAMMPS simulations in each CompXX_* folder.
     Make sure each folder has a valid results_summary.txt.

  2) Collect lattice parameters:
       python lattice_all.py
     → creates lattice_all_compositions_all_temps.csv
       and saves FCC/HCP/DHCP a(T) plots.

  3) Compute SFEs (if not precomputed):
       python <SFE_processing_script>.py
     → generates SFE_results_AlCrCo.csv

  4) Produce SFE vs T plots and ternary maps:
       python <plot_sfe_vs_T.py>
       python <plot_ternary_sfe_350K.py>

  5) Phase stability map:
       python <phase_stability_ternary.py>

  6) Ising comparisons:
       python compare_DMLF_Ising_ISF.py
       python compare_DMLF_Ising_ESF.py
       python compare_DMLF_Ising_Twin.py

  7) Compile the three LaTeX reports:
       pdflatex Report_ISF.tex
       pdflatex Report_ESF.tex
       pdflatex Report_Twin.tex



4. GROUP MEMBERS AND CONTRIBUTIONS
----------------------------------

Group Members:
  - Mohit Garhewal    (Roll No: 230005028)
  - Pranav Munurkar   (Roll No: 230005026)
  - Ansh Kyal         (Roll No: 230005005)
  - Nawed Ashraf      (Roll No: 230005033)

Individual contributions:

  • Mohit Garhewal (230005028)
    - Led the analysis of the effect of alloying and temperature on the
      lattice parameter.
    - Wrote and refined the lattice-related Python scripts (e.g. lattice_all.py).
    - Collected and cleaned the lattice data from all compositions and temperatures.
    - Generated the full set of lattice parameter plots (a vs T, composition trends).
    - Implemented and analysed the comparison between ISF energies from DMLF
      and the axial Ising model (DMLF_vs_Ising_ISF plots).
    - Helped organise the repository structure and overall workflow documentation.

  • Pranav Munurkar (230005026)
    - Initiated and designed the full stacking fault energy calculation
      strategy for each composition, temperature and structure.
    - Wrote and optimised the LAMMPS input scripts used to generate
      the dataset (FCC/HCP/DHCP simulations, NPT runs).
    - Ran and monitored the MD simulations, ensuring convergence and
      consistency of results.
    - Processed outputs into SFE_results_AlCrCo.csv and produced the
      SFE vs temperature curves.
    - Generated and interpreted the comparison between ESF energies from
      DMLF and the axial Ising model (DMLF_vs_Ising_ESF plots).

  • Ansh Kyal (230005005)
    - Assisted in preparing and polishing the graphs for all sections
      (including lattice and SFE plots).
    - Led the comparison for twin boundary energies between DMLF and the
      axial Ising approach (DMLF_vs_Ising_Twin plots).
    - Worked on report generation and editing, particularly the twin-focused
      report and integrating figures into LaTeX.

  • Nawed Ashraf (230005033)
    - Took the lead in assembling the final reports from the processed data.
    - Integrated all tables and figures into the LaTeX documents and ensured
      formatting consistency (alignment of figures, captions, references).

End of README
