# Dynamical System and Schrödinger Simulations for Quantum Tunneling in an Asymmetric Double-Well Potential

This repository contains Python codes and Jupyter notebooks accompanying the manuscript:

- Swetamber Das (ORCID: 0000-0002-2669-0842)  swetamber.p@srmap.edu.in
- Arghya Dutta (ORCID: 0000-0003-2116-6475)  arghya.d@srmap.edu.in
- Department of Physics, SRM University-AP, Amaravati 522240, Andhra Pradesh, India  

## Manuscript Status

This repository accompanies a manuscript currently under review. Once the paper is accepted, this README will be updated to include the journal citation and DOI.


## Repository Structure

```
.
├── data/                              # Data files used in simulations
│
├── notebooks/                         # Jupyter notebooks for figure generation
│   ├── double_well_plot.ipynb               # Generates Figure 1
│   ├── asym-double-well-fixed-point.ipynb   # Generates Figure 2
│   ├── time_series_plots_from_data.ipynb    # Generates Figures 3 and 4
│   ├── time_series_schroedinger_sims.ipynb  # Generates Figures 5 and 6
│   └── skewness_kurtosis_plots.ipynb        # Generates Figures 7 and 8
│
├── dynamical-systems-model-solution/
│   ├── flow_asym_double_well.py            # 4D dynamical system model (mean-variance evolution)
│   ├── time_series_left_to_right_well.ipynb # Time evolution: left -> right well transition
│   └── time_series_right_to_left_well.ipynb # Time evolution: right -> left well transition
│
├── schroedinger-simulations/
│   ├── left-to-right-tunneling.py     # Crank-Nicolson solver for left -> right tunneling
│   ├── right-to-left-tunneling.py     # Crank-Nicolson solver for right -> left tunneling
│   └── run_four_cases_together.py     # Runs the four Schrödinger tunneling cases together
│
├── requirements.txt                   # Python dependencies
└── README.md                          # Project description
```

## How to Run

Download this package and install the dependencies:

```bash
cd quantum-dynamics
pip install -r requirements.txt
```

To reproduce figures from the manuscript, open and execute the notebooks in the `notebooks/` directory. The supporting dynamical-systems notebooks are in `dynamical-systems-model-solution/`, and the Schrödinger simulation scripts are in `schroedinger-simulations/`.

## Description

The repository contains two complementary numerical approaches to study tunneling in an asymmetric double-well potential:

1. **Dynamical Systems Approach**  
   Implements the reduced four-dimensional dynamical system derived from the Ehrenfest equations to model the time evolution of the mean and variance.

2. **Time-Dependent Schrödinger Simulations**  
   Solves the Schrödinger equation using the **Crank–Nicolson method** to visualize wave packet tunneling between the two wells.


## License

This work is distributed under the **MIT License**.

## Citation

If you use this repository, please cite the associated manuscript. Once the paper is accepted, this README will be updated with the journal reference and DOI.

For archival citation, please cite the Zenodo record associated with the GitHub release you used.


## Keywords

Quantum tunneling, asymmetric double-well potential, Crank–Nicolson method,  
Ehrenfest dynamics, dynamical systems, wavepacket dynamics, nonlinear physics.
