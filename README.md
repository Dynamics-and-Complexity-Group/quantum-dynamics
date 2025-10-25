## Dynamical System and Schrödinger Simulations for Quantum Tunneling in an Asymmetric Double-Well Potential

This repository accompanies the manuscript titled “Dynamical system and Schrödinger simulations for quantum tunneling in an asymmetric double-well potential,” currently under review.

Once accepted, this README will be updated with the journal citation and DOI.

Last updated: October 25, 2025

**Authors:**  
> **Swetamber Das** (ORCID: [0000-0002-2669-0842]) — swetamber.p@srmap.edu.in  
> **Arghya Dutta** (ORCID: [0000-0003-2116-6475]) — arghya.d@srmap.edu.in  
> Department of Physics, SRM University-AP, Amaravati 522240, Andhra Pradesh, India.  

We welcome any comments or suggestions for improvement.

---

## Repository Structure

```
.
├── data/                              # Data files used in simulations
│
├── notebooks/                         # Jupyter notebooks for figure generation
│   └── [various notebooks]            # Reproduce all figures from the manuscript
│
├── dynamical-systems-model-simulations/
│   ├── flow_asym_double_well.py       # 4D dynamical system model (time evolution of mean and variance)
│   ├── time_series_left_to_right.ipynb# Time evolution: left → right well tunneling
│   └── time_series_right_to_left.ipynb# Time evolution: right → left well tunneling
│
├── schroedinger-simulations/
│   ├── left-to-right-tunneling.py     # Crank–Nicolson solver for left → right tunneling
│   └── right-to-left-tunneling.py     # Crank–Nicolson solver for right → left tunneling
│
├── requirements.txt                   # Python dependencies
└── README.md                          # Project description
```

---

## How to Run

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/Dynamics-and-Complexity-Group/quantum-dynamics.git
cd quantum-dynamics/asym_double-well-tunneling
pip install -r requirements.txt
```

To reproduce the figures from the manuscript, open and execute the notebooks in the `notebooks/` directory.

---

## Description

This repository presents two complementary numerical approaches to study quantum tunneling in an asymmetric double-well potential:

1. **Dynamical Systems Approach**  
   Implements a reduced four-dimensional dynamical system derived from the Ehrenfest equations to model the time evolution of the mean and variance of the wavepacket.

2. **Time-Dependent Schrödinger Simulations**  
   Solves the time-dependent Schrödinger equation using the **Crank–Nicolson method** to visualize wavepacket tunneling between the two wells.

---

## License

This project is distributed under the terms of the **MIT License**.

---

## Keywords

Quantum tunneling, asymmetric double-well potential, Crank–Nicolson method,  
Ehrenfest dynamics, dynamical systems, wavepacket dynamics, nonlinear physics.
