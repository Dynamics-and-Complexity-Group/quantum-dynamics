# Dynamical System and Schrödinger Simulations for Quantum Tunneling in an Asymmetric Double-Well Potential
# Last updated on October 24, 2025

This repository contains Python codes and Jupyter notebooks accompanying the manuscript:

> **Swetamber Das** (ORCID: [0000-0002-2669-0842])  swetamber.p@srmap.edu.in
> **Arghya Dutta** (ORCID: [0000-0003-2116-6475])  arghya.d@srmap.edu.in
> Department of Physics, SRM University-AP, Amaravati 522240, Andhra Pradesh, India  

---

## 🧩 Manuscript Status

*This repository accompanies a manuscript currently under review.*  
Once the paper is accepted, this README will be updated to include the journal citation and DOI.

---

## 📁 Repository Structure

```
.
├── data/                              # Data files used in simulations
│
├── notebooks/                         # Jupyter notebooks for figure generation
│   └── [various notebooks]            # Reproduces all figures from the manuscript
│
├── dynamical-systems-model-simulations/
│   ├── flow_asym_double_well.py       # 4D dynamical system model (mean–variance evolution)
│   ├── time_series_left_to_right.ipynb# Time evolution: left → right well transition
│   └── time_series_right_to_left.ipynb# Time evolution: right → left well transition
│
├── schroedinger-simulations/
│   ├── left-to-right-tunneling.py     # Crank–Nicolson solver for left → right tunneling
│   └── right-to-left-tunneling.py     # Crank–Nicolson solver for right → left tunneling
│
├── requirements.txt                   # Python dependencies
└── README.md                          # Project description
```

---

## ▶️ How to Run

Clone the repository and install dependencies:

```bash
git clone https://github.com/Dynamics-and-Complexity-Group/quantum-dynamics.git
cd quantum-dynamics/asym_double-well-tunneling
pip install -r requirements.txt
```

To reproduce figures from the manuscript, open and execute the notebooks in the `notebooks/` directory.

---

## 🧠 Description

The repository contains two complementary numerical approaches to study tunneling in an asymmetric double-well potential:

1. **Dynamical Systems Approach**  
   Implements the reduced four-dimensional dynamical system derived from the Ehrenfest equations to model the time evolution of the mean and variance.

2. **Time-Dependent Schrödinger Simulations**  
   Solves the Schrödinger equation using the **Crank–Nicolson method** to visualize wavepacket tunneling between the two wells.

---

## 📘 License

This work is distributed under the **MIT License**.

---

## 🏷️ Keywords

Quantum tunneling, asymmetric double-well potential, Crank–Nicolson method,  
Ehrenfest dynamics, dynamical systems, wavepacket dynamics, nonlinear physics.
