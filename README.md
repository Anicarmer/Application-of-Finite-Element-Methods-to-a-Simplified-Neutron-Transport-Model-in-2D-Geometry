# Finite Element Approximation of a Simplified Two-Dimensional Neutron Transport Equation

## Overview

This repository contains the Python implementation developed for my Master's dissertation entitled:

**Finite Element Approximation of a Simplified Two-Dimensional Neutron Transport Equation**

The project presents a finite element solver for a simplified stationary neutron transport equation in two dimensions. The implementation is used to validate the numerical method against an analytical solution, study mesh convergence, investigate the influence of the Breit–Wigner radiative capture cross-section, and analyze the effect of different transport directions on the neutron flux.

---

## Features

* Finite element approximation using continuous piece-wise linear ((P_1)) triangular elements.
* Structured triangular mesh generation over the unit square.
* Assembly of the global transport and reaction matrices.
* Implementation of homogeneous inflow boundary conditions.
* Validation against an analytical solution obtained by the method of characteristics.
* Computation of:

  * Maximum error,
  * (L^2) error,
  * (H^1) seminorm error.
* Mesh convergence study.
* Breit–Wigner radiative capture cross-section analysis.
* Investigation of different neutron transport directions.
* Automatic generation of tables and figures used in the dissertation.

---

## Requirements

The code requires Python 3 and the following packages:

* NumPy
* Matplotlib

Install the dependencies using:

```bash
pip install numpy matplotlib
```

---

## Running the Code

Simply execute:

```bash
Neutron Transport Simulation Python.py
```

The program will:

1. Solve the reference transport problem.
2. Compute the error norms.
3. Perform the convergence study.
4. Evaluate the Breit–Wigner cross-section.
5. Perform the energy-dependent numerical study.
6. Analyze different transport directions.
7. Display all figures generated during the simulations.

---

## Repository Structure

```text
.
├── Neutron Transport Simulation Python.py    # Main Python implementation
├── README.md              # Project documentation
└── LICENSE                # License (optional)
```

---

## Dissertation

This implementation accompanies the Master's dissertation:

> **Finite Element Approximation of a Simplified Two-Dimensional Neutron Transport Equation**

Faculty of Sciences of El Jadida
Chouaïb Doukkali University
Master's Program in Computational Physics

---

## License

This project is released for academic and educational purposes. Please cite the corresponding dissertation if you use or reference this work.
