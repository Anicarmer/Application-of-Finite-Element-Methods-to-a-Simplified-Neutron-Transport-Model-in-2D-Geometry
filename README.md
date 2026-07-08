# Application of Finite Element Methods to a Simplified Neutron Transport Model in Two-Dimensional Geometry

## Overview

This repository contains the Python implementation developed for my Master's dissertation entitled:

**Application of Finite Element Methods to a Simplified Neutron Transport Model in Two-Dimensional Geometry**

The project presents a finite-element solver for a simplified, stationary, two-dimensional neutron transport equation. The implementation is used to validate the numerical method against an analytical solution, study mesh convergence, investigate the influence of the Breit–Wigner radiative capture cross-section, and analyze the effect of different transport directions on the neutron flux.

<img width="700" height="794" alt="fem_solution" src="https://github.com/user-attachments/assets/f91b27bb-ac74-4df1-964a-0fead3bab83f" />
<img width="1536" height="802" alt="Convergence Study" src="https://github.com/user-attachments/assets/a0f39425-976d-4b9b-bd19-1a9aa9612a1c" />
<img width="1536" height="802" alt="Absorption cross section" src="https://github.com/user-attachments/assets/8a28884d-3972-4b81-8a5a-2fe3a3c4611e" />
<img width="1536" height="802" alt="Error norms vs energy" src="https://github.com/user-attachments/assets/7fcd515d-a9cc-4d55-b26f-df0a647f0292" />
<img width="1536" height="802" alt="Maximum nodal value vs energy" src="https://github.com/user-attachments/assets/555bc620-4e85-49e2-8ad5-ef4e4a9cf484" />
<img width="1536" height="802" alt="FEM solution for varying transport direction" src="https://github.com/user-attachments/assets/19cef68b-e130-4d53-9702-75e08fd0e0fa" />
<img width="1536" height="802" alt="Exact solution profiles for varying transport direction" src="https://github.com/user-attachments/assets/4ea50e34-994e-4df8-bedc-7ab1433899f5" />
<img width="1536" height="802" alt="Error norms vs transport direction" src="https://github.com/user-attachments/assets/3737e7f7-af7a-422d-8683-4fc8cedca33a" />



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
* Pandas

Install the dependencies using:

```bash
pip install numpy matplotlib pandas
```

---

## Running the Code

Execute:

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

## Dissertation

This implementation accompanies the Master's dissertation:

> **Application of Finite Element Methods to a Simplified Neutron Transport Model in Two-Dimensional Geometry**

Faculty of Sciences of El Jadida
Chouaïb Doukkali University
Master's Program in Computational Physics

---

## Citation

@misc{haho2026github,
  author       = {Nizar Haho},
  title        = {Application of Finite Element to a Simplified Neutron Transport Model in Two-Dimensional Geometry},
  year         = {2026},
  howpublished = {\url{https://github.com/Anicarmer/Application-of-Finite-Element-Methods-to-a-Simplified-Neutron-Transport-Model-in-2D-Geometry}},
  note         = {Python implementation accompanying the Master's dissertation}
}

## License

This project is released for academic and educational purposes. Please cite the corresponding dissertation if you use or reference this work.
