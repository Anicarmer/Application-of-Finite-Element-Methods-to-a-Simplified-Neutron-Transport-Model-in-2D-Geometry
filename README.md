# Application of Finite Element Methods to a Simplified Neutron Transport Model in Two-Dimensional Geometry

## Overview

This repository contains the Python implementation developed for my Master's dissertation entitled:

**Application of Finite Element Methods to a Simplified Neutron Transport Model in Two-Dimensional Geometry**

The project presents a finite-element solver for a simplified, stationary, two-dimensional neutron transport equation. The implementation is used to validate the numerical method against an analytical solution, study mesh convergence, investigate the influence of the Breit–Wigner radiative capture cross-section, and analyze the effect of different transport directions on the neutron flux.

<img width="800" height="794" alt="FEM Approximation" src="https://github.com/user-attachments/assets/1a08eccb-80fc-4094-8f63-afc8b4037a7e" />
<img width="1536" height="802" alt="Convergence Study" src="https://github.com/user-attachments/assets/eee8eb2b-a0e2-4646-807d-0277b9420115" />
<img width="1536" height="802" alt="Absorption cross section" src="https://github.com/user-attachments/assets/7fae6bb6-ddf1-4b97-bf4b-a984307fcd10" />
<img width="1536" height="802" alt="Error norms vs energy" src="https://github.com/user-attachments/assets/e646ed60-4860-48a7-b979-0af12d4df14b" />
<img width="1536" height="802" alt="Maximum nodal value vs energy" src="https://github.com/user-attachments/assets/c439002d-99e1-454e-8327-402391161f19" />
<img width="1536" height="802" alt="FEM solution for varying transport direction" src="https://github.com/user-attachments/assets/34877771-9d55-4ffb-aa88-0c8cc9334887" />
<img width="1536" height="802" alt="Exact solution profiles for varying transport direction" src="https://github.com/user-attachments/assets/2d6b890a-2fa5-4b16-a2a4-3588a3f8b36b" />
<img width="1536" height="802" alt="Error norms vs transport direction" src="https://github.com/user-attachments/assets/74818b24-3c0b-46ed-81fa-c9a14d12af89" />



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
