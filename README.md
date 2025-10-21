# CISM Udine 2024 – Computational Mechanics Exercises

This repository contains MATLAB scripts and resources developed for the [CISM course on the acoustics of musical instruments](https://cism.it/en/activities/courses/C2404/), with a focus on numerical methods eigenvalue problems.

##  Overview

The project includes a collection of examples, exercises, and numerical experiments covering topics such as:
- Euler–Bernoulli beam theory
- Finite difference methods
- Numerical error convergence
- Polynomial interpolation and Vandermonde systems
- Poisson equation solvers

These materials are intended for educational and research purposes in computational mechanics, numerical methods, and applied mathematics.

##  Repository Structure

```
CISM_Udine2024-main/
│
├── LICENSE                     # License information
├── README.md                   # Project documentation (this file)
└── src/                        # MATLAB source files
    ├── EulerBernoulliBars/     # Beam vibration and frequency analysis
    └── FiniteDiffCoefficients/ # Finite difference stencils and Poisson solver examples
```

##  Requirements

- **MATLAB R2020b** or newer (or GNU Octave for open-source users)
- Basic knowledge of MATLAB syntax and numerical methods

##  Usage

1. Clone or download this repository:
   ```bash
   git clone https://github.com/yourusername/CISM_Udine2024.git
   ```
   or download the ZIP archive and extract it.

2. Open MATLAB and set the working directory to the project root:
   ```matlab
   cd('path/to/CISM_Udine2024-main/src')
   ```

3. Run any of the example scripts, e.g.:
   ```matlab
   BarVariableThick
   PoissonSolver
   ```

4. Modify parameters or code as desired to explore the numerical behavior.

##  Educational Purpose

This material supports lectures and workshops in **computational mechanics** and **numerical analysis**. It demonstrates practical implementations of classical PDE solvers and structural models used in engineering simulations.

##  License

This project is distributed under the terms of the license found in the `LICENSE` file.

##  Acknowledgments

Developed for the **CISM course: Physics of Musical Instruments Applied to Instrument Making, Udine, Italy, May 2024**.
