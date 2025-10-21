[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17406808.svg)](https://doi.org/10.5281/zenodo.17406808)

# CISM Udine 2024 – Physics of Musical Instruments


**Author:** Michele Ducceschi  
**Affiliation:** NEMUS Lab, University of Bologna  
**Email:** michele.ducceschi@unibo.it  
**Course:** CISM Course "Physics of Musical Instruments" (Udine, 2024)

This repository contains MATLAB scripts and resources developed for the course [Physics of Musical Instruments](https://cism.it/en/activities/courses/C2404/), held at Udine in May 2024 (by CISM). It focuses on numerical methods, eigenvalue problems, and their application to acoustics and instrument‑making.

---

## Overview

The project includes a collection of examples, exercises, and numerical experiments covering topics such as:

- Grid-based methods for differential equations in acoustics
- Polynomial interpolation (Lagrange, Vandermonde) for grid coefficient generation
- Poisson equation solvers on uniform and non-uniform grids
- Euler–Bernoulli beam theory (vibration and frequency analysis)  
- Use of the *MAGPIE* toolbox for orthotropic plate modal analysis and time domain integration

---

## Repository Structure

```
CISM_Udine2024/
│
├── LICENSE                    # License information
├── README.md                  # You are here
├── CISM_LectureNotes.pdf      # The reference lecture notes 
├── .gitmodules                # If using submodules (for magpie‑matlab)
├── external/                  # External dependencies
│   └── magpie‑matlab/         # The magpie‑matlab toolbox (see Dependencies)
└── src/                       # MATLAB source files
    ├── EulerBernoulliBars/    # Beam vibration and frequency analysis
    ├── FiniteDiffCoefficients/ # Finite difference stencils and Poisson solver examples
    └── …                      # (other folders as required)
```

---

## Dependencies

- **MATLAB** (R2020b or newer) or **GNU Octave** (for open‑source users)  
- Basic knowledge of MATLAB syntax and numerical methods  
- The toolbox magpie‑matlab is optionally included in `external/magpie‑matlab/` via a Git submodule. If you’ve not initialized it yet, run:

  ```bash
  git submodule update --init --recursive
  ```

  Then in MATLAB add its path:

  ```matlab
  addpath(genpath('external/magpie‑matlab'));
  savepath;
  ```

  This toolbox provides additional functionality for modal and matrix‐analysis tasks used in the course.

---

## Usage

1. Clone or download the repository:

   ```bash
   git clone https://github.com/mdphys/CISM_Udine2024.git
   cd CISM_Udine2024
   git submodule update --init --recursive
   ```

2. In MATLAB, set the working directory to `src`, e.g.:

   ```matlab
   cd('path/to/CISM_Udine2024/src');
   ```

3. Run an example script, for instance:

   ```matlab
   BarVariableThick;
   PoissonSolver;
   ```

   Then experiment by modifying parameters or the code to explore the numerical behaviour.

---

## Educational Purpose

This material supports lectures and workshops in the acoustics of musical instruments and, particularly, numerical analysis and eigenvalue problems.  
Grid‑based methods are employed to solve problems in one spatial dimension. Non‑uniform grids are introduced to address varying geometrical or material parameters, showing that finite‐difference methods are not restricted to uniform meshes.

---

## License

This project is distributed under the terms of the [MIT License](LICENSE).  
You are free to use, modify, and distribute the material. Attribution and including the original license are appreciated.

---

## Acknowledgments

Developed for the course *Physics of Musical Instruments Applied to Instrument Making*, Udine, Italy — May 2024.  
Special thanks to all teaching staff, workshop participants, and contributors.

---

## About

This repository holds the scripts used to produce the case studies and numerical examples in the lecture notes.

---
