# EadyModel_3D

Spectral and Finite Difference 3D Code for Solving the Eady Problem in Potential Vorticity Formulation

---

## Overview

EadyModel_3D is a Python-based repository that implements numerical techniques to solve the Eady problem, a fundamental model in geophysical fluid dynamics. This model describes baroclinic instability in the atmosphere and ocean, providing insight into the dynamics of mid-latitude weather systems.

The code utilizes spectral methods and finite-difference schemes to address the problem in a 3D potential vorticity framework, ensuring both precision and computational efficiency.

---

## Features

- Spectral Transforms: Efficient Fourier-based algorithms for spatial derivatives.
- Finite Difference Schemes: Methods to compute derivatives in the vertical direction.
- Dealiasing: Implementation of spectral dealiasing to avoid aliasing errors.
- Configurable Parameters: Adjustable grid resolution, domain size, and initial conditions.

---

## Code Structure

### eady.py
The main script that:
- Sets up the Eady problem.
- Defines the grid and physical parameters.
- Solves the equations (potential vorticity formulation) using spectral and finite-difference methods.
- 
---

### functions.py
A module with utility functions for Fourier transforms, derivative calculations, and wavenumber generation.

Key Functions:
- transform and inverse_transform: Perform forward and inverse Fourier transforms.
- X_derivative and Y_derivative: Compute spatial derivatives in spectral space.
- wavenumbers: Generate wavenumbers for Fourier transforms.
-

---
