# EadyModel_3D

Spectral and Finite Difference 3D Code for Solving the Eady Problem in Potential Vorticity Formulation

---

## Overview

EadyModel_3D is a Python-based repository that implements numerical techniques to solve the Eady problem, a fundamental model in geophysical fluid dynamics.  
This model describes baroclinic instability in the atmosphere and ocean, providing insight into the dynamics of mid-latitude weather systems.

The code utilizes spectral methods and finite-difference schemes to address the problem in a 3D potential vorticity framework, ensuring both precision and computational efficiency.

---

## Features

- Spectral Transforms: Efficient Fourier-based algorithms for spatial derivatives.
- Finite Difference Schemes: Methods to compute derivatives in the vertical direction.
- Dealiasing: Implementation of spectral dealiasing to avoid aliasing errors.
- Configurable Parameters: Adjustable grid resolution, domain size, and initial conditions.
- 3D Compatibility: Solves the Eady problem in a 3D setup with variable vertical layers.

---

## Code Structure

### eady.py
The main script that:
- Sets up the Eady problem.
- Defines the grid and physical parameters.
- Solves the equations using spectral and finite-difference methods.

Key Functions:
- Initialization: Define grid dimensions and initial conditions.
- Numerical Solution: Solve potential vorticity equations iteratively.

---

### functions.py
A module with utility functions for Fourier transforms, derivative calculations, and wavenumber generation.

Key Functions:
- transform and inverse_transform: Perform forward and inverse Fourier transforms.
- X_derivative and Y_derivative: Compute spatial derivatives in spectral space.
- wavenumbers: Generate wavenumbers for Fourier transforms.

---

## Mathematical Formulation

The model is governed by the following key equations:

1. **Potential Vorticity Definition**:
   \[
   q = \nabla^2_h \psi + f + \frac{1}{\rho} \frac{\partial}{\partial z} 
   \left( \frac{\rho}{N^2} \frac{\partial \psi}{\partial z} \right)
   \]

2. **Conservation of Potential Vorticity**:
   \[
   \frac{D(q - f)}{Dt} = -v \beta
   \]

Where:
- \( q \): Potential vorticity.
- \( \psi \): Streamfunction.
- \( f \): Coriolis parameter.
- \( \rho \): Density.
- \( N^2 \): Brunt–Väisälä frequency.
- \( \beta \): Meridional gradient of the Coriolis parameter.
- \( v \): Meridional velocity.
- \( \nabla^2_h \): Horizontal Laplacian.
- \( \frac{D}{Dt} \): Material derivative.

These equations form the backbone of the Eady problem's dynamics and are solved using the numerical techniques implemented in this repository.

---
