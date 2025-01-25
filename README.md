# EadyModel_3D

Spectral and Finite Difference 3D Code for Solving the Eady Problem in Potential Vorticity Formulation.
Supervisor: Antonio Segalini

---

## Overview

EadyModel_3D is a Python-based repository that implements numerical techniques to solve the Eady problem, a fundamental model in geophysical fluid dynamics. This model describes baroclinic instability in the atmosphere and ocean, providing insight into the dynamics of mid-latitude weather systems.

The code utilizes spectral methods and finite-difference schemes to address the problem in a 3D potential vorticity framework, ensuring both precision and computational efficiency.

![Fun Animation](https://gifyu.com/image/SeZLj)

---

## Features

- **Spectral Transforms**: Efficient Fourier-based algorithms for spatial derivatives.
- **Finite Difference Schemes**: Methods to compute derivatives in the vertical direction.
- **Dealiasing**: Implementation of spectral dealiasing to avoid aliasing errors.
- **Configurable Parameters**: Adjustable grid resolution, domain size, and initial conditions.
- **Comprehensive Analysis**: Includes tools for growth rate computation, visualization, and velocity profile studies.

---

## Code Structure

### `eady.py`
The main module that:
- Sets up the Eady problem.
- Defines the grid and physical parameters.
- Solves the equations in the potential vorticity formulation using spectral and finite-difference methods.

---

### `functions.py`
A utility module with functions for Fourier transforms, derivative calculations, and wavenumber generation.

Key Functions:
- `transform` and `inverse_transform`: Perform forward and inverse Fourier transforms.
- `X_derivative` and `Y_derivative`: Compute spatial derivatives in spectral space.
- `wavenumbers`: Generate wavenumbers for Fourier transforms.

---

### `main.py`
The primary script for running simulations, performing analyses, and visualizing results. It integrates with `eady.py` and `functions.py` for exploring the Eady problem.

#### Key Functions:
- `eady_analysis`: Core function for running the Eady model and computing results.
- `plot_growth_rate_setups`: Visualizes growth rates across multiple setups.
- `analyze_and_plot_growth_rates`: Examines growth rates for various perturbations or random initial conditions.
- `plot_contour_of_v`: Generates contour plots of the velocity field.
- `analyze_growth_rate_vs_perturbation`: Studies the relationship between growth rates and perturbation modes.
---
