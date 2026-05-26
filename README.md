# phd-rd-solver

A numerical solver for a reaction-diffusion PDE system modeling the spatial spread of waterborne infectious disease. This is a work in progress developed under PhD dissertation research.

---

## Overview

This project implements a finite difference method to solve a 1D reaction-diffusion system that models waterborne disease transmission — such as cholera — through both environment-to-human and human-to-human pathways. The spatial domain represents a river or aquatic environment containing a pathogen.

The model tracks four population/concentration variables over space and time:

| Variable | Description |
|---|---|
| `S(x, t)` | Susceptible individuals |
| `I(x, t)` | Infected individuals |
| `R(x, t)` | Recovered individuals |
| `B(x, t)` | Waterborne pathogen concentration in the aquatic environment |

---

## Model

The system solves the following reaction-diffusion equations:

$$\frac{\partial S}{\partial t} - c_1 \frac{\partial^2 S}{\partial x^2} = \Lambda - \alpha SB - \beta SI - \mu S$$

$$\frac{\partial I}{\partial t} - c_2 \frac{\partial^2 I}{\partial x^2} = \alpha SB + \beta SI - (\mu + \omega + \gamma)I$$

$$\frac{\partial R}{\partial t} - c_3 \frac{\partial^2 R}{\partial x^2} = \gamma I - \mu R$$

$$\frac{\partial B}{\partial t} + v\frac{\partial B}{\partial x} - d\frac{\partial^2 B}{\partial x^2} = gB\left(1 - \frac{B}{K}\right) + \xi I - \delta B$$

for $0 < x < a$ and $t > 0$.

### Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Population influx rate | `Λ` | 0.01 | Recruitment rate into susceptible class |
| Environment-to-human transmission | `α` | 0.2 | Rate of transmission from pathogen to human |
| Human-to-human transmission | `β` | 0.1 | Direct transmission rate between individuals |
| Natural death rate | `μ` | 5e-5 | Background mortality rate |
| Disease-induced death rate | `ω` | 0.01 | Additional mortality due to infection |
| Recovery rate | `γ` | 0.2 | Rate of recovery from infection |
| Pathogen removal rate | `δ` | 0.1 | Rate of pathogen removal from aquatic environment |
| Shedding rate | `ξ` | 0.3 | Rate infected hosts shed pathogen into environment |
| Pathogen growth rate | `g` | 0.3 | Intrinsic growth rate of waterborne pathogen |
| Carrying capacity | `K` | 2.0 | Maximum pathogen concentration |
| Pathogen diffusion rate | `d` | 1.0 | Diffusion coefficient for pathogen in water |
| Advection rate | `v` | 2.0 | Speed of river flow |
| S diffusion rate | `c1` | 1.0 | Diffusion coefficient for susceptible individuals |
| I diffusion rate | `c2` | 2.0 | Diffusion coefficient for infected individuals |
| R diffusion rate | `c3` | 1.5 | Diffusion coefficient for recovered individuals |

### Boundary Conditions

- No-flux boundary conditions for S, I, R
- Robin-type boundary conditions for B accounting for advection at the river boundaries

### Initial Conditions

- Susceptible population uniformly initialized at 100
- Single infected individual introduced at the midpoint of the spatial domain
- All other variables initialized at zero

---

## Numerical Method

- **Method:** Explicit finite difference (Forward Euler in time, 2nd order central difference in space)
- **Spatial discretization:** 2nd order central difference
- **Time discretization:** 1st order forward Euler
- **Domain:** 1D spatial domain $[0, a]$

---

## Files

| File | Description |
|---|---|
| `epi_rd_solver.py` | Main solver — defines parameters, reaction terms, and finite difference scheme |
| `environment.yml` | Conda environment for reproducibility |

---

## Requirements

Recreate the environment with:

```bash
conda env create -f environment.yml
conda activate rds1
```

Or install manually:

```bash
pip install numpy
```

---

## Status

🚧 **Work in progress.** This solver is an introductory implementation developed to build familiarity with reaction-diffusion PDE systems before moving to more complex models.

---

## Author

Ethan Reyes
PhD Student — University of Tennessee at Chattanooga
Advisor: Dr. Jin Wang
