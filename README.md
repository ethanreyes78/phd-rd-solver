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

The system solves the following reaction-diffusion equations, incorporating spatial heterogeneity and seasonality (Wang 2022):

$$ \frac{\partial S}{\partial t} - c_1 \frac{\partial^2 S}{\partial x^2} = \Lambda - \alpha(x,t) SB - \beta(x,t) SI - \mu S $$

$$ \frac{\partial I}{\partial t} - c_2 \frac{\partial^2 I}{\partial x^2} = \alpha(x,t) SB + \beta(x,t) SI - (\mu + \omega + \gamma)I $$

$$ \frac{\partial R}{\partial t} - c_3 \frac{\partial^2 R}{\partial x^2} = \gamma I - \mu R $$

$$ \frac{\partial B}{\partial t} + v\frac{\partial B}{\partial x} - d\frac{\partial^2 B}{\partial x^2} = gB\left(1 - \frac{B}{K}\right) + \xi(x,t) I - \delta B $$

for $$ 0 < x < L $$ and $$ t > 0 $$. 

The transmission rates fluctuate based on location and time of year using the multiplier:
$$ H(x)T(t) = (0.5 + 0.25\cos(2\pi x))(0.5 + 0.25\sin(2\pi t / 12)) $$

### Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Domain Length | `L` | $$3\pi$$ | Length of the theoretical river |
| Population influx rate | `Λ` | 1.0 | Recruitment rate into susceptible class |
| Env-to-human transmission | `α` | 0.0002 | Base rate of transmission from pathogen to human |
| Human-to-human transmission | `β` | 0.0001 | Base direct transmission rate between individuals |
| Natural death rate | `μ` | 0.0001 | Background mortality rate |
| Disease-induced death rate | `ω` | 0.01 | Additional mortality due to infection |
| Recovery rate | `γ` | 0.2 | Rate of recovery from infection |
| Pathogen removal rate | `δ` | 0.033 | Rate of pathogen removal from aquatic environment |
| Shedding rate | `ξ` | 0.3 | Base rate infected hosts shed pathogen |
| Pathogen growth rate | `g` | 0.05 | Intrinsic growth rate of waterborne pathogen |
| Carrying capacity | `K` | 100000.0 | Maximum pathogen concentration |
| Pathogen diffusion rate | `d` | 0.1 | Diffusion coefficient for pathogen in water |
| Advection rate | `v` | 0.5 | Speed of river flow |
| S diffusion rate | `c1` | 1.0 | Diffusion coefficient for susceptible individuals |
| I diffusion rate | `c2` | 0.5 | Diffusion coefficient for infected individuals |
| R diffusion rate | `c3` | 0.8 | Diffusion coefficient for recovered individuals |

### Boundary Conditions

- **Human Populations (S, I, R):** No-flux Neumann boundary conditions ($$ \frac{\partial U}{\partial x} = 0 $$), discretized using 2nd-order 3-point one-sided differences.
- **Pathogen (B):** Robin-type boundary conditions accounting for advection at the river boundaries. The downstream boundary incorporates a relative bacterial loss rate of $$ c = 2 $$ (Wang 2022).

### Initial Conditions

Populations are initialized with spatial clustering using periodic functions (Wang 2022):
- $$ S(x, 0) = 10000 - 500\cos(2x) $$
- $$ I(x, 0) = 1 - \cos(2x) $$
- $$ R(x, 0) = 0 $$
- $$ B(x, 0) = 0.5 - 0.3\cos(2x) $$

---

## Numerical Method

- **Method:** Explicit finite difference (Forward Time, Centered Space - FTCS).
- **Spatial discretization:** 2nd-order central difference for interior nodes; 2nd-order 3-point one-sided differences for boundaries.
- **Time discretization:** 1st-order Forward Euler.
- **Accuracy Verification:** Includes a Richardson Extrapolation module to empirically verify $$ O(\Delta x^2) $$ spatial accuracy.

---

## Files

| File | Description |
|---|---|
| `epi_rd_solver.py` | Main modular solver — includes parameters, FTCS scheme, 3D plotting, and accuracy checks. |
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
pip install numpy matplotlib
```

---

## Status

**Work in progress.** This solver is an introductory implementation developed to build familiarity with reaction-diffusion PDE systems before moving to more complex models.

---

## Author

Ethan Reyes  
PhD Student — University of Tennessee at Chattanooga  
Advisor: Dr. Jin Wang
