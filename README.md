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

The system solves the following reaction-advection-diffusion equations, incorporating spatial heterogeneity and seasonality based on a bilinear incidence framework:

$$
\frac{\partial S}{\partial t} - c_1 \frac{\partial^2 S}{\partial x^2} = \Lambda - \alpha(x,t) SB - \beta(x,t) SI - \mu S
$$

$$
\frac{\partial I}{\partial t} - c_2 \frac{\partial^2 I}{\partial x^2} = \alpha(x,t) SB + \beta(x,t) SI - (\mu + \omega + \gamma)I
$$

$$
\frac{\partial R}{\partial t} - c_3 \frac{\partial^2 R}{\partial x^2} = \gamma I - \mu R
$$

$$
\frac{\partial B}{\partial t} + v\frac{\partial B}{\partial x} - d\frac{\partial^2 B}{\partial x^2} = g(t)B\left(1 - \frac{B}{K}\right) + \xi(x,t) I - \delta B
$$

for $0 < x < L$ and $t > 0$. 

### Spatio-Temporal Dynamics

To simulate realistic environmental conditions, the transmission and growth rates fluctuate based on location and time of year (seasonality):

- **Spatial Heterogeneity:** $H(x) = 0.5 + 0.25\cos(2x)$
- **Transmission Seasonality:** $T(t) = 0.5 + 0.25\sin(2\pi t / 12)$
- **Pathogen Growth Seasonality:** $g(t) = g_0 + g_1 \sin(\pi t / 6)$

The transmission rates ($\alpha, \beta, \xi$) are scaled by $H(x)T(t)$ to reflect both location-based risk (e.g., clustered villages) and seasonal weather changes. The pathogen growth rate ($g$) fluctuates on a 12-month cycle to simulate the effect of seasonal water temperature changes uniformly across the river.

### Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Domain Length | `L` | $3\pi$ | Length of the theoretical river |
| Population influx rate | `Λ` | 19.0 | Recruitment rate into susceptible class |
| Env-to-human transmission | `α` | 0.000033 | Base rate of transmission from pathogen to human |
| Human-to-human transmission | `β` | 0.00047 | Base direct transmission rate between individuals |
| Natural death rate | `μ` | 0.0019 | Background mortality rate (43.5 years) |
| Disease-induced death rate | `ω` | 0.001 | Additional mortality due to infection |
| Recovery rate | `γ` | 6.0 | Rate of recovery from infection (approx. 5 days) |
| Pathogen removal rate | `δ` | 1.0 | Rate of pathogen natural death in aquatic environment |
| Shedding rate | `ξ` | 300.0 | Base rate infected hosts shed pathogen |
| Base pathogen growth rate | `g_0` | 0.5 | Baseline intrinsic growth rate of waterborne pathogen |
| Seasonal growth amplitude | `g_1` | 0.05 | Amplitude of seasonal pathogen growth fluctuations |
| Carrying capacity | `K` | 2,000,000.0 | Maximum pathogen concentration in water |
| Pathogen diffusion rate | `d` | 0.1 | Diffusion coefficient for pathogen in water |
| Advection rate | `v` | 0.0 | Speed of river flow (Tested at 0.0, 0.1, 0.7, 1.2) |
| S diffusion rate | `c1` | 1.0 | Diffusion coefficient for susceptible individuals |
| I diffusion rate | `c2` | 0.5 | Diffusion coefficient for infected individuals |
| R diffusion rate | `c3` | 0.8 | Diffusion coefficient for recovered individuals |

### Boundary Conditions

- **Human Populations (S, I, R):** Zero-flux Neumann boundary conditions ($\frac{\partial U}{\partial x} = 0$), discretized using 2nd-order 3-point one-sided differences to preserve global $O(\Delta x^2)$ accuracy.
- **Pathogen (B):** Zero-flux boundary conditions accounting for advection at the river boundaries ($vB - dB_x = 0$). Also discretized using 2nd-order 3-point one-sided differences.

### Initial Conditions

Populations are initialized with spatial clustering to simulate distinct outbreak hotspots:
- $S(x, 0) = (\Lambda / \mu) - 500\cos(2x)$
- $I(x, 0) = 1 - \cos(2x)$
- $R(x, 0) = 0$
- $B(x, 0) = 0.5 - 0.3\cos(2x)$

---

## Numerical Method

- **Method:** Explicit finite difference (Forward Time, Centered Space - FTCS).
- **Spatial discretization:** 2nd-order central difference for interior nodes; 2nd-order 3-point one-sided differences for boundaries.
- **Time discretization:** 1st-order Forward Euler.
- **Accuracy Verification:** Includes a Richardson Extrapolation module to empirically verify $O(\Delta x^2)$ spatial accuracy using the $L_2$ norm.

---

## Files

| File | Description |
|---|---|
| `epi_rd_solver.py` | Main modular solver — includes parameters, vectorized FTCS scheme, 3D plotting, and accuracy checks. |
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



## Author

Ethan Reyes  
PhD Student — University of Tennessee at Chattanooga  
Advisor: Dr. Jin Wang
