# Finite Difference Option Pricing – Black–Scholes PDE (Implicit Euler)
Overview

This project implements a finite difference framework in C++ to solve the Black–Scholes partial differential equation (PDE) under the risk-neutral measure. The solver prices:

 1. European put options

 2. Down-and-out European barrier put options

The PDE is discretized using the Implicit Euler (Backward Time, Central Space – BTCS) scheme, and the resulting tridiagonal systems are solved efficiently via the Thomas algorithm.

The primary objective of this project was to construct a numerically stable, computationally efficient PDE-based pricing engine capable of handling both vanilla and barrier-style derivatives within a consistent risk-neutral framework.

# Mathematical Framework

Under the risk-neutral measure, the Black–Scholes PDE for a non-dividend-paying asset is:

∂V/∂t + (1/2) σ² S² ∂²V/∂S² + rS ∂V/∂S − rV = 0

Terminal condition (European Put):

V(S,T) = max(K − S, 0)

Boundary conditions (Vanilla Put):

V(0,t) = K e^(−r(T−t))

V(S_max,t) = 0

For the down-and-out barrier option with barrier level B:

V(S ≤ B, t) = 0 (absorbing boundary condition)

# Model Parameters

The numerical experiments were conducted using:

Initial stock price: S0 = 40

Strike price: K = 40

Time to maturity: T = 1 year

Risk-free rate: r = 10%

Volatility: σ = 30%

Maximum stock price boundary: S_max = 100

Spatial grid steps: M = 400

Time steps: N = 400

Barrier level (down-and-out case): B = 30

The spatial and temporal discretization parameters were chosen to ensure numerical stability and convergence.

# Numerical Methodology
Discretization Scheme

The PDE was discretized using:

Backward difference in time (first-order accurate)

Central difference in space (second-order accurate)

This yields the Implicit Euler (BTCS) scheme:

A V^n = V^(n+1)

At each time step, the discretized system forms a tridiagonal matrix, which is solved using the Thomas algorithm.

# Linear System Solver

The Thomas algorithm (specialized Gaussian elimination for tridiagonal matrices) was implemented with:

O(M) computational complexity per time step

O(M) memory usage

Forward elimination and backward substitution

Total computational complexity of the solver:

O(M × N)

This ensures scalability for finer grid refinements.

# European Put Pricing

The European put option was priced by:

Initializing the terminal payoff at maturity.

Iterating backward in time using the implicit scheme.

Enforcing boundary conditions at S = 0 and S = S_max.

Extracting the price at S0 from the spatial grid.

The numerical price converges to the analytical Black–Scholes solution as:

ΔS → 0

Δt → 0

Grid refinement confirms convergence and validates implementation correctness.

# Down-and-Out Barrier Extension

The solver was extended to handle a down-and-out barrier option by:

Imposing an absorbing boundary condition at S ≤ B

Forcing V(S ≤ B, t) = 0 at each time step

Ensuring the barrier does not destabilize the tridiagonal system

The barrier option price is strictly lower than the vanilla put due to knockout risk. Care was taken to maintain stability near the discontinuity introduced by the barrier.

# Stability and Convergence Analysis

The Implicit Euler scheme used in this project is:

Unconditionally stable (no CFL restriction)

First-order accurate in time

Second-order accurate in space

Stability arises from the diagonally dominant coefficient matrix. Numerical experiments show monotonic convergence toward the analytical Black–Scholes price for the European case.

Observed numerical effects:

Discretization bias for coarse grids

Numerical diffusion near the strike

Sensitivity to far-field truncation (choice of S_max)

Grid alignment sensitivity near barrier level

These effects diminish with grid refinement.
