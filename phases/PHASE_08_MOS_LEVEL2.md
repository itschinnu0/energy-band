# Phase 8 — Numerical MOS Surface-Potential Solver

## Objective
Add a bounded numerical MOS charge/surface-potential solution after Level-1 is validated.

## Implement
- physically defined Qs(psi_s)
- bounded root function
- bracket search where necessary
- Brent-style bracketing solver
- residual calculation
- convergence diagnostics
- model-validity errors

## Requirements
The result must expose:
- converged
- iterations
- residual
- bounds
- warning/error

Do not silently clamp psi_s to force convergence.

## Acceptance
Numerical results agree with Level-1 in regimes where the approximation is valid, within documented tolerance.
