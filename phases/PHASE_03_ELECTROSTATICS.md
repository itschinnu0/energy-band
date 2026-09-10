# Phase 3 — Shared Electrostatics

## Objective
Build reusable 1-D electrostatic utilities.

## Implement
- piecewise charge-density representation
- Poisson-based electric-field construction
- potential integration
- boundary-condition handling
- consistency checks

## Required checks
- `dE/dx = rho/eps` in the modeled regions.
- `E = -dphi/dx` numerically.
- Potential difference matches the specified electrostatic drop.
- No invalid array shapes or non-finite values.

## Acceptance
The electrostatic layer can be used by PN and BJT without duplicating field/potential math.
