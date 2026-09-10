# Phase 4 — PN Junction Equilibrium

## Objective
Implement the abrupt 1-D silicon PN equilibrium model.

## Lock
`VA = VP - VN`.

At equilibrium `VA=0`, one flat `EF` must be returned.

## Implement
- validated PN parameters
- built-in potential
- depletion width
- xp/xn
- electric field
- electrostatic potential
- EC/EV
- equilibrium EF
- region labels

## Mandatory identities
- `Vbi = VT ln(NA ND / ni^2)`
- `W = sqrt((2 eps/q)(1/NA+1/ND)Vbi)`
- `xp+xn=W`
- `NA*xp=ND*xn`
- `Emax` equality on both sides
- depletion potential drop = `Vbi`
- `EC-EV=Eg`
- EF is constant

## Acceptance
Reference PN case and limiting tests pass.
