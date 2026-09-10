# Phase 6 — NPN BJT Electrostatics and Classification

## Objective
Implement the analytical two-junction NPN BJT model.

## Lock
`VBE = VB - VE`
`VBC = VB - VC`

## Implement
- emitter/base/collector parameter validation
- EB built-in potential
- CB built-in potential
- EB/CB depletion widths
- complete EC/EV profile using one consistent electrostatic potential
- automatic operating-region classifier
- base depletion overlap detection
- region/terminal quasi-Fermi reference support only

## Operating regions
- EB reverse + CB reverse → cutoff
- EB forward + CB reverse → forward active
- EB forward + CB forward → saturation
- EB reverse + CB forward → reverse active

Near-zero bias must use a documented tolerance.

## Acceptance
No independent band shifting that violates `EC-EV=Eg`.
Overlap produces a controlled validity warning/error.
