# Phase 7 — NMOS/MOS Level-1 Electrostatics

## Objective
Implement vertical 1-D MOS electrostatics for p-type silicon.

## Lock
- x=0 at Si/SiO2 interface
- x>0 into silicon bulk
- bulk phi=0
- psi_s=phi(0)
- psi_s>0 → downward band bending

## Implement
- oxide capacitance
- phiF
- VFB
- VT
- Level-1 depletion approximation
- accumulation/depletion/inversion classification
- vertical EC/EV/Ei bending

## Equations
`Cox = eps_ox/tox`
`phiF = VT ln(NA/ni)`
`VFB = PhiMS - Qox/Cox`
`VT = VFB + 2 phiF + sqrt(4 q eps_si NA phiF)/Cox`

Before strong inversion:
`Wd = sqrt(2 eps_si psi_s / (q NA))`
`Qs = -q NA Wd`

## Acceptance
- tox increase lowers Cox
- VG progression produces accumulation → flat band → depletion → inversion
- positive psi_s gives downward band bending
- silicon bands are not drawn through oxide
