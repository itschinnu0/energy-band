# Phase 5 — PN Bias and Quasi-Fermi References

## Objective
Extend PN to forward/reverse bias without falsely claiming equilibrium.

## Lock
`VD = Vbi - VA`.

Reject `VD <= 0` for the depletion approximation.

## Energy reference
Use:
- `EFn_n = 0 eV`
- `EFp_p = -VA eV`

Therefore:
`EFn_n - EFp_p = VA`.

## Implement
- biased depletion width/field/potential
- barrier variation
- operating condition
- neutral-region quasi-Fermi references
- voltage sweep support at the model layer

## Explicitly forbidden
- one flat global equilibrium EF for biased PN
- invented straight quasi-Fermi curves through depletion
- absolute-value-only quasi-Fermi separation

## Acceptance
- forward bias narrows depletion and lowers barrier
- reverse bias widens depletion and raises barrier
- zero-bias limit converges to flat EF
- signed quasi-Fermi separation matches VA
- invalid VD is rejected
