# Phase 11 — Voltage Sweep, Comparison, and Export

## Objective
Make voltage dependence the central interactive feature.

## Implement
- voltage slider
- voltage range/step
- multi-voltage comparison
- optional animation
- CSV export
- PNG export

## CSV columns
- position_m
- EC_eV
- EV_eV
- Ei_eV
- EF_eV
- EFn_eV
- EFp_eV
- potential_V
- electric_field_V_per_m
- region

Unavailable values must be blank/missing, not zero-filled.

Metadata should include device, model, temperature, applied voltage, assumptions.

## Acceptance
Sweep outputs preserve sign conventions and units.
