# Phase 9 — Energy-Band Visualization

## Objective
Build Matplotlib visualization independent of Streamlit.

## Implement
- EC
- EV
- Ei where supported
- EF at equilibrium
- EFn/EFp only where model supports them
- region shading/labels
- depletion markers
- deterministic energy reference
- clear legends
- device-specific axis labels

## Integrity rules
- no fabricated zero curves for unavailable quantities
- no silicon bands inside oxide
- EC-EV=Eg
- quasi-Fermi lines only where supported
- biased PN must not visually imply global equilibrium

## Acceptance
Plot functions accept standardized simulation results and contain no device physics equations.
