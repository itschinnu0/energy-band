# Phase 2 — Semiconductor Physics Foundation

## Objective
Implement shared semiconductor calculations without any GUI dependency.

## Implement
- `constants.py`
- `units.py`
- `validation.py`
- `semiconductor.py`
- initial `results.py`
- domain-specific exceptions

## Required calculations
- thermal voltage
- temperature-dependent bandgap
- effective density of states
- intrinsic carrier concentration
- equilibrium Fermi level for n-type and p-type silicon
- Boltzmann quasi-Fermi carrier relationships
- energy/potential conversions

## Requirements
- Material parameters live in one place.
- Unit conversions are explicit.
- Use SI internally.
- Do not silently clip invalid inputs.
- Test representative and edge cases.

## Acceptance
- Unit tests pass.
- `EC-EV=Eg` conversion logic is testable.
- Known limiting/sign behavior is covered.
