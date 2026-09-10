# Phase 10 — Streamlit GUI

## Objective
Expose the validated physics through a clean interactive interface.

## Implement
- device selector
- device-specific inputs
- voltage controls
- simulation/update/reset
- calculated result cards
- model name
- assumptions
- warnings/errors
- energy-band plot
- operating-region display
- theory/equation explanation panel

## Hard boundary
`app.py` collects inputs and renders results. It must not contain device equations.

## Acceptance
Changing an input calls the physics model and updates the plot/results.
Invalid inputs produce readable messages rather than crashes.
