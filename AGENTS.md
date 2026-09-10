# Energy Band Diagram Simulator — Coding Agent Instructions

## Mission

Build a physics-first Python application that interactively calculates and visualizes semiconductor energy-band diagrams for:

1. Abrupt 1-D PN junction
2. NPN BJT using a two-junction analytical electrostatic approximation
3. NMOS as a 1-D MOS electrostatic model

The application must be academically defensible. Correct physics, sign conventions, units, validation, and model limitations take priority over visual polish.

## Source of truth

`IMPLEMENTATION_PLAN_SOURCE.md` is the project specification supplied by the user.

The final authoritative corrections in Sections 36–50 take precedence over earlier sections if there is any conflict.

Do not silently replace the specified conventions with a different textbook convention. If a real inconsistency is discovered, stop and report it before implementing the affected model.

## Non-negotiable rules

### Physics
- Use `EC = EC_ref - q*phi`.
- Use `EV = EC - Eg`.
- Keep potential (V), energy (J), and plotted energy (eV) explicitly distinct.
- A common energy offset is allowed, but it must be deterministic and must not change energy differences.
- Homogeneous silicon must satisfy `EC - EV = Eg` within tolerance.
- Never fabricate a physical curve merely to make a plot look right.
- Never silently clip invalid physical states.

### PN junction
- Lock `VA = VP - VN`.
- `VA > 0` = forward bias; `VA < 0` = reverse bias.
- `VD = Vbi - VA`.
- The depletion approximation is valid only for `VD > 0`.
- At equilibrium, use one flat `EF`.
- Under bias, do not draw a single global equilibrium `EF`.
- Use neutral-region `EFp` and `EFn` references only.
- Enforce, in eV plotting units, `EFn_n - EFp_p = VA`.
- Do not invent a spatial quasi-Fermi profile through the depletion region.

### BJT
- Model NPN as N+ emitter / P base / N collector.
- Lock `VBE = VB - VE` and `VBC = VB - VC`.
- Classify cutoff, forward-active, saturation, and reverse-active from junction bias.
- Treat the two junctions as coupled electrostatic approximations.
- Detect base depletion overlap; do not silently generate an invalid diagram.
- Do not claim full spatial quasi-Fermi transport without transport equations.

### MOS/NMOS
- First model is vertical 1-D MOS electrostatics, not a complete lateral MOSFET.
- Use p-type substrate and x=0 at the Si/SiO2 interface, x>0 into bulk.
- `phi_bulk = 0`, `psi_s = phi(0)`.
- `psi_s > 0` means downward semiconductor band bending under the selected convention.
- `VFB = PhiMS - Qox/Cox`.
- Level-1 threshold:
  `VT = VFB + 2*phiF + sqrt(4*q*eps_si*NA*phiF)/Cox`.
- Do not use depletion-only charge as a complete accumulation/inversion model.
- Level-2 surface-potential solving must be bounded and must expose convergence diagnostics.

## Architecture

Keep the physics layer independent of Streamlit.

Expected structure:

```text
energy-band/
├── app.py
├── launcher.py
├── physics/
│   ├── __init__.py
│   ├── constants.py
│   ├── units.py
│   ├── validation.py
│   ├── results.py
│   ├── semiconductor.py
│   ├── electrostatics.py
│   ├── pn_junction.py
│   ├── bjt.py
│   ├── mos.py
│   └── solver.py
├── visualization/
│   ├── __init__.py
│   └── band_plot.py
├── tests/
│   ├── test_units.py
│   ├── test_semiconductor.py
│   ├── test_electrostatics.py
│   ├── test_pn.py
│   ├── test_bjt.py
│   ├── test_mos.py
│   ├── test_results.py
│   └── test_validation.py
├── pyproject.toml
├── uv.lock
├── README.md
└── Implementation_Plan.md
```

## Coding standards

- Python 3.12.
- Use type hints.
- Prefer dataclasses for parameter/result objects.
- Prefer pure functions for physics calculations.
- Avoid hidden global state.
- Use NumPy vectorization for spatial profiles.
- Use SciPy only where numerical solving is actually required.
- Use bounded/bracketed root finding for physical scalar roots.
- Raise domain-specific exceptions:
  - `InvalidParameterError`
  - `ModelValidityError`
  - `SolverConvergenceError`
  - `UnitConversionError`
- Do not put device equations inside `app.py`.
- Do not duplicate semiconductor constants across modules.
- Do not duplicate unit conversions.
- Avoid magic numbers; document material parameters and tolerances.
- Keep functions small enough to test independently.

## Result contract

Device simulations should return a common result object containing, as applicable:

- device
- model
- position
- EC
- EV
- Ei
- EF
- EFn
- EFp
- potential
- electric_field
- region
- calculated_parameters
- operating_condition
- warnings
- assumptions
- model_valid

Use `None` or masked values for unavailable physical quantities. Never use fabricated zero curves.

## Validation gate

A phase is not complete merely because the code runs.

Before moving to the next phase:
1. Run the relevant tests.
2. Check analytical identities.
3. Check limiting behavior.
4. Check units.
5. Check model-validity handling.
6. Report any unresolved discrepancy.

If a test fails because the specification appears physically inconsistent, do not weaken the test to make it pass. Explain the discrepancy.

## Development workflow

For each phase:
1. Read the corresponding phase task file.
2. Inspect existing code.
3. Implement the smallest coherent change.
4. Add/update tests at the same time.
5. Run the phase test command.
6. Run the full test suite if feasible.
7. Review for duplicated equations, unit mistakes, and sign mistakes.
8. Produce a concise completion report with:
   - files changed
   - physics implemented
   - tests run
   - test result
   - known limitations
   - next phase

## Explicit non-goals

Do not turn this into:
- full TCAD
- process simulation
- quantum device simulation
- full BJT transport simulator
- SPICE replacement
- cloud service
- database-backed application
- authentication system

The target is a focused, 1-D analytical energy-band visualization tool.
