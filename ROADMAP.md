# Build Roadmap

## Phase order

| Phase | Goal | Gate |
|---|---|---|
| 0 | Repository/agent setup | Clean project conventions |
| 1 | Environment + packaging spike | Minimal EXE starts |
| 2 | Semiconductor foundation | Physics/unit tests pass |
| 3 | Shared electrostatics | Field/potential identities pass |
| 4 | PN equilibrium | PN equilibrium tests pass |
| 5 | PN bias + quasi-Fermi references | Signed separation + zero-bias limit pass |
| 6 | BJT electrostatics/classification | Four operating regions + overlap detection pass |
| 7 | MOS Level-1 | VFB/VT/regions/bending tests pass |
| 8 | MOS Level-2 solver | Bounded solver diagnostics pass |
| 9 | Visualization | Plot integrity tests pass |
| 10 | Streamlit GUI | GUI calls physics only |
| 11 | Sweeps/comparison/export | Outputs preserve units/missing values |
| 12 | Full validation/polish | All tests pass |
| 13 | Final Windows packaging | Clean-machine EXE test passes |

## Critical rule

Do not build the polished GUI first.

The dependency chain is:

physics → validation → plotting → GUI → analysis → packaging

## Current starting point

Start with **Phase 1**. It is intentionally a deployment feasibility spike so PyInstaller/Streamlit problems are discovered before the project is large.
