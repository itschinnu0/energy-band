# Energy Band Diagram Simulator

## Project
**Energy Band Diagram Simulator** — A physics-first semiconductor energy-band diagram calculator and interactive visualization tool.

## Purpose
Educational 1-D semiconductor energy-band visualization and electrostatics modeling. Provides academically defensible calculations and interactive visualization of band bending, potential, electric field, and quasi-Fermi levels under bias.

## Supported Devices
1. **1-D Abrupt PN Junction**: Abrupt depletion approximation, built-in barrier, forward/reverse bias, charge neutrality, and spatial band bending.
2. **NPN BJT**: Coupled two-junction analytical electrostatic model (N+ emitter, P base, N collector), automated operating-region classification (Cutoff, Forward active, Saturation, Reverse active), and base punch-through detection.
3. **NMOS / MOS Capacitor**: 1-D vertical MOS electrostatics on p-type silicon substrate, oxide capacitance, flat-band voltage, Level-1 analytical threshold model, Level-2 numerical Poisson-Boltzmann charge solving (`brentq`), and accumulation/depletion/inversion classification.

## Technology
- **Python 3.12**
- **NumPy**: Vectorized spatial profiles and array operations
- **SciPy**: Bounded numerical root finding (`brentq`)
- **Matplotlib**: Energy-band and voltage sweep plotting
- **Streamlit**: Web GUI and interactive controls
- **PyInstaller**: Standalone Windows executable packaging

## Supported Features
- Static and voltage-dependent energy-band diagrams ($E_C$, $E_V$, $E_i$, $E_F$, $E_{Fn}$, $E_{Fp}$)
- Region shading and interface boundaries
- Voltage-dependent band bending under bias
- PN depletion approximation and barrier variation
- BJT junction barrier lowering/raising and operating mode detection
- Level-1 analytical and Level-2 numerical MOS electrostatics
- Voltage sweeps (PN $V_A$, BJT $V_{BE}$, MOS $V_G$) with auxiliary plots
- PNG figure export
- CSV spatial profile and sweep data export
- Explicit error handling for unphysical parameter regimes without silent clipping

## Physics Limitations
- 1-D analytical and semi-analytical electrostatic models only.
- Not a full TCAD or process simulator.
- Does not calculate full drift-diffusion carrier transport equations or spatial quasi-Fermi profiles inside depletion regions.
- Does not model BJT transistor currents, recombination-generation, or gain ($\beta$, $\alpha$).
- Does not model lateral MOSFET transport ($I_D-V_{DS}$) or 2-D short-channel effects.
- Assumes ideal oxide with zero leakage and homogeneous silicon without quantum confinement.

## Development & Running

### Environment Setup
```powershell
uv sync
```

### Run Test Suite
```powershell
uv run pytest -v
```

### Launch Interactive GUI
```powershell
uv run streamlit run app.py
```

## Packaging & Windows Executable

Build a standalone single-file Windows executable:
```powershell
uv run pyinstaller --noconfirm --onefile --clean --copy-metadata streamlit --collect-data streamlit --collect-submodules streamlit --add-data "app.py;." launcher.py
```
The resulting executable is generated in `dist/launcher.exe`.
