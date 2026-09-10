# Energy Band Diagram Simulator — Implementation Plan

## 1. Project Overview

**Project:** Python-based semiconductor energy-band diagram simulator with a graphical web interface.

**Primary objective:** Build a physically meaningful simulator that models and plots the energy-band diagrams of:

1. PN junction
2. NPN BJT
3. NMOS/MOS structure

The simulator will accept semiconductor/device parameters and applied voltages, calculate the corresponding electrostatic potential and energy-band positions, and interactively display the variation of:

\[
E_C,\quad E_V,\quad E_i,\quad E_F
\]

and, where required under non-equilibrium bias:

\[
E_{Fn},\quad E_{Fp}
\]

with position and applied voltage.

**Final deliverable:** A single Windows executable:

```text
Energy-Band-Simulator.exe
```

The application will provide an interactive interface where the user can select a device, modify its electrical parameters and applied voltage, and observe the corresponding energy-band diagram.

---

# 2. Model Reference and Locked Conventions

## 2.1 Primary reference

The implementation will use semiconductor-device physics references for the analytical models, with **S. M. Sze and Kwok K. Ng, _Physics of Semiconductor Devices_, 3rd ed.** as the primary device-physics reference, consistent with the MOS-capacitor project.

Where a numerical material parameter or approximation is required, the implementation will document its source.

If the course textbook uses a different convention, the course convention takes precedence and the implementation/report will explicitly state the difference.

---

# 2.2 Energy reference convention

The simulator will use electron energy in **electron-volts (eV)** for visualization.

For a semiconductor electrostatic potential \(\phi(x)\):

\[
E_C(x)=E_{C,ref}-q\phi(x)
\]

and:

\[
E_V(x)=E_C(x)-E_g
\]

For plotting in eV:

\[
E_C^{(eV)}(x)
=
E_{C,ref}^{(eV)}
-
\phi(x)
\]

when \(\phi\) is expressed in volts.

Similarly:

\[
E_V^{(eV)}(x)
=
E_C^{(eV)}(x)-E_g^{(eV)}
\]

The implementation must maintain a clear distinction between:

- electrical potential: V
- energy: J
- plotted energy: eV

The conversion must never be mixed implicitly.

An arbitrary constant may be added to all plotted energy levels because only energy differences are physically meaningful.

---

# 2.3 Fermi-level and quasi-Fermi-level convention

This section is **locked** and must be followed by all device models.

## Thermal equilibrium

At thermal equilibrium:

\[
E_F=\text{constant}
\]

throughout a connected semiconductor structure.

For an equilibrium PN junction:

```text
P region          depletion          N region

EC ─────────╲________________╱─────────
EF ─────────────────────────────────────
EV ─────────╲________________╱─────────
```

The bands bend while the Fermi level remains flat.

---

## Non-equilibrium bias

When an external bias is applied:

\[
V_A\neq0
\]

the PN junction is no longer in global thermal equilibrium.

Therefore:

\[
\boxed{\text{Do not represent a biased PN junction with one constant }E_F.}
\]

The implementation shall instead use quasi-Fermi levels:

- electron quasi-Fermi level \(E_{Fn}\)
- hole quasi-Fermi level \(E_{Fp}\)

These describe the non-equilibrium electron and hole populations.

---

## PN-junction voltage convention

For the PN-junction module:

\[
\boxed{V_A=V_P-V_N}
\]

where:

- \(V_P\) = externally applied P-side voltage
- \(V_N\) = externally applied N-side voltage

Therefore:

\[
V_A>0
\]

means forward bias, and:

\[
V_A<0
\]

means reverse bias.

The applied voltage is related to the quasi-Fermi-level separation by:

\[
\boxed{
E_{Fn,n}-E_{Fp,p}=qV_A
}
\]

where:

- \(E_{Fn,n}\) is the electron quasi-Fermi level in the quasi-neutral N region;
- \(E_{Fp,p}\) is the hole quasi-Fermi level in the quasi-neutral P region.

For \(V_A>0\):

\[
E_{Fn,n}>E_{Fp,p}
\]

and the separation is:

\[
E_{Fn,n}-E_{Fp,p}=qV_A
\]

When energies are plotted in eV, the numerical separation is:

\[
E_{Fn,n}^{(eV)}-E_{Fp,p}^{(eV)}=V_A
\]

---

## Spatial treatment of quasi-Fermi levels

The initial analytical PN model does **not** solve the full drift-diffusion transport equations.

Therefore:

- \(E_{Fp}\) will be treated as approximately constant in the quasi-neutral P region.
- \(E_{Fn}\) will be treated as approximately constant in the quasi-neutral N region.
- The depletion-region quasi-Fermi-level variation will **not** be artificially drawn as a straight line.
- If the depletion-region quasi-Fermi behavior is not explicitly solved, the plot will show the neutral-region quasi-Fermi levels and label the depletion region appropriately.

This prevents the simulator from presenting an unsupported quasi-Fermi profile as a physically exact result.

A future advanced transport model may calculate the spatial variation of \(E_{Fn}(x)\) and \(E_{Fp}(x)\) from carrier transport equations.

---

## Equilibrium limit

When:

\[
V_A=0
\]

the quasi-Fermi levels must converge to the equilibrium Fermi level:

\[
E_{Fn}=E_{Fp}=E_F
\]

in the neutral regions.

The simulator must therefore recover the standard flat-Fermi-level equilibrium diagram continuously as:

\[
V_A\rightarrow0
\]

---

# 3. Final Technology Stack

| Component | Technology | Purpose |
|---|---|---|
| Language | Python 3.12 | Core implementation |
| Package/environment | `uv` | Dependency and environment management |
| Numerical arrays | NumPy | Vectorized calculations |
| Scientific computing | SciPy | Numerical solving |
| Plotting | Matplotlib | Energy-band visualization |
| GUI/web interface | Streamlit | Interactive simulator UI |
| Testing | pytest | Automated validation |
| Packaging | PyInstaller | Single Windows executable |

## Core principle

The physics/model layer must remain independent of Streamlit.

The GUI is only a presentation and interaction layer.

---

# 4. Project Architecture

```text
energy-band/
│
├── app.py
│
├── physics/
│   ├── __init__.py
│   ├── constants.py
│   ├── parameters.py
│   ├── semiconductor.py
│   ├── pn_junction.py
│   ├── bjt.py
│   ├── mosfet.py
│   └── solver.py
│
├── visualization/
│   ├── __init__.py
│   └── band_plot.py
│
├── tests/
│   ├── test_semiconductor.py
│   ├── test_pn.py
│   ├── test_bjt.py
│   └── test_mosfet.py
│
├── launcher.py
├── pyproject.toml
├── uv.lock
├── README.md
└── Implementation.md
```

## Responsibility of each module

### `app.py`

- Streamlit application entry point.
- Collect user inputs.
- Select device.
- Call the physics model.
- Display calculated values.
- Display energy-band plots.
- Display operating-region information.
- Must not contain device equations.

### `physics/constants.py`

Contains:

- elementary charge \(q\)
- Boltzmann constant \(k\)
- vacuum permittivity \(\epsilon_0\)
- silicon relative permittivity
- silicon bandgap parameters
- intrinsic carrier concentration parameters
- documented semiconductor material parameters

### `physics/parameters.py`

Defines validated device parameters.

Examples:

- doping concentrations
- temperature
- applied voltage
- junction dimensions
- oxide thickness
- work-function difference
- material properties

### `physics/semiconductor.py`

Shared semiconductor calculations:

- thermal voltage
- intrinsic carrier concentration
- Fermi potential
- bandgap
- carrier concentrations
- equilibrium Fermi-level position
- energy/potential conversion

### `physics/pn_junction.py`

PN-junction model.

Responsibilities:

- calculate built-in potential
- calculate depletion width
- calculate electric field
- calculate electrostatic potential
- calculate \(E_C\)
- calculate \(E_V\)
- calculate equilibrium \(E_F\)
- calculate biased quasi-Fermi-level references
- determine bias condition
- generate energy-band data

### `physics/bjt.py`

NPN-BJT model.

Responsibilities:

- define emitter/base/collector structure
- calculate equilibrium junction barriers
- calculate emitter-base bias effects
- calculate collector-base bias effects
- calculate approximate band profile
- calculate appropriate quasi-Fermi-level references
- classify BJT operating condition

### `physics/mosfet.py`

NMOS/MOS structure model.

Responsibilities:

- calculate flat-band voltage
- calculate Fermi potential
- calculate threshold voltage
- calculate surface potential
- calculate semiconductor band bending
- classify accumulation/depletion/inversion
- generate semiconductor energy bands

### `physics/solver.py`

Numerical methods.

Responsibilities:

- root solving
- bounded numerical calculations
- surface-potential solving
- convergence/error handling
- numerical stability

### `visualization/band_plot.py`

Matplotlib-specific visualization.

Responsibilities:

- plot \(E_C\)
- plot \(E_V\)
- plot \(E_i\)
- plot equilibrium \(E_F\)
- plot \(E_{Fn}\)
- plot \(E_{Fp}\)
- mark material/device regions
- mark junction/depletion regions
- display voltage-dependent diagrams
- optional multiple-voltage comparison

### `tests/`

Automated validation of:

- physical constants
- analytical equations
- numerical calculations
- limiting behavior
- voltage dependence
- Fermi/quasi-Fermi conventions
- parameter sensitivity

### `launcher.py`

Packaging/runtime launcher.

Responsibilities:

- start Streamlit
- run locally
- open browser automatically
- support PyInstaller packaging

---

# 5. Physics Scope

The initial simulator will use **one-dimensional analytical semiconductor models**.

The project will contain three device modules.

```text
Device
│
├── PN Junction
│
├── NPN BJT
│
└── NMOS
```

The initial models will assume:

- silicon semiconductor
- one-dimensional geometry
- uniform doping within each region
- abrupt PN junctions
- non-degenerate semiconductor statistics
- Boltzmann approximation where applicable
- ideal interfaces initially
- temperature configurable
- no interface traps initially
- no tunneling
- no quantum confinement
- no advanced TCAD transport model

These assumptions must be documented in the application and final report.

---

# 6. Unit Convention

## GUI units

The GUI may use semiconductor-engineering units:

- doping: \(cm^{-3}\)
- length: nm or \(\mu m\)
- voltage: V
- temperature: K
- energy: eV

## Internal units

All physics calculations must use SI units where applicable.

Examples:

\[
N[cm^{-3}]
\rightarrow
N[m^{-3}]
\]

\[
L[\mu m]
\rightarrow
L[m]
\]

\[
q = 1.602\times10^{-19}\ C
\]

Energy will be converted consistently between joules and electron-volts.

---

# 7. Fundamental Semiconductor Equations

## 7.1 Thermal voltage

\[
V_T=\frac{kT}{q}
\]

This quantity is used throughout the semiconductor calculations.

---

# 7.2 Intrinsic carrier concentration

The temperature-dependent intrinsic carrier concentration will follow the same documented material-parameter strategy as the MOS-capacitance project.

The silicon bandgap is modeled using:

\[
E_g(T)
=
E_{g0}
-
\frac{\alpha T^2}{T+\beta}
\]

The effective density of states is:

\[
N_C(T)
=
N_C(300)
\left(\frac{T}{300}\right)^{3/2}
\]

\[
N_V(T)
=
N_V(300)
\left(\frac{T}{300}\right)^{3/2}
\]

and:

\[
n_i(T)
=
\sqrt{N_C(T)N_V(T)}
\exp
\left(
-\frac{E_g(T)}{2kT}
\right)
\]

The material parameters must be stored in `physics/constants.py` rather than scattered throughout the source code.

---

# 7.3 Fermi level in doped silicon

For non-degenerate p-type silicon:

\[
E_i-E_F
=
kT\ln\left(\frac{N_A}{n_i}\right)
\]

Therefore:

\[
E_F
=
E_i
-
kT\ln\left(\frac{N_A}{n_i}\right)
\]

For n-type silicon:

\[
E_F
=
E_i
+
kT\ln\left(\frac{N_D}{n_i}\right)
\]

The exact plotted energy reference may be shifted by an arbitrary constant because only energy differences are physically relevant.

---

# 7.4 Quasi-Fermi levels

Under non-equilibrium conditions:

\[
E_{Fn}\neq E_{Fp}
\]

in general.

Under the Boltzmann approximation:

\[
n=n_i
\exp
\left(
\frac{E_{Fn}-E_i}{kT}
\right)
\]

and:

\[
p=n_i
\exp
\left(
\frac{E_i-E_{Fp}}{kT}
\right)
\]

Therefore:

\[
np
=
n_i^2
\exp
\left(
\frac{E_{Fn}-E_{Fp}}{kT}
\right)
\]

This relationship will be used when interpreting non-equilibrium carrier populations.

---

# 7.5 Band-edge relationship

The conduction and valence bands satisfy:

\[
E_C-E_V=E_g
\]

Therefore:

\[
E_V(x)=E_C(x)-E_g
\]

The energy-band calculation will primarily derive the band bending from electrostatic potential.

---

# 8. PN Junction Model

## 8.1 Initial device

The initial PN-junction model will be:

```text
P-type silicon       N-type silicon

     NA                   ND
     │                    │
     │     depletion      │
─────┤<──────────────────>├─────
     │                    │
```

Assumptions:

- abrupt junction
- uniformly doped P and N regions
- one-dimensional structure
- depletion approximation
- silicon at configurable temperature
- equilibrium and externally biased conditions

---

# 8.2 Built-in potential

At equilibrium:

\[
V_{bi}
=
\frac{kT}{q}
\ln
\left(
\frac{N_A N_D}{n_i^2}
\right)
\]

This determines the equilibrium electrostatic barrier.

---

# 8.3 Applied voltage

Define:

\[
\boxed{V_A=V_P-V_N}
\]

with positive \(V_A\) representing forward bias.

The effective depletion potential is:

\[
V_D=V_{bi}-V_A
\]

For forward bias:

\[
V_A>0
\]

the barrier decreases.

For reverse bias:

\[
V_A<0
\]

the barrier increases.

The model must reject a depletion-approximation calculation when:

\[
V_D\le0
\]

rather than silently producing an invalid square root.

---

# 8.4 Depletion width

The total depletion width is:

\[
W
=
\sqrt{
\frac{2\epsilon_{si}}{q}
\left(
\frac{1}{N_A}+\frac{1}{N_D}
\right)
V_D
}
\]

The depletion widths on each side are:

\[
x_p=
\frac{N_D}{N_A+N_D}W
\]

\[
x_n=
\frac{N_A}{N_A+N_D}W
\]

Therefore:

\[
W=x_p+x_n
\]

The depletion region must extend farther into the more lightly doped side.

---

# 8.5 Electric field

Under the depletion approximation, the electric field is calculated from the charge density.

For the P side:

\[
\rho=-qN_A
\]

For the N side:

\[
\rho=+qN_D
\]

Using Poisson's equation:

\[
\frac{dE}{dx}
=
\frac{\rho}{\epsilon_{si}}
\]

The electric field is integrated piecewise across the depletion region.

---

# 8.6 Electrostatic potential

The electrostatic potential is obtained by:

\[
E(x)=-\frac{d\phi}{dx}
\]

The potential profile is integrated across the depletion region.

The resulting potential difference across the depletion region must satisfy:

\[
\Delta\phi=V_D
\]

within numerical tolerance.

This provides an important validation test.

---

# 8.7 Equilibrium energy bands

At:

\[
V_A=0
\]

the PN junction is in thermal equilibrium.

Therefore:

\[
\boxed{E_F=\text{constant}}
\]

The equilibrium band diagram is generated from:

\[
E_C(x)=E_{C,ref}-q\phi(x)
\]

and:

\[
E_V(x)=E_C(x)-E_g
\]

The simulator shall display one flat \(E_F\).

---

# 8.8 Biased energy bands

When:

\[
V_A\neq0
\]

the PN junction is non-equilibrium.

The electrostatic band profile is still obtained from the applied-bias depletion model:

\[
E_C(x)=E_{C,ref}-q\phi(x)
\]

\[
E_V(x)=E_C(x)-E_g
\]

However, the Fermi-level representation changes.

The simulator shall **not** display one flat equilibrium \(E_F\) across the biased junction.

Instead, it will display:

\[
E_{Fp,p}
\]

in the quasi-neutral P region and:

\[
E_{Fn,n}
\]

in the quasi-neutral N region.

The applied voltage must satisfy:

\[
\boxed{
E_{Fn,n}-E_{Fp,p}=qV_A
}
\]

or, in the eV plotting scale:

\[
\boxed{
E_{Fn,n}^{(eV)}-E_{Fp,p}^{(eV)}=V_A
}
\]

The absolute energy reference may be shifted arbitrarily.

---

# 8.9 Quasi-Fermi-level plotting rule

The initial PN model will not solve the complete semiconductor continuity and drift-diffusion equations.

Therefore the simulator shall:

- show \(E_{Fp}\) in the quasi-neutral P region;
- show \(E_{Fn}\) in the quasi-neutral N region;
- maintain the correct applied-bias separation;
- avoid inventing a spatial quasi-Fermi-level profile through the depletion region.

The depletion region may be visually marked as:

```text
Quasi-Fermi variation not explicitly solved
```

or the quasi-Fermi curves may be shown only in the regions where the approximation is defined.

This limitation must be documented.

---

# 8.10 Equilibrium limit of the biased model

As:

\[
V_A\rightarrow0
\]

the quasi-Fermi levels must satisfy:

\[
E_{Fn,n}-E_{Fp,p}\rightarrow0
\]

and therefore:

\[
E_{Fn,n}\rightarrow E_{Fp,p}\rightarrow E_F
\]

The simulator must pass this limit test.

---

# 8.11 PN-junction operating states

The GUI will classify:

```text
VA < 0
    → Reverse bias

VA = 0
    → Equilibrium

VA > 0
    → Forward bias
```

The simulator will show:

- barrier height
- depletion width
- electric field
- band bending
- applied voltage
- appropriate Fermi/quasi-Fermi representation

---

# 9. BJT Model

## 9.1 Initial device

The initial BJT model will be an **NPN transistor**:

```text
Emitter          Base          Collector

   N+              P              N
┌────────┐      ┌──────┐      ┌────────┐
│        │      │      │      │        │
│   N+   │──────│  P   │──────│   N    │
│        │      │      │      │        │
└────────┘      └──────┘      └────────┘
```

The model will treat the device as two coupled PN junctions:

```text
Emitter ── E-B junction ── Base ── B-C junction ── Collector
```

This is an analytical educational model rather than a full BJT transport simulator.

---

# 9.2 BJT inputs

The GUI will provide:

- emitter doping \(N_E\)
- base doping \(N_B\)
- collector doping \(N_C\)
- temperature
- \(V_{BE}\)
- \(V_{BC}\)

Optional parameters:

- emitter width
- base width
- collector width

---

# 9.3 Junction potentials

The emitter-base built-in potential is:

\[
V_{bi,EB}
=
\frac{kT}{q}
\ln
\left(
\frac{N_E N_B}{n_i^2}
\right)
\]

The collector-base built-in potential is:

\[
V_{bi,CB}
=
\frac{kT}{q}
\ln
\left(
\frac{N_C N_B}{n_i^2}
\right)
\]

Under applied junction voltages:

\[
V_{barrier,EB}
=
V_{bi,EB}-V_{BE}
\]

and:

\[
V_{barrier,CB}
=
V_{bi,CB}-V_{BC}
\]

The voltage polarity must be defined explicitly in the GUI.

---

# 9.4 BJT depletion widths

Each junction can initially use the abrupt-junction depletion approximation.

For the emitter-base junction:

\[
W_{EB}
=
\sqrt{
\frac{2\epsilon_{si}}{q}
\left(
\frac{1}{N_E}+\frac{1}{N_B}
\right)
V_{barrier,EB}
}
\]

For the collector-base junction:

\[
W_{CB}
=
\sqrt{
\frac{2\epsilon_{si}}{q}
\left(
\frac{1}{N_C}+\frac{1}{N_B}
\right)
V_{barrier,CB}
}
\]

The implementation must ensure:

\[
V_{barrier}>0
\]

before applying the depletion approximation.

---

# 9.5 BJT band construction

The BJT band diagram will be constructed by joining the energy profiles of:

```text
Emitter → E-B depletion → Base → B-C depletion → Collector
```

For each region:

1. establish doping-dependent equilibrium band reference;
2. calculate the local electrostatic potential;
3. apply the corresponding band shift;
4. maintain:

\[
E_C-E_V=E_g
\]

through homogeneous silicon.

The resulting plot will show:

\[
E_C(x)
\]

and:

\[
E_V(x)
\]

across the complete transistor.

---

# 9.6 BJT quasi-Fermi levels

Under equilibrium:

\[
E_{Fn}=E_{Fp}=E_F
\]

and the Fermi level is constant.

Under bias, the BJT is a non-equilibrium device.

Therefore the simulator shall not force one constant equilibrium Fermi level through the entire transistor.

The initial analytical model will use appropriate electron and hole quasi-Fermi-level references for the biased emitter, base, and collector regions.

The model will not claim to calculate the complete spatial variation of \(E_{Fn}(x)\) and \(E_{Fp}(x)\) unless the carrier transport equations are explicitly implemented.

---

# 9.7 BJT operating regions

The simulator will classify the BJT according to the bias of its two junctions.

| E-B junction | B-C junction | Operating region |
|---|---|---|
| Reverse | Reverse | Cutoff |
| Forward | Reverse | Forward active |
| Forward | Forward | Saturation |
| Reverse | Forward | Reverse active |

The GUI should display the detected region.

Example:

```text
VBE = 0.70 V
VBC < 0

Operating region:
FORWARD ACTIVE
```

The classification is based on junction bias, not merely on a manually selected label.

---

# 9.8 BJT limitation

The initial BJT model will focus on **energy-band electrostatics, junction barriers and operating-region visualization**.

It will not attempt to reproduce:

- complete Ebers-Moll transport
- recombination-generation distributions
- Early effect
- high-level injection
- avalanche breakdown
- detailed carrier transit-time behavior
- heterojunction effects

These are outside the initial scope.

---

# 10. MOSFET Model

## 10.1 Initial device

The initial MOSFET will be:

**NMOS with p-type silicon substrate.**

The basic structure is:

```text
              Gate
        ┌──────────────┐
        │    Metal     │
        └──────────────┘
             SiO₂
        ────────────────
          p-type Si
        ────────────────
        n+            n+
       Source        Drain
```

The energy-band calculation will initially focus on the vertical gate-to-semiconductor direction.

---

# 10.2 MOSFET inputs

The GUI will provide:

- substrate doping \(N_A\)
- oxide thickness \(t_{ox}\)
- gate area
- temperature
- gate voltage \(V_G\)
- metal-semiconductor work-function difference
- oxide charge

Initially:

\[
Q_{ox}=0
\]

---

# 10.3 Oxide capacitance

The oxide capacitance per unit area is:

\[
C'_{ox}
=
\frac{\epsilon_{ox}}{t_{ox}}
\]

where:

\[
\epsilon_{ox}=3.9\epsilon_0
\]

---

# 10.4 Fermi potential

For the p-type substrate:

\[
\phi_F
=
\frac{kT}{q}
\ln
\left(
\frac{N_A}{n_i}
\right)
\]

The sign convention must be explicitly maintained throughout the model.

---

# 10.5 Flat-band voltage

For the ideal model:

\[
V_{FB}=\Phi_{MS}
\]

For the extended model:

\[
V_{FB}
=
\Phi_{MS}
-
\frac{Q_{ox}}{C'_{ox}}
\]

The initial default is:

\[
Q_{ox}=0
\]

---

# 10.6 Threshold voltage

For the p-type substrate NMOS model:

\[
V_T
=
V_{FB}
+
2\phi_F
+
\frac{
\sqrt{4q\epsilon_{si}N_A\phi_F}
}{
C'_{ox}
}
\]

The exact convention will be verified against the selected reference before implementation.

---

# 10.7 Surface potential

The surface potential \(\psi_s\) describes semiconductor band bending.

The MOS voltage relationship is:

\[
V_G
=
V_{FB}
+
\psi_s
-
\frac{Q_s(\psi_s)}{C'_{ox}}
\]

The Level 1 implementation may use the depletion approximation.

The Level 2 implementation will solve the semiconductor charge equation numerically.

---

# 10.8 MOS energy bands

The surface conduction-band shift is related to surface potential by:

\[
\Delta E_C=-q\psi_s
\]

and:

\[
\Delta E_V=-q\psi_s
\]

The intrinsic level follows the same electrostatic shift.

Therefore:

\[
E_C(x)=E_{C,bulk}-q\phi(x)
\]

\[
E_V(x)=E_C(x)-E_g
\]

The band bending will be plotted from the bulk semiconductor toward the oxide interface.

---

# 10.9 MOS operating regions

For a p-type substrate:

```text
VG < VFB
    → Accumulation

VG ≈ VFB
    → Flat band

VFB < VG < VT
    → Depletion

VG > VT
    → Inversion
```

The simulator will automatically classify the operating state.

---

# 10.10 MOS inversion

When:

\[
V_G>V_T
\]

the surface enters strong inversion.

The energy bands will show stronger downward bending toward the semiconductor surface for the NMOS/p-type-substrate convention.

The simulator will explicitly display:

```text
Surface condition:
Strong inversion
```

and the calculated:

- \(V_{FB}\)
- \(V_T\)
- \(\psi_s\)
- depletion width
- band bending

---

# 11. Numerical Model Strategy

The simulator will use two conceptual levels where useful.

## Level 1 — Analytical model

Purpose:

- fast calculation
- easy debugging
- clear textbook comparison
- educational visualization

Level 1 will use:

- depletion approximation
- analytical junction electrostatics
- analytical MOS approximations

---

## Level 2 — Numerical model

The numerical model will be introduced where analytical approximations become insufficient.

Examples:

- MOS surface-potential solving
- smooth potential generation
- numerical integration
- bounded root solving

SciPy's bracketing methods such as `brentq` are preferred where a scalar physical root can be bracketed.

Unconstrained numerical iteration must not be used in production when it can wander into nonphysical regions.

---

# 12. Numerical Stability and Solver Bounds

The simulator must not allow invalid numerical operations.

Potential problems include:

- negative depletion potentials
- invalid square roots
- exponential overflow
- non-finite numerical values
- impossible voltage combinations

The production solver will:

1. establish physical bounds;
2. verify input parameters;
3. reject invalid physical states;
4. use bounded root solving;
5. check numerical residuals;
6. report controlled errors.

The implementation must not silently clip physically meaningful values merely to make a graph appear valid.

---

# 13. Voltage Sweep

A central feature of the project will be voltage-dependent visualization.

The user will be able to select:

```text
Voltage minimum
Voltage maximum
Voltage step
```

For example:

```text
-2.0 V
-1.5 V
-1.0 V
-0.5 V
 0.0 V
+0.5 V
+1.0 V
+1.5 V
+2.0 V
```

For every voltage:

```text
Applied voltage
      ↓
Device physics
      ↓
Electrostatic potential
      ↓
Energy bands
      ↓
Fermi/quasi-Fermi levels
      ↓
Band diagram
```

---

# 14. Interactive Voltage Control

The GUI should provide a voltage slider:

```text
Voltage

-2.0 V ─────────●──────── +2.0 V
                 ↑
              V = 0.4 V
```

Moving the slider must:

1. update the applied voltage;
2. recalculate the physical model;
3. regenerate the band data;
4. redraw the energy-band diagram;
5. update calculated parameters;
6. update operating-region classification;
7. update the Fermi/quasi-Fermi representation.

The diagram must therefore be generated dynamically rather than being a collection of manually prepared images.

---

# 15. Animation

After the core simulator is stable, an optional voltage-sweep animation may be added.

Example:

```text
V = -2 V
    ↓
V = -1 V
    ↓
V = 0 V
    ↓
V = +1 V
    ↓
V = +2 V
```

The animation should demonstrate the physical evolution of the bands.

Animation is an enhancement and must not be implemented before the underlying voltage-dependent model is validated.

---

# 16. GUI Requirements

## 16.1 Device selection

The UI should provide:

```text
Device:

[ PN Junction ▼ ]
```

Options:

```text
PN Junction
NPN BJT
NMOS
```

---

# 16.2 Device parameter controls

Parameters should change according to the selected device.

### PN Junction

- \(N_A\)
- \(N_D\)
- temperature
- applied voltage

### NPN BJT

- emitter doping
- base doping
- collector doping
- \(V_{BE}\)
- \(V_{BC}\)
- temperature

### NMOS

- \(N_A\)
- oxide thickness
- gate area
- temperature
- \(V_G\)
- \(\Phi_{MS}\)
- \(Q_{ox}\)

---

# 16.3 Simulation controls

Controls should include:

- Simulate / Update
- Reset
- Applied-voltage input
- Voltage slider
- Voltage range
- Voltage step
- Optional animation

---

# 16.4 Calculated results

The result panel should display device-specific quantities.

### PN Junction

- \(V_{bi}\)
- depletion width
- maximum electric field
- effective barrier
- bias condition
- \(E_{Fp,p}\)
- \(E_{Fn,n}\)
- quasi-Fermi separation

### BJT

- E-B barrier
- B-C barrier
- E-B depletion width
- B-C depletion width
- operating region
- appropriate quasi-Fermi-level references

### MOSFET

- \(C_{ox}\)
- \(\phi_F\)
- \(V_{FB}\)
- \(V_T\)
- \(\psi_s\)
- depletion width
- operating region

---

# 17. Energy-Band Plot Requirements

The main plot must show:

- X-axis: position
- Y-axis: energy (eV)
- \(E_C\)
- \(E_V\)
- \(E_i\) where appropriate
- \(E_F\) at equilibrium
- \(E_{Fn}\) where appropriate
- \(E_{Fp}\) where appropriate
- material/device region labels
- depletion region
- selected voltage

---

# 17.1 PN Junction plot

### Equilibrium

At:

\[
V_A=0
\]

the plot should show:

```text
Energy
  ↑
  │
  │ EC ───────╲____________╱──────
  │
  │ EF ───────────────────────────
  │
  │ EV ───────╲____________╱──────
  │
  └────────────────────────────────→ Position
          P       depletion       N
```

### Biased condition

At:

\[
V_A\neq0
\]

the plot must not show one flat \(E_F\).

Instead, it should show the appropriate quasi-Fermi levels:

```text
Energy
  ↑
  │ EC ───────╲____________╱──────
  │
  │ EFn ────────────────────────   N-side reference
  │
  │ EFp ────────────────────────   P-side reference
  │
  │ EV ───────╲____________╱──────
  │
  └────────────────────────────────→ Position
          P       depletion       N
```

The actual placement and region of the quasi-Fermi lines must follow the locked convention in Section 2.3.

The simulator must verify:

\[
E_{Fn,n}^{(eV)}
-
E_{Fp,p}^{(eV)}
=
V_A
\]

within numerical tolerance.

The depletion region shall not contain an invented linear quasi-Fermi profile unless a transport model has been implemented.

---

# 17.2 BJT plot

```text
Energy
  ↑
  │
  │ EC ────╲________╱────────╲____
  │
  │
  │ EFn / EFp where applicable
  │
  │ EV ────╲________╱────────╲____
  │
  └────────────────────────────────→ Position
       Emitter    Base      Collector
```

The plot must clearly identify:

```text
Emitter | Base | Collector
```

and:

```text
E-B junction
B-C junction
```

---

# 17.3 MOSFET plot

The vertical MOS energy-band diagram will show:

```text
Energy
  ↑

EC ───────────────╲
                   ╲
Ei ─────────────────╲
                     ╲
EV ────────────────╲
                    Surface
                      │
                      ↓
              Semiconductor depth
```

The plot should clearly indicate:

```text
Bulk semiconductor
        ↓
Surface
        ↓
Oxide interface
```

The exact orientation will be fixed in the visualization module and documented.

---

# 18. Multiple-Voltage Comparison

A comparison mode should allow the user to display multiple applied voltages simultaneously.

Example:

```text
PN Junction:

V = -1 V
V =  0 V
V = +0.5 V
```

This makes the voltage-dependent band movement immediately visible.

For PN junctions, the comparison must also show the corresponding Fermi/quasi-Fermi representation:

```text
V = -1 V
→ reverse-biased quasi-Fermi separation

V = 0 V
→ single equilibrium EF

V = +0.5 V
→ forward-biased quasi-Fermi separation
```

The plot must use a clear legend.

This feature is particularly important because the project requirement explicitly concerns **variation of the bands with applied voltage**.

---

# 19. Parameter Sweep Features

A parameter comparison mode should eventually support:

## PN junction doping

Compare:

\[
N_A,\ N_D
\]

Expected effects should be verified analytically.

---

## PN junction applied voltage

Compare:

\[
V_A=-1,\ 0,\ +0.5\ V
\]

Expected:

```text
Reverse bias
    ↓
larger barrier/depletion width

Equilibrium
    ↓
built-in barrier
flat EF

Forward bias
    ↓
smaller barrier/depletion width
non-equilibrium quasi-Fermi levels
```

---

## BJT base doping

Compare different base doping values and observe changes in:

- junction barriers
- depletion widths
- band profile

---

## MOS oxide thickness

Compare:

\[
t_{ox}=5,\ 10,\ 20\ nm
\]

and observe its effect on:

\[
C'_{ox}
\]

and voltage distribution across the oxide/semiconductor system.

---

## MOS gate voltage

Compare:

```text
Accumulation
Flat band
Depletion
Inversion
```

and display the corresponding band bending.

---

# 20. Export Features

Add:

- Export band data as CSV.
- Export current plot as PNG.

CSV structure should contain at minimum:

```text
Position,EC,EV,Ei,EF,EFn,EFp
...
```

Only applicable columns need to be populated for a particular device/mode.

For biased PN junctions, additional columns may include:

```text
Potential
ElectricField
Region
```

Device-specific calculated quantities may also be included.

---

# 21. Input Validation

The application must reject physically invalid values.

Examples:

- \(N_A\le0\)
- \(N_D\le0\)
- temperature \(\le0\)
- negative oxide thickness
- invalid voltage range
- invalid voltage step
- invalid junction barrier
- numerical solver failure

Errors should be displayed clearly in the GUI.

The application must not crash because of a user-entered invalid parameter.

---

# 22. Testing Strategy

Testing is a major part of making the simulator trustworthy.

## 22.1 Semiconductor tests

Test:

- thermal voltage
- intrinsic carrier concentration
- bandgap
- equilibrium Fermi-level position
- energy conversion
- quasi-Fermi carrier relationships

---

# 22.2 PN-junction tests

Test:

- built-in potential
- depletion width
- depletion-width partition
- electric-field maximum
- potential difference
- bandgap separation
- equilibrium flat Fermi level
- biased quasi-Fermi-level separation

Physical trends:

```text
Reverse bias ↑
→ depletion width ↑
→ barrier ↑

Forward bias ↑
→ depletion width ↓
→ barrier ↓
```

---

# 22.3 PN Fermi-level validation tests

These tests are mandatory.

### Equilibrium

For:

\[
V_A=0
\]

verify:

\[
E_{Fn,n}=E_{Fp,p}=E_F
\]

within numerical tolerance.

### Forward bias

For:

\[
V_A>0
\]

verify:

\[
E_{Fn,n}-E_{Fp,p}=qV_A
\]

or, in the eV plotting scale:

\[
E_{Fn,n}^{(eV)}-E_{Fp,p}^{(eV)}=V_A
\]

### Reverse bias

For:

\[
V_A<0
\]

the same signed relationship must hold:

\[
E_{Fn,n}-E_{Fp,p}=qV_A
\]

The implementation must not replace this signed relationship with:

\[
|E_{Fn}-E_{Fp}|=q|V_A|
\]

because doing so would lose the defined voltage polarity.

### No false equilibrium Fermi level

A biased PN junction must not return a single constant:

\[
E_F(x)=\text{constant}
\]

as its physical Fermi-level representation.

---

# 22.4 BJT tests

Test:

- E-B built-in potential
- B-C built-in potential
- depletion widths
- junction bias classification
- operating-region classification
- appropriate quasi-Fermi-level representation

Verify:

```text
EB forward + BC reverse
→ forward active
```

and:

```text
EB forward + BC forward
→ saturation
```

---

# 22.5 MOS tests

Test:

- oxide capacitance
- Fermi potential
- flat-band voltage
- threshold voltage
- surface potential
- accumulation/depletion/inversion classification

Verify:

```text
VG < VFB
→ accumulation

VFB < VG < VT
→ depletion

VG > VT
→ inversion
```

---

# 22.6 Energy-band tests

Every homogeneous silicon region must satisfy:

\[
E_C-E_V=E_g
\]

within numerical tolerance.

The calculated band bending must correspond to the calculated electrostatic potential:

\[
\Delta E_C=-q\Delta\phi
\]

This is one of the most important consistency tests in the project.

---

# 23. Limiting-Behavior Tests

Verify expected trends.

### PN junction

\[
|V_A|\uparrow\text{ in reverse bias}
\Rightarrow
W\uparrow
\]

### PN junction Fermi-level limit

\[
V_A\rightarrow0
\Rightarrow
E_{Fn}-E_{Fp}\rightarrow0
\]

and:

\[
E_{Fn},E_{Fp}\rightarrow E_F
\]

### Doping

Increasing one side's doping should cause a smaller fraction of the depletion region to extend into that side.

### MOS oxide

\[
t_{ox}\uparrow
\Rightarrow
C_{ox}\downarrow
\]

### MOS gate voltage

Increasing positive gate voltage on a p-type substrate should move the surface through:

```text
accumulation
→ flat band
→ depletion
→ inversion
```

The simulator must test these behaviors rather than relying on visual inspection.

---

# 24. Reference Parameter Sets

Use fixed reference configurations throughout development.

## PN Junction

```text
Material:       Silicon
NA:             1 × 10^16 cm^-3
ND:             1 × 10^16 cm^-3
Temperature:    300 K
Applied voltage: 0 V
```

Additional validation cases:

```text
VA = -1.0 V
VA = +0.2 V
VA = +0.5 V
```

For each case, validate:

- depletion width
- barrier
- band bending
- Fermi/quasi-Fermi representation

---

## NPN BJT

```text
Emitter doping:   1 × 10^18 cm^-3
Base doping:      1 × 10^17 cm^-3
Collector doping: 1 × 10^16 cm^-3
Temperature:      300 K
VBE:              0 V
VBC:              0 V
```

Forward-active test:

```text
VBE > 0
VBC < 0
```

---

## NMOS

```text
Substrate:      p-type silicon
NA:             1 × 10^16 cm^-3
tox:            10 nm
Temperature:    300 K
Phi_MS:         0 V
Qox:            0 C
```

Voltage sweep:

```text
-2 V to +3 V
```

These are development defaults, not universal physical constants.

---

# 25. Validation Philosophy

The simulator must not be validated by visual appearance alone.

Primary theoretical reference:

**S. M. Sze and Kwok K. Ng, _Physics of Semiconductor Devices_, 3rd ed.**

For each device:

1. Calculate key quantities independently.
2. Compare simulator results against analytical equations.
3. Verify voltage dependence.
4. Verify energy/potential consistency.
5. Verify bandgap separation.
6. Verify physical operating-region classification.
7. Verify Fermi/quasi-Fermi conventions.
8. Compare the overall shape with trusted textbook/reference diagrams.

The report should include validation tables.

Example:

| Quantity | Reference | Simulator | Relative Error |
|---|---:|---:|---:|
| \(V_{bi}\) | ... | ... | ... |
| \(W\) | ... | ... | ... |
| \(V_T\) | ... | ... | ... |
| \(\psi_s\) | ... | ... | ... |
| \(E_{Fn}-E_{Fp}\) | ... | ... | ... |

---

# 26. Development Phases

## Phase 1 — Environment and Packaging Feasibility Spike

- Initialize project with `uv`.
- Pin Python 3.12.
- Add runtime dependencies.
- Add development dependencies.
- Confirm imports.
- Create initial project structure.
- Create minimal Streamlit application.
- Package the minimal application using PyInstaller as a single executable.
- Verify that the frozen executable starts Streamlit.
- Verify browser auto-launch.
- Test outside the active development environment.
- Record required PyInstaller hidden imports and assets.

Do not wait until the final stage to discover that the Streamlit + PyInstaller deployment strategy has problems.

---

# Phase 2 — Semiconductor Physics Foundation

Implement:

- physical constants
- material parameters
- unit conversion
- thermal voltage
- bandgap
- intrinsic carrier concentration
- equilibrium Fermi-level calculations
- quasi-Fermi-level relationships
- energy/potential conversion

No Streamlit physics logic.

---

# Phase 3 — PN Junction Equilibrium

Implement:

- built-in potential
- depletion width
- depletion charge
- electric field
- electrostatic potential
- \(E_C\)
- \(E_V\)
- equilibrium \(E_F\)

Generate the equilibrium band diagram.

Validate:

\[
E_F=\text{constant}
\]

---

# Phase 4 — PN Junction Bias

Implement:

- forward bias
- reverse bias
- voltage sweep
- band bending changes
- depletion-width changes
- barrier changes
- biased quasi-Fermi-level references

Implement the locked relationship:

\[
E_{Fn,n}-E_{Fp,p}=qV_A
\]

Validate:

- forward bias sign
- reverse bias sign
- zero-bias limit
- no false flat equilibrium \(E_F\)

---

# Phase 5 — BJT Energy-Band Model

Implement:

- NPN structure
- emitter/base/collector doping
- E-B junction
- B-C junction
- depletion widths
- band profiles
- junction bias
- quasi-Fermi-level representation where supported
- operating-region classification

Validate:

```text
Cutoff
Forward active
Saturation
Reverse active
```

---

# Phase 6 — MOSFET Energy-Band Model

Implement:

- p-type substrate
- oxide
- gate
- \(V_{FB}\)
- \(V_T\)
- surface potential
- depletion width
- accumulation
- depletion
- inversion
- energy-band bending

Reuse validated semiconductor and electrostatic utilities rather than duplicating calculations.

---

# Phase 7 — Numerical Surface-Potential Solver

Implement:

- \(Q_s(\psi_s)\)
- bounded surface-potential solving
- analytical derivative where required
- numerical stability
- convergence checking
- comparison against analytical approximation

The numerical model must not be introduced before the analytical model is validated.

---

# Phase 8 — Physics Validation

Run:

- unit tests
- reference parameter tests
- limiting-behavior tests
- energy/potential consistency tests
- bandgap consistency tests
- PN Fermi/quasi-Fermi tests
- voltage-sweep validation
- BJT operating-region validation
- MOS operating-region validation

Fix sign and unit errors before GUI development.

---

# Phase 9 — Streamlit GUI

Implement:

- device selection
- parameter controls
- voltage slider
- simulation controls
- result cards
- energy-band plot
- Fermi/quasi-Fermi display
- operating-region display
- theory/equation information

The GUI must call the physics layer rather than duplicate it.

---

# Phase 10 — Voltage Sweep and Animation

Add:

- voltage sweep
- multiple-voltage comparison
- interactive slider
- optional animation
- selected-voltage marker
- Fermi/quasi-Fermi update during voltage changes

This phase directly addresses the requirement of showing band variation with applied voltage.

---

# Phase 11 — Analysis Features

Add:

- parameter sweeps
- doping comparison
- oxide-thickness comparison
- voltage comparison
- multiple curves
- CSV export
- PNG export

---

# Phase 12 — Polish

Add:

- theory/equation panel
- unit explanations
- assumptions
- validation information
- clear operating-region descriptions
- Fermi/quasi-Fermi explanation
- professional layout
- example configurations

---

# Phase 13 — Testing

Run:

- full pytest suite
- edge cases
- invalid-input tests
- numerical solver tests
- Fermi-level convention tests
- export tests
- voltage-sweep tests
- device switching tests

---

# Phase 14 — Packaging

- Add launcher.
- Configure PyInstaller.
- Resolve hidden imports.
- Resolve Streamlit assets.
- Build single executable.
- Test on clean Windows environment.
- Verify browser auto-launch.
- Verify simulator works without Python installed.

---

# 27. Final Application Feature Set

## Required

- PN-junction energy-band model.
- NPN-BJT energy-band model.
- NMOS energy-band model.
- Configurable doping.
- Configurable temperature.
- Applied-voltage control.
- \(E_C\) plot.
- \(E_V\) plot.
- \(E_F\) at equilibrium.
- \(E_{Fn}\) and \(E_{Fp}\) where appropriate.
- Correct biased PN-junction quasi-Fermi-level representation.
- Band bending.
- Voltage sweep.
- Operating-region identification.
- Streamlit GUI.
- Input validation.
- Automated tests.
- Single Windows executable.

## Strong enhancements

- Multiple-voltage comparison.
- Voltage animation.
- Parameter sweeps.
- CSV export.
- PNG export.
- Theory/equation panel.
- Calculated electric field.
- Calculated electrostatic potential.
- Quasi-Fermi-level analysis.

## Optional, only after the core is stable

- advanced non-ideal effects
- temperature sweeps
- material selection
- additional semiconductor materials
- advanced carrier transport
- full spatial quasi-Fermi-level calculation

---

# 28. Non-Goals

Do not expand the initial project into:

- full TCAD simulation
- process simulation
- fabrication simulation
- quantum mechanical device simulation
- complete BJT transport simulation
- SPICE replacement
- semiconductor manufacturing simulation
- cloud backend
- database
- authentication
- multi-user infrastructure

The project should remain focused on:

> **Physics-based energy-band visualization and its variation under applied voltage.**

---

# 29. Packaging Plan

The final target is:

```text
Energy-Band-Simulator.exe
```

The executable will package:

- Python runtime
- NumPy
- SciPy
- Matplotlib
- Streamlit
- application code
- required runtime assets

Expected runtime:

```text
Double-click EXE
       ↓
Launcher starts
       ↓
Local Streamlit server starts
       ↓
Browser opens automatically
       ↓
Energy Band Simulator
```

The executable must be tested independently from the development environment.

---

# 30. Documentation and Academic Deliverables

The project should eventually include:

1. Project README.
2. Implementation documentation.
3. Semiconductor theory reference.
4. Equation reference.
5. GUI screenshots.
6. Validation results.
7. Voltage-sweep results.
8. Parameter-sweep results.
9. Example energy-band diagrams.
10. Fermi/quasi-Fermi-level explanation.
11. Executable.
12. Source code.
13. Dependency lock file.

The final report must explain not only **what** was implemented, but **why the applied voltage changes the energy bands and Fermi/quasi-Fermi levels in the observed way**.

---

# 31. Quality Criteria

The project is considered complete only when:

- The semiconductor equations are implemented correctly.
- Sign conventions are explicitly documented.
- Units are handled consistently.
- Energy and electrostatic potential are correctly related.
- PN-junction band bending behaves physically.
- Equilibrium PN junction shows a single flat \(E_F\).
- Biased PN junction does not incorrectly show a single equilibrium \(E_F\).
- Biased PN junction uses the appropriate quasi-Fermi-level representation.
- The signed quasi-Fermi separation agrees with the applied voltage.
- BJT junction bias changes the band profile correctly.
- MOS gate voltage produces the expected accumulation/depletion/inversion behavior.
- Bandgap remains physically consistent.
- Numerical results are independently validated.
- The GUI rejects invalid inputs.
- Automated tests pass.
- Voltage sweeps work correctly.
- The single executable launches successfully.
- The project can be explained clearly during a viva.

---

# 32. Implementation Principle

The order is intentionally:

```text
Semiconductor physics
        ↓
Analytical device models
        ↓
Numerical model
        ↓
Fermi/quasi-Fermi validation
        ↓
Energy-band visualization
        ↓
Interactive GUI
        ↓
Voltage sweep / animation
        ↓
Analysis features
        ↓
Testing
        ↓
Packaging
```

Not:

```text
Pretty graph
    ↓
Manually move bands
    ↓
Random equations
    ↓
Graph that looks correct
```

The credibility of the simulator comes from the physics model first.

---

# 33. Central Mathematical Principle

All three device models ultimately connect **electrostatic potential** to **energy bands**.

The central relationship is:

\[
\boxed{
E_C(x)=E_{C,ref}-q\phi(x)
}
\]

\[
\boxed{
E_V(x)=E_C(x)-E_g
}
\]

For non-equilibrium semiconductor regions:

\[
\boxed{
np=n_i^2
\exp
\left(
\frac{E_{Fn}-E_{Fp}}{kT}
\right)
}
\]

For the biased PN junction under the locked voltage convention:

\[
\boxed{
V_A=V_P-V_N
}
\]

and:

\[
\boxed{
E_{Fn,n}-E_{Fp,p}=qV_A
}
\]

Therefore:

```text
Applied voltage
      ↓
Charge distribution
      ↓
Electric field
      ↓
Electrostatic potential
      ↓
Energy-band shift
      ↓
EC / EV / Ei
      ↓
Fermi or quasi-Fermi levels
      ↓
Energy-band diagram
```

This relationship is the core of the entire project.

The PN junction, BJT, and MOSFET modules differ primarily in **how their charge distribution and electrostatic potential are established**.

---

# 34. Recommended Final GUI Structure

```text
┌──────────────────────────────────────────────────┐
│          ENERGY BAND DIAGRAM SIMULATOR            │
├──────────────────────────────────────────────────┤
│ Device                                           │
│ [ PN Junction ▼ ]                                │
│                                                  │
│ Parameters                                       │
│                                                  │
│ NA:              1 × 10¹⁶ cm⁻³                 │
│ ND:              1 × 10¹⁶ cm⁻³                 │
│ Temperature:     300 K                          │
│                                                  │
│ Applied Voltage                                 │
│ -2 V ─────────────●──────────── +2 V            │
│                    ↑                             │
│                  0.4 V                           │
├──────────────────────────────────────────────────┤
│                                                  │
│              ENERGY BAND DIAGRAM                 │
│                                                  │
│ Energy (eV)                                      │
│   ↑                                              │
│   │ EC ───────╲____________╱────                │
│   │                                              │
│   │ EFn ───────────────────────                 │
│   │ EFp ───────────────────────                 │
│   │                                              │
│   │ EV ───────╲____________╱────                │
│   │                                              │
│   └──────────────────────────────→ Position      │
│        P      Depletion       N                  │
├──────────────────────────────────────────────────┤
│ Operating condition: Forward Bias                │
│ Built-in potential:    ... V                    │
│ Depletion width:       ... μm                   │
│ Barrier:               ... eV                   │
│ Quasi-Fermi separation: ... eV                  │
└──────────────────────────────────────────────────┘
```

The controls and calculated results must change automatically when the selected device changes.

For equilibrium PN junctions, the GUI shall replace the quasi-Fermi-level display with the single equilibrium \(E_F\).

For biased PN junctions, the GUI shall clearly label \(E_{Fn}\) and \(E_{Fp}\) and must not label either one as the global equilibrium Fermi level.

---

# 35. Final Project Definition

The completed project is:

> **A Python-based interactive semiconductor energy-band simulator that calculates and visualizes the variation of conduction-band, valence-band, intrinsic-level and appropriate Fermi/quasi-Fermi levels for PN junctions, NPN BJTs and NMOS structures as a function of applied electrical bias.**

The project is fundamentally a **physics simulator**, with Streamlit serving only as the interface.

The implementation must therefore prioritize:

\[
\boxed{
\text{Correct physics}
>
\text{Numerical correctness}
>
\text{Validation}
>
\text{Visualization}
>
\text{GUI polish}
}
\]

This ordering is mandatory for the project to be academically defensible.
---

# 36. FINAL PHYSICS AND IMPLEMENTATION CORRECTIONS

This section is **authoritative**. If any earlier section conflicts with this section, this section takes precedence. The purpose is to remove ambiguities before programming begins.

## 36.1 Energy-band convention

Use the electrostatic-potential relationship everywhere:

\[
E_C(x)=E_{C,ref}-q\phi(x)
\]

\[
E_V(x)=E_C(x)-E_g
\]

For plotting in eV:

\[
E_C^{(eV)}(x)=E_{C,ref}^{(eV)}-\phi(x)
\]

when \(\phi\) is in volts.

The code shall distinguish explicitly between:

- potential: V
- energy: J
- plotted energy: eV

A common vertical energy offset is permitted because only energy differences are physical. The offset must never change any energy difference or voltage-dependent result.

Every homogeneous silicon region must satisfy:

\[
E_C-E_V=E_g
\]

within numerical tolerance.

---

# 37. GLOBAL ENERGY-REFERENCE RULE

Every generated band diagram must have a deterministic energy reference.

For the initial PN analytical model, use the quasi-neutral N-side electron quasi-Fermi level as the plotting reference:

\[
\boxed{E_{Fn,n}=0\ \text{eV}}
\]

For the locked PN voltage convention:

\[
V_A=V_P-V_N
\]

use:

\[
\boxed{E_{Fp,p}=-V_A\ \text{eV}}
\]

so that:

\[
E_{Fn,n}^{(eV)}-E_{Fp,p}^{(eV)}=V_A
\]

At equilibrium:

\[
E_{Fn,n}=E_{Fp,p}=E_F=0\ \text{eV}
\]

A different common offset may be used internally, but it must be applied to all energy curves equally.

---

# 38. PN JUNCTION — FINAL LOCKED THEORY

## 38.1 Applied-voltage convention

Lock:

\[
\boxed{V_A=V_P-V_N}
\]

where \(V_P\) and \(V_N\) are the externally applied P- and N-side terminal potentials.

Therefore:

- \(V_A>0\): forward bias
- \(V_A=0\): equilibrium
- \(V_A<0\): reverse bias

This convention must be used consistently in equations, GUI labels, plots, tests, and documentation.

## 38.2 Built-in and depletion potential

\[
V_{bi}=V_T\ln\left(\frac{N_AN_D}{n_i^2}\right)
\]

\[
\boxed{V_D=V_{bi}-V_A}
\]

The abrupt depletion model is valid only while:

\[
V_D>0
\]

If \(V_D\le0\), the program must return a controlled model-validity error. It must not take the square root of a negative value and must not silently clip \(V_D\) to zero.

## 38.3 Depletion width

\[
W=
\sqrt{\frac{2\epsilon_{si}}{q}
\left(\frac{1}{N_A}+\frac{1}{N_D}\right)V_D}
\]

\[
x_p=\frac{N_D}{N_A+N_D}W
\]

\[
x_n=\frac{N_A}{N_A+N_D}W
\]

and:

\[
W=x_p+x_n
\]

The depletion region must extend farther into the more lightly doped side.

## 38.4 Electric field and potential

Use Poisson's equation:

\[
\frac{dE}{dx}=\frac{\rho}{\epsilon_{si}}
\]

with depletion charge:

\[
\rho=-qN_A \quad \text{on the P side}
\]

\[
\rho=+qN_D \quad \text{on the N side}
\]

and:

\[
E=-\frac{d\phi}{dx}
\]

The code must derive the electric field and potential from the selected coordinate convention rather than manually shaping curves.

Validate:

\[
E_{max}=\frac{qN_Ax_p}{\epsilon_{si}}
=\frac{qN_Dx_n}{\epsilon_{si}}
\]

and the potential drop across depletion must equal \(V_D\), with the sign matching the chosen coordinate/potential convention.

## 38.5 Equilibrium PN band diagram

At \(V_A=0\):

\[
E_F=\text{constant}
\]

The simulator shall show one flat equilibrium Fermi level.

The bands are obtained from the electrostatic potential:

\[
E_C(x)=E_{C,ref}-q\phi(x)
\]

\[
E_V(x)=E_C(x)-E_g
\]

## 38.6 Biased PN band diagram

For \(V_A\ne0\), the junction is not in global thermal equilibrium. Do **not** plot one constant equilibrium \(E_F\) across the entire junction.

The initial model shall display:

- \(E_{Fp}\) only in the quasi-neutral P region;
- \(E_{Fn}\) only in the quasi-neutral N region;
- no invented straight quasi-Fermi curve through the depletion region.

Use:

\[
\boxed{E_{Fn,n}^{(eV)}-E_{Fp,p}^{(eV)}=V_A}
\]

as the locked quasi-neutral reference relationship.

The model does not solve the drift-diffusion/continuity equations, so it must not claim to calculate the complete spatial \(E_{Fn}(x)\) and \(E_{Fp}(x)\) profiles.

## 38.7 Terminal voltage versus depletion voltage

Do not equate the terminal voltage directly with the depletion-region potential drop.

The depletion electrostatics use:

\[
V_D=V_{bi}-V_A
\]

to determine depletion charge, field, width, and band bending.

The terminal voltage is represented separately through the terminal/quasi-Fermi reference convention.

---

# 39. BJT — FINAL LOCKED THEORY

## 39.1 Scope

The first BJT implementation is an **analytical electrostatic energy-band model**, not a full BJT transport simulator.

Structure:

```text
N+ emitter | P base | N collector
```

Treat the device initially as two coupled abrupt PN junctions for electrostatic calculations.

## 39.2 Junction voltage conventions

Lock:

\[
\boxed{V_{BE}=V_B-V_E}
\]

\[
\boxed{V_{BC}=V_B-V_C}
\]

Positive \(V_{BE}\) is the forward-bias direction of the E-B junction.

Positive \(V_{BC}\) is the forward-bias direction of the B-C junction.

Then:

\[
V_{barrier,EB}=V_{bi,EB}-V_{BE}
\]

\[
V_{barrier,CB}=V_{bi,CB}-V_{BC}
\]

where:

\[
V_{bi,EB}=V_T\ln\left(\frac{N_EN_B}{n_i^2}\right)
\]

\[
V_{bi,CB}=V_T\ln\left(\frac{N_CN_B}{n_i^2}\right)
\]

The depletion approximation is valid only when the corresponding effective barrier is positive.

## 39.3 BJT operating-region classifier

| E-B junction | B-C junction | Operating region |
|---|---|---|
| Reverse | Reverse | Cutoff |
| Forward | Reverse | Forward active |
| Forward | Forward | Saturation |
| Reverse | Forward | Reverse active |

Near-zero junction voltages shall use a configurable tolerance and may be reported as `Transition / near-zero bias`.

## 39.4 BJT band construction

Construct:

```text
Emitter → E-B depletion → Base → B-C depletion → Collector
```

using one consistent electrostatic potential.

Never independently move \(E_C\) and \(E_V\). Calculate:

\[
E_C=E_{C,ref}-q\phi
\]

then:

\[
E_V=E_C-E_g
\]

This guarantees bandgap consistency.

## 39.5 Base depletion overlap

If the independent-junction approximation predicts:

\[
W_{EB,base}+W_{CB,base}\ge W_B
\]

then the simple two-independent-junction model is no longer valid.

The program shall issue a model-validity warning/error rather than silently producing an apparently valid diagram.

## 39.6 BJT quasi-Fermi limitation

Do not generate full spatial \(E_{Fn}(x)\) or \(E_{Fp}(x)\) profiles without solving carrier transport.

The initial model may show region/terminal quasi-Fermi **reference values** only.

The GUI shall state:

> Quasi-Fermi spatial variation is not solved in the analytical BJT model.

---

# 40. MOS / NMOS — FINAL LOCKED THEORY

## 40.1 Scope correction

The first "NMOS" model is a **1-D MOS electrostatic / MOS-capacitor model**.

It models vertical gate–oxide–semiconductor band bending. It does not model complete lateral MOSFET operation.

It does not solve:

- drain-current characteristics;
- channel-length effects;
- pinch-off;
- saturation;
- source/drain transport;
- \(V_{DS}\)-dependent channel potential.

The GUI should therefore label the device:

> **NMOS (MOS electrostatic model)**

## 40.2 MOS coordinate convention

Set:

\[
x=0
\]

at the Si/SiO₂ interface and:

\[
x>0
\]

into the p-type silicon bulk.

Set:

\[
\phi(x\rightarrow\infty)=0
\]

and:

\[
\boxed{\psi_s=\phi(0)}
\]

For the selected p-substrate NMOS convention:

\[
\boxed{\psi_s>0\Rightarrow\text{downward semiconductor band bending}}
\]

## 40.3 Fermi potential

Define the p-substrate Fermi-potential magnitude as:

\[
\boxed{\phi_F=V_T\ln\left(\frac{N_A}{n_i}\right)>0}
\]

At the strong-inversion threshold:

\[
\boxed{\psi_s=2\phi_F}
\]

## 40.4 Oxide capacitance and flat-band voltage

\[
C'_{ox}=\frac{\epsilon_{ox}}{t_{ox}}
\]

\[
\epsilon_{ox}=3.9\epsilon_0
\]

Use:

\[
\boxed{V_{FB}=\Phi_{MS}-\frac{Q_{ox}}{C'_{ox}}}
\]

with the signs and units of \(\Phi_{MS}\) and \(Q_{ox}\) explicitly documented.

Default:

\[
Q_{ox}=0
\]

## 40.5 Threshold voltage

Lock the Level-1 long-channel expression:

\[
\boxed{
V_T=V_{FB}+2\phi_F+
\frac{\sqrt{4q\epsilon_{si}N_A\phi_F}}{C'_{ox}}
}
\]

## 40.6 Semiconductor charge sign convention

Define \(Q_s\) as net semiconductor charge per unit area.

For the p-substrate convention:

- accumulation: \(Q_s>0\);
- depletion: \(Q_s<0\);
- inversion: \(Q_s<0\).

Use:

\[
\boxed{V_G=V_{FB}+\psi_s-\frac{Q_s(\psi_s)}{C'_{ox}}}
\]

with this sign convention.

## 40.7 Level-1 depletion approximation

Before strong inversion:

\[
W_d=\sqrt{\frac{2\epsilon_{si}\psi_s}{qN_A}}
\]

\[
Q_s=-qN_AW_d
\]

Therefore:

\[
V_G=V_{FB}+\psi_s+\frac{qN_AW_d}{C'_{ox}}
\]

At threshold:

\[
\psi_s=2\phi_F
\]

which recovers the threshold expression.

Do not use the depletion-only charge formula as the complete accumulation or strong-inversion charge model.

## 40.8 Level-2 numerical model

The numerical MOS model shall solve a physically defined \(Q_s(\psi_s)\) relation using a bounded/bracketed root solver such as Brent's method.

The solver must report:

- convergence status;
- residual;
- iteration count;
- bracket/bounds;
- warning/error state.

It must not silently clip \(\psi_s\) to force convergence.

## 40.9 MOS operating states

Use:

```text
VG < VFB     → Accumulation
VG = VFB     → Flat band
VFB < VG < VT → Depletion
VG ≥ VT      → Strong inversion
```

Use a numerical tolerance around equality boundaries.

Weak inversion may be displayed only if the implemented model actually supports that distinction.

## 40.10 MOS energy bands

Use:

\[
E_C(x)=E_{C,bulk}-q\phi(x)
\]

\[
E_V(x)=E_C(x)-E_g
\]

and shift \(E_i\) using the same electrostatic potential.

For \(\psi_s>0\), the semiconductor bands bend downward toward the surface.

The plot must distinguish:

```text
Metal | Oxide | Semiconductor
```

and must not draw silicon energy bands inside the oxide.

---

# 41. STANDARD PROGRAMMING DATA MODEL

All device models shall return a common result object.

Recommended structure:

```text
SimulationResult
├── device
├── model
├── position
├── EC
├── EV
├── Ei
├── EF
├── EFn
├── EFp
├── potential
├── electric_field
├── region
├── calculated_parameters
├── operating_condition
├── warnings
├── assumptions
└── model_valid
```

Use `None` or masked values for physically unavailable quantities. Never use fabricated zero-valued curves to fill missing data.

---

# 42. PROGRAMMING RULES

## 42.1 Physics/GUI separation

`app.py` shall contain no device equations.

It may only:

1. collect inputs;
2. create validated parameter objects;
3. call the physics layer;
4. pass results to visualization;
5. display calculated values and errors.

## 42.2 Strong parameter validation

Use `dataclass`-based parameter objects with validation for:

- doping;
- temperature;
- geometry;
- oxide thickness;
- voltage ranges;
- work-function parameters;
- solver settings.

## 42.3 Explicit unit conversion

Create a dedicated `physics/units.py` module.

Conversions must include:

```text
cm^-3 → m^-3
nm → m
um → m
eV ↔ J
```

Prefer explicit names such as:

```text
NA_cm3
NA_m3
tox_nm
tox_m
energy_eV
energy_J
potential_V
```

## 42.4 No silent physical clipping

Never use clipping merely to make a graph render.

Bad:

```python
value = max(value, 0)
```

when the negative value represents an invalid physical state.

Correct:

```text
validate → reject invalid state → report reason
```

## 42.5 Controlled exceptions

Implement:

```text
InvalidParameterError
ModelValidityError
SolverConvergenceError
UnitConversionError
```

The Streamlit layer shall catch these and display readable messages.

## 42.6 Solver diagnostics

Numerical functions shall return or expose:

```text
converged
iterations
residual
bounds
warning
```

---

# 43. REVISED PROJECT ARCHITECTURE

```text
energy-band/
│
├── app.py
├── launcher.py
│
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
│
├── visualization/
│   ├── __init__.py
│   └── band_plot.py
│
├── tests/
│   ├── test_units.py
│   ├── test_semiconductor.py
│   ├── test_electrostatics.py
│   ├── test_pn.py
│   ├── test_bjt.py
│   ├── test_mos.py
│   ├── test_results.py
│   └── test_validation.py
│
├── pyproject.toml
├── uv.lock
├── README.md
└── Implementation_Plan.md
```

### Shared module responsibilities

`units.py` — all unit conversions.

`validation.py` — parameter and model-validity checks.

`results.py` — standardized simulation-result dataclasses.

`electrostatics.py` — shared electric-field/potential utilities.

`mos.py` — the initial MOS electrostatic model.

---

# 44. REVISED DEVELOPMENT ORDER

The implementation shall follow this order:

```text
1. Lock conventions
        ↓
2. Units + constants
        ↓
3. Semiconductor foundation
        ↓
4. Standard result objects
        ↓
5. Shared electrostatics
        ↓
6. PN equilibrium
        ↓
7. PN bias + quasi-Fermi references
        ↓
8. PN validation
        ↓
9. BJT two-junction electrostatics
        ↓
10. BJT classification + validation
        ↓
11. MOS Level-1 electrostatics
        ↓
12. MOS Level-1 validation
        ↓
13. MOS Level-2 surface-potential solver
        ↓
14. Band visualization
        ↓
15. Streamlit GUI
        ↓
16. Voltage sweep/comparison
        ↓
17. Export
        ↓
18. Full testing
        ↓
19. Packaging
```

Do not build the GUI first.

---

# 45. MANDATORY PHYSICS VALIDATION

## 45.1 PN validation

Verify:

\[
V_{bi}=V_T\ln(N_AN_D/n_i^2)
\]

\[
W=\sqrt{\frac{2\epsilon_{si}}q
\left(\frac1{N_A}+\frac1{N_D}\right)V_D}
\]

\[
x_p+x_n=W
\]

\[
N_Ax_p=N_Dx_n
\]

\[
E_{max}=\frac{qN_Ax_p}{\epsilon_{si}}
=\frac{qN_Dx_n}{\epsilon_{si}}
\]

\[
\Delta\phi=V_D
\]

\[
E_C-E_V=E_g
\]

and, in eV plotting units:

\[
E_{Fn,n}-E_{Fp,p}=V_A
\]

## 45.2 PN limiting cases

```text
VA = 0
→ flat EF

VA > 0
→ smaller depletion width and smaller depletion barrier

VA < 0
→ larger depletion width and larger depletion barrier

VA → 0
→ EFn - EFp → 0
```

## 45.3 BJT validation

Test:

- E-B built-in potential;
- B-C built-in potential;
- both depletion widths;
- junction polarity;
- operating-region classification;
- base depletion overlap detection;
- equilibrium flat Fermi level;
- absence of fabricated spatial quasi-Fermi profiles.

## 45.4 MOS validation

Test:

\[
C'_{ox}=\epsilon_{ox}/t_{ox}
\]

\[
\phi_F=V_T\ln(N_A/n_i)
\]

\[
V_{FB}=\Phi_{MS}-Q_{ox}/C'_{ox}
\]

\[
V_T=V_{FB}+2\phi_F+
\frac{\sqrt{4q\epsilon_{si}N_A\phi_F}}{C'_{ox}}
\]

and:

```text
VG < VFB      → accumulation
VG = VFB      → flat band
VFB < VG < VT → depletion
VG ≥ VT       → strong inversion
```

Also verify:

\[
t_{ox}\uparrow\Rightarrow C'_{ox}\downarrow
\]

and:

\[
\psi_s>0\Rightarrow\text{downward band bending}
\]

---

# 46. MANDATORY PLOT-INTEGRITY TESTS

Every valid result must satisfy:

1. EC, EV and Ei use the same electrostatic potential.
2. EC − EV = Eg within tolerance.
3. No NaN/Inf values appear in valid regions.
4. Region boundaries match the model geometry.
5. Depletion boundaries match calculated widths.
6. Equilibrium EF is flat.
7. Biased PN has no fabricated global equilibrium EF.
8. Quasi-Fermi curves are drawn only where supported by the model.
9. Positive MOS surface potential produces downward band bending.
10. The vertical energy reference is deterministic.

---

# 47. GUI MODEL-TRANSPARENCY RULE

The GUI shall show the model being used.

PN:

```text
Model: Abrupt 1-D depletion approximation
```

BJT:

```text
Model: Two-junction analytical electrostatic approximation
```

NMOS:

```text
Model: 1-D MOS electrostatic model
```

For biased PN/BJT cases, show:

```text
Quasi-Fermi representation:
Neutral-region references only
```

unless a future transport model explicitly calculates spatial quasi-Fermi levels.

---

# 48. EXPORT SPECIFICATION

CSV export shall use explicit units in column names:

```text
position_m
EC_eV
EV_eV
Ei_eV
EF_eV
EFn_eV
EFp_eV
potential_V
electric_field_V_per_m
region
```

Unavailable quantities must be blank or represented by the documented missing-value convention, not misleading zeros.

Metadata shall include:

```text
device
model
temperature_K
applied_voltage_V
assumptions
```

PNG export shall include the device name, applied voltage, and operating condition in the figure metadata/title where practical.

---

# 49. REVISED ACCEPTANCE CRITERIA

The simulator is complete only when:

- all sign conventions are explicit;
- all energy references are deterministic;
- all internal units are consistent;
- PN depletion calculations pass analytical validation;
- PN equilibrium has one flat EF;
- biased PN does not show a false global equilibrium EF;
- PN quasi-Fermi references have the correct signed separation;
- no unsupported depletion-region quasi-Fermi profile is drawn;
- BJT junction voltage polarity is explicit;
- BJT operating-region classification is automatic;
- BJT depletion overlap is detected;
- BJT quasi-Fermi levels are not falsely presented as transport solutions;
- the initial NMOS model is explicitly identified as MOS electrostatic rather than full MOSFET transport;
- MOS Qs, ψs, VFB and VT conventions are locked;
- positive ψs produces downward band bending under the selected convention;
- Level-1 MOS results are validated before Level-2 solving;
- invalid physical states are rejected rather than silently clipped;
- numerical diagnostics are available;
- bandgap and potential/energy consistency tests pass;
- voltage sweeps preserve all sign conventions;
- GUI contains no device equations;
- exports preserve units and unavailable values correctly;
- the packaged Windows executable runs independently of the development environment.

---

# 50. FINAL IMPLEMENTATION PRINCIPLES

When a visually attractive result conflicts with calculated physics:

\[
\boxed{\text{Physics wins.}}
\]

When an approximation cannot calculate a quantity:

\[
\boxed{\text{Do not invent it.}}
\]

When an input is outside the model's validity:

\[
\boxed{\text{Reject it explicitly.}}
\]

When a sign convention is ambiguous:

\[
\boxed{\text{Do not implement until it is defined.}}
\]

The project priority is:

\[
\boxed{
\text{Correct physics}
>
\text{Correct conventions}
>
\text{Numerical validation}
>
\text{Visualization}
>
\text{GUI polish}
}
\]

The completed application is a **physics-first, one-dimensional analytical semiconductor energy-band visualization tool** for PN junctions, NPN BJTs, and NMOS/MOS electrostatic structures. It does not claim to be a TCAD replacement or a full carrier-transport simulator unless those equations are explicitly implemented.
