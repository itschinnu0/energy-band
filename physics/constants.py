"""Physical constants for semiconductor physics.

Reference: S. M. Sze and Kwok K. Ng, 'Physics of Semiconductor Devices', 3rd ed.
All constants in SI units unless explicitly noted with units in variable name.
"""

# Fundamental Physical Constants (CODATA / Sze)
Q: float = 1.602176634e-19  # Elementary charge [C]
K_B: float = 1.380649e-23  # Boltzmann constant [J/K]
K_B_EV: float = 8.617333262e-5  # Boltzmann constant [eV/K]
EPSILON_0: float = 8.8541878128e-12  # Vacuum permittivity [F/m]
H_PLANCK: float = 6.62607015e-34  # Planck constant [J*s]
M_0: float = 9.1093837015e-31  # Free electron rest mass [kg]

# Silicon Material Parameters (Sze, 3rd ed., Appendix G & Chap 1)
EPS_R_SI: float = 11.7  # Silicon relative dielectric permittivity
EPS_SI: float = EPS_R_SI * EPSILON_0  # Silicon permittivity [F/m]

# Silicon Bandgap Varshni parameters (Sze): Eg(T) = Eg0 - alpha * T^2 / (T + beta)
EG_0_EV: float = 1.166  # Bandgap at 0 K [eV]
VARSHNI_ALPHA_EV: float = 4.73e-4  # Varshni parameter alpha [eV/K]
VARSHNI_BETA_K: float = 636.0  # Varshni parameter beta [K]

# Silicon Effective Density of States at 300 K (Sze, App G)
# In cm^-3 and m^-3:
NC_300_CM3: float = 2.86e19  # Effective conduction band DOS at 300 K [cm^-3]
NV_300_CM3: float = 2.66e19  # Effective valence band DOS at 300 K [cm^-3]
NC_300_M3: float = NC_300_CM3 * 1e6  # [m^-3]
NV_300_M3: float = NV_300_CM3 * 1e6  # [m^-3]

# Reference room temperature
T_REF_K: float = 300.0  # Room temperature [K]

# Silicon dioxide parameters (Sze, for MOS reference)
EPS_R_OX: float = 3.9  # SiO2 relative dielectric permittivity
EPS_OX: float = EPS_R_OX * EPSILON_0  # SiO2 permittivity [F/m]
