"""1-D abrupt PN junction physics model under the depletion approximation.

Conventions:
    V_A = V_P - V_N (V_A > 0 is forward bias, V_A < 0 is reverse bias)
    V_bi = V_T * ln(N_A * N_D / n_i^2)
    V_D = V_bi - V_A
    Valid only when V_D > 0. If V_D <= 0, raises ModelValidityError.

Depletion parameters:
    W = sqrt( (2 * eps_si / q) * (1/N_A + 1/N_D) * V_D )
    x_p = N_D / (N_A + N_D) * W
    x_n = N_A / (N_A + N_D) * W
    Charge neutrality: N_A * x_p = N_D * x_n

Depletion region coordinates:
    Junction metallurgical interface at x = 0.
    Depletion extends from -x_p to +x_n.
    p-side neutral: x < -x_p
    n-side neutral: x > x_n

Charge density:
    rho(x) = 0           for x < -x_p
    rho(x) = -q * N_A    for -x_p <= x <= 0
    rho(x) = +q * N_D    for 0 < x <= x_n
    rho(x) = 0           for x > x_n

Electric field:
    E(x) = - (q * N_A / eps_si) * (x + x_p)   for -x_p <= x <= 0
    E(x) = - (q * N_D / eps_si) * (x_n - x)   for 0 <= x <= x_n
    E(x) = 0 in neutral regions.
    Peak magnitude at x = 0: E_max = q * N_A * x_p / eps_si = q * N_D * x_n / eps_si.
    E(x) <= 0 throughout depletion region with this coordinate convention.

Electrostatic potential phi(x):
    Choosing phi(-x_p) = 0 V (p-neutral reference):
    For x <= -x_p: phi(x) = 0
    For -x_p <= x <= 0: phi(x) = (q * N_A / (2 * eps_si)) * (x + x_p)^2
    For 0 <= x <= x_n: phi(x) = V_D - (q * N_D / (2 * eps_si)) * (x_n - x)^2
    For x >= x_n: phi(x) = V_D
    Total potential drop across depletion: Delta phi = phi(x_n) - phi(-x_p) = V_D.

Bands:
    E_C(x) = E_C_ref - phi(x)
    E_V(x) = E_C(x) - E_g
    E_i(x) = E_C(x) - E_g / 2

Fermi / Quasi-Fermi references:
    Equilibrium (V_A == 0):
        Flat global E_F across entire device:
        E_F = E_i(p-bulk) - V_T * ln(N_A / n_i) = E_i(n-bulk) + V_T * ln(N_D / n_i)
    Biased (V_A != 0):
        No global E_F.
        EFp in neutral p region (x < -x_p): EFp = -V_A eV (or relative offset matching EFn_n - EFp_p = V_A).
        EFn in neutral n region (x > x_n): EFn = 0.0 eV.
        Depletion region has no fabricated quasi-Fermi profile (NaN/None).
"""

import math
from typing import Optional, Tuple
import numpy as np

from physics.constants import Q, EPS_SI, K_B_EV
from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.validation import validate_temperature, validate_doping
from physics.units import cm3_to_m3
from physics.semiconductor import (
    thermal_voltage,
    silicon_bandgap_ev,
    intrinsic_carrier_concentration,
    bands_from_potential,
)
from physics.results import SimulationResult


def calculate_pn_built_in_potential(
    na_cm3: float,
    nd_cm3: float,
    temp_k: float = 300.0,
) -> float:
    """Calculate PN junction built-in potential Vbi in Volts [V].

    Vbi = V_T * ln(NA * ND / ni^2)
    """
    validate_temperature(temp_k)
    validate_doping(na_cm3, doping_name="Acceptor (NA)")
    validate_doping(nd_cm3, doping_name="Donor (ND)")

    vt = thermal_voltage(temp_k)
    ni = intrinsic_carrier_concentration(temp_k)
    vbi = vt * math.log((na_cm3 * nd_cm3) / (ni**2))
    return vbi


def calculate_pn_depletion_widths(
    na_cm3: float,
    nd_cm3: float,
    vd_v: float,
    epsilon: float = EPS_SI,
) -> Tuple[float, float, float]:
    """Calculate total depletion width W and regional widths (xp, xn) in meters [m].

    W = sqrt( (2 * eps / q) * (1/NA + 1/ND) * VD )
    xp = ND / (NA + ND) * W
    xn = NA / (NA + ND) * W

    Args:
        na_cm3: Acceptor doping in cm^-3.
        nd_cm3: Donor doping in cm^-3.
        vd_v: Effective barrier voltage VD = Vbi - VA in Volts [V]. Must be > 0.
        epsilon: Permittivity in F/m.

    Returns:
        (W, xp, xn) in meters.

    Raises:
        ModelValidityError: If vd_v <= 0 (abrupt depletion approximation invalid).
    """
    validate_doping(na_cm3, doping_name="Acceptor (NA)")
    validate_doping(nd_cm3, doping_name="Donor (ND)")

    if vd_v <= 0.0:
        raise ModelValidityError(
            f"Effective barrier voltage VD = {vd_v:.4e} V <= 0. "
            "The abrupt depletion approximation is invalid for VD <= 0 (flat-band or strong forward bias)."
        )

    na_m3 = cm3_to_m3(na_cm3)
    nd_m3 = cm3_to_m3(nd_cm3)

    w = math.sqrt((2.0 * epsilon / Q) * (1.0 / na_m3 + 1.0 / nd_m3) * vd_v)
    xp = (nd_m3 / (na_m3 + nd_m3)) * w
    xn = (na_m3 / (na_m3 + nd_m3)) * w

    return w, xp, xn


def simulate_pn_junction(
    na_cm3: float = 1.0e16,
    nd_cm3: float = 1.0e16,
    va_v: float = 0.0,
    temp_k: float = 300.0,
    lp_m: Optional[float] = None,
    ln_m: Optional[float] = None,
    num_points: int = 500,
) -> SimulationResult:
    """Simulate 1-D abrupt PN junction and return SimulationResult.

    Args:
        na_cm3: Acceptor doping concentration in cm^-3.
        nd_cm3: Donor doping concentration in cm^-3.
        va_v: Applied voltage VA = VP - VN in Volts [V].
        temp_k: Temperature in Kelvin [K].
        lp_m: Optional neutral p-region thickness in meters. If None, defaults to 3 * xp.
        ln_m: Optional neutral n-region thickness in meters. If None, defaults to 3 * xn.
        num_points: Number of spatial grid points.

    Returns:
        SimulationResult containing band diagram and electrostatics profiles.
    """
    validate_temperature(temp_k)
    validate_doping(na_cm3, doping_name="Acceptor (NA)")
    validate_doping(nd_cm3, doping_name="Donor (ND)")

    if not math.isfinite(va_v):
        raise InvalidParameterError(f"Applied voltage VA must be finite, got: {va_v}")

    # Physics constants and properties
    eg_ev = silicon_bandgap_ev(temp_k)
    vbi = calculate_pn_built_in_potential(na_cm3, nd_cm3, temp_k)
    vd = vbi - va_v

    # Depletion approximation validation
    if vd <= 0.0:
        raise ModelValidityError(
            f"Effective barrier voltage VD = Vbi - VA = {vbi:.4f} V - {va_v:.4f} V = {vd:.4e} V <= 0. "
            "Abrupt depletion approximation is outside its physical validity range."
        )

    w, xp, xn = calculate_pn_depletion_widths(na_cm3, nd_cm3, vd, epsilon=EPS_SI)

    # Geometry boundaries
    # Junction is at x = 0. Depletion is [-xp, xn].
    # Neutral regions extend to -wp_total and +wn_total
    wp_neutral = lp_m if lp_m is not None else 3.0 * xp
    wn_neutral = ln_m if ln_m is not None else 3.0 * xn

    if wp_neutral <= 0.0 or wn_neutral <= 0.0:
        raise InvalidParameterError("Neutral region thicknesses must be positive.")

    x_min = - (wp_neutral + xp)
    x_max = wn_neutral + xn

    # Construct spatial grid
    # Include explicit points right at -xp, 0, xn for sharp depletion boundary transitions
    grid_p_neutral = np.linspace(x_min, -xp, max(int(num_points * 0.25), 20), endpoint=False)
    grid_p_depletion = np.linspace(-xp, 0.0, max(int(num_points * 0.25), 20), endpoint=False)
    grid_n_depletion = np.linspace(0.0, xn, max(int(num_points * 0.25), 20), endpoint=True)
    grid_n_neutral = np.linspace(xn, x_max, max(int(num_points * 0.25), 20) + 1)[1:]

    x = np.concatenate([grid_p_neutral, grid_p_depletion, grid_n_depletion, grid_n_neutral])
    x = np.unique(x)

    # Analytical Charge density rho(x) [C/m^3]
    na_m3 = cm3_to_m3(na_cm3)
    nd_m3 = cm3_to_m3(nd_cm3)
    rho = np.zeros_like(x)
    p_dep_mask = (x >= -xp) & (x < 0.0)
    n_dep_mask = (x >= 0.0) & (x <= xn)

    rho[p_dep_mask] = -Q * na_m3
    rho[n_dep_mask] = Q * nd_m3

    # Analytical Electric field E(x) [V/m]
    # E(-xp) = 0, E(xn) = 0. Peak magnitude at x = 0
    electric_field = np.zeros_like(x)
    electric_field[p_dep_mask] = -(Q * na_m3 / EPS_SI) * (x[p_dep_mask] + xp)
    electric_field[n_dep_mask] = -(Q * nd_m3 / EPS_SI) * (xn - x[n_dep_mask])
    peak_e_field = float(np.min(electric_field))

    # Analytical Electrostatic potential phi(x) [V]
    # Reference: phi(-xp) = 0 V in the p-neutral region
    potential = np.zeros_like(x)
    potential[x < -xp] = 0.0
    potential[p_dep_mask] = (Q * na_m3 / (2.0 * EPS_SI)) * (x[p_dep_mask] + xp)**2
    potential[n_dep_mask] = vd - (Q * nd_m3 / (2.0 * EPS_SI)) * (xn - x[n_dep_mask])**2
    potential[x > xn] = vd

    potential_drop = float(potential[-1] - potential[0])

    # Region classification array
    region = np.empty(x.shape, dtype=object)
    region[x < -xp] = "p-neutral"
    region[p_dep_mask] = "p-depletion"
    region[n_dep_mask] = "n-depletion"
    region[x > xn] = "n-neutral"

    # Energy bands [eV]
    # Reference convention: EC = EC_ref - phi
    vt_ev = K_B_EV * temp_k
    ni = intrinsic_carrier_concentration(temp_k)
    delta_ef_p = vt_ev * math.log(na_cm3 / ni)  # Ei - EF in p-type

    is_equilibrium = math.isclose(va_v, 0.0, abs_tol=1e-9)

    # Locked quasi-Fermi references:
    # EFn_n = 0.0 eV
    # EFp_p = -VA eV
    # EFn_n - EFp_p = VA eV
    # In equilibrium (VA = 0), EF = 0.0 eV flat everywhere.
    # In p-bulk, EFp = -VA. Since EFp = Ei_p - delta_ef_p, Ei_p = -VA + delta_ef_p.
    # Since Ei = EC - Eg/2, EC_p = Ei_p + Eg/2 = -VA + delta_ef_p + Eg/2.
    # Since phi(p) = 0 and EC(p) = EC_ref - phi(p), we set:
    ec_ref_ev = -va_v + delta_ef_p + 0.5 * eg_ev

    ec, ev, ei = bands_from_potential(potential, eg_ev=eg_ev, ec_ref_ev=ec_ref_ev)

    # Fermi / Quasi-Fermi levels
    ef: Optional[np.ndarray] = None
    efn: Optional[np.ndarray] = None
    efp: Optional[np.ndarray] = None

    if is_equilibrium:
        # At equilibrium, single flat EF = 0.0 eV across the whole device
        ef = np.zeros_like(x)
        operating_condition = "Equilibrium (VA = 0 V)"
    else:
        # Biased PN junction: DO NOT fabricate a global flat EF or depletion profile
        ef = None
        efn = np.full_like(x, fill_value=np.nan)
        efp = np.full_like(x, fill_value=np.nan)

        # Only define in neutral regions where model is physically justified
        efp[x < -xp] = -va_v
        efn[x > xn] = 0.0

        if va_v > 0:
            operating_condition = f"Forward bias (VA = {va_v:.3f} V > 0)"
        else:
            operating_condition = f"Reverse bias (VA = {va_v:.3f} V < 0)"

    calculated_parameters = {
        "Vbi_V": vbi,
        "VD_V": vd,
        "VA_V": va_v,
        "W_m": w,
        "xp_m": xp,
        "xn_m": xn,
        "peak_electric_field_V_per_m": peak_e_field,
        "potential_drop_V": potential_drop,
        "temperature_K": temp_k,
        "NA_cm3": na_cm3,
        "ND_cm3": nd_cm3,
        "bandgap_eV": eg_ev,
        "built_in_potential_V": vbi,
        "effective_barrier_V": vd,
        "depletion_width_m": w,
    }

    assumptions = [
        "1-D abrupt junction approximation",
        "Full depletion approximation with sharp boundaries at -xp and +xn",
        "Complete dopant ionization",
        "Boltzmann carrier statistics (non-degenerate silicon)",
        "Uniform doping in P and N regions",
        "Neutral regions have zero electric field and flat potential",
        "Quasi-Fermi levels are specified only in neutral regions (no full transport solution)",
    ]

    warnings = []
    if va_v > 0.5 * vbi:
        warnings.append(
            f"Forward bias VA = {va_v:.3f} V is approaching built-in potential Vbi = {vbi:.3f} V. "
            "High-level injection and series resistance effects are not included in the abrupt depletion model."
        )

    return SimulationResult(
        device="PN Junction",
        model="Abrupt 1-D depletion approximation",
        position=x,
        EC=ec,
        EV=ev,
        Ei=ei,
        EF=ef,
        EFn=efn,
        EFp=efp,
        potential=potential,
        electric_field=electric_field,
        region=region,
        calculated_parameters=calculated_parameters,
        operating_condition=operating_condition,
        warnings=warnings,
        assumptions=assumptions,
        model_valid=True,
    )
