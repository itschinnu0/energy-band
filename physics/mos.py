"""1-D vertical MOS electrostatics and NMOS band model.

Structure:
    Metal | SiO2 (oxide) | p-type Silicon

Coordinate:
    x = 0 at the Si/SiO2 interface
    x > 0 into silicon bulk

Conventions:
    phi_bulk = 0 V
    psi_s = phi(0) (surface potential)
    psi_s > 0 means downward band bending in p-type silicon.
    Band bending: EC(x) = EC_bulk - phi(x)
                  EV(x) = EC(x) - Eg
                  Ei(x) = EC(x) - Eg/2
    Semiconductor bands are NOT drawn inside the oxide (x < 0).

Parameters & Equations:
    Oxide capacitance:
        Cox = epsilon_ox / tox [F/m^2]
        epsilon_ox = 3.9 * epsilon_0

    Fermi potential for p-type silicon:
        phi_F = VT * ln(NA / ni) [V]
        where VT = kT/q. Require phi_F > 0.

    Flat-band voltage:
        VFB = PhiMS - Qox / Cox [V]
        PhiMS: metal-semiconductor work function difference in Volts [V].
        Qox: equivalent oxide fixed charge per unit area in C/m^2.

    Level-1 Threshold voltage:
        threshold_voltage = VFB + 2*phi_F + sqrt(4*q*epsilon_si*NA*phi_F) / Cox [V]

    Semiconductor charge convention:
        Qs > 0 : accumulation (psi_s < 0)
        Qs < 0 : depletion and inversion (psi_s > 0)
        Depletion approximation (0 < psi_s <= 2*phi_F):
            Wd = sqrt(2 * epsilon_si * psi_s / (q * NA))
            Qs = - q * NA * Wd

    Gate-voltage relation:
        VG = VFB + psi_s - Qs / Cox

    Operating state classification:
        tol = 1e-3 V (1 mV) around equality
        VG < VFB - tol                  -> Accumulation
        |VG - VFB| <= tol               -> Flat band
        VFB + tol < VG < VT - tol       -> Depletion
        |VG - VT| <= tol                -> Threshold / transition
        VG > VT + tol                   -> Strong inversion

Level-2 Numerical Solving:
    Exact analytical semiconductor charge Qs(psi_s) from Poisson-Boltzmann (Sze):
        F(psi_s, phi_F, VT) = sqrt(2) * ( (exp(-psi_s/VT) + psi_s/VT - 1)
                                      + (ni/NA)^2 * (exp(psi_s/VT) - psi_s/VT - 1) )^0.5
        Qs(psi_s) = - sgn(psi_s) * (2*q*eps_si*NA*VT)^0.5 * F(...)
    Bounded root solving for psi_s:
        f(psi_s) = VFB + psi_s - Qs(psi_s)/Cox - VG == 0
"""

import math
from typing import Optional, Tuple, Dict, Any
import numpy as np

from physics.constants import Q, EPS_SI, EPS_OX, EPSILON_0, K_B_EV
from physics.exceptions import InvalidParameterError, ModelValidityError, SolverConvergenceError
from physics.validation import validate_temperature, validate_doping
from physics.units import cm3_to_m3
from physics.semiconductor import (
    thermal_voltage,
    silicon_bandgap_ev,
    intrinsic_carrier_concentration,
    bands_from_potential,
)
from physics.results import SimulationResult
from physics.solver import solve_bracketed_root


def calculate_oxide_capacitance(tox_m: float, epsilon_ox: float = EPS_OX) -> float:
    """Calculate oxide capacitance per unit area Cox in Farads per square meter [F/m^2].

    Cox = epsilon_ox / tox_m

    Args:
        tox_m: Oxide physical thickness in meters [m].
        epsilon_ox: Oxide dielectric permittivity in F/m (default: 3.9 * eps0).

    Returns:
        Cox in F/m^2.

    Raises:
        InvalidParameterError: If tox_m <= 0 or not finite.
    """
    if not math.isfinite(tox_m):
        raise InvalidParameterError(f"Oxide thickness must be finite, got: {tox_m}")
    if tox_m <= 0.0:
        raise InvalidParameterError(f"Oxide thickness must be strictly positive (tox > 0), got: {tox_m} m")
    return epsilon_ox / tox_m


def calculate_fermi_potential(na_cm3: float, temp_k: float = 300.0) -> float:
    """Calculate bulk Fermi potential phi_F in Volts [V] for p-type silicon.

    phi_F = (kT/q) * ln(NA / ni)

    Args:
        na_cm3: Substrate acceptor doping in cm^-3.
        temp_k: Temperature in Kelvin [K].

    Returns:
        phi_F in Volts (> 0 for p-type).
    """
    validate_temperature(temp_k)
    validate_doping(na_cm3, doping_name="Substrate acceptor (NA)")

    vt = thermal_voltage(temp_k)
    ni = intrinsic_carrier_concentration(temp_k)
    phi_f = vt * math.log(na_cm3 / ni)

    if phi_f <= 0.0:
        raise ModelValidityError(
            f"Calculated Fermi potential phi_F = {phi_f:.4e} V <= 0. "
            "P-type substrate requires NA > ni so phi_F > 0."
        )
    return phi_f


def calculate_flat_band_voltage(phims_v: float, qox_c_m2: float, cox_f_m2: float) -> float:
    """Calculate flat-band voltage VFB in Volts [V].

    VFB = PhiMS - Qox / Cox

    Args:
        phims_v: Metal-semiconductor work function difference in Volts [V].
        qox_c_m2: Equivalent oxide fixed charge per unit area in C/m^2.
        cox_f_m2: Oxide capacitance per unit area in F/m^2.

    Returns:
        VFB in Volts.
    """
    if not math.isfinite(phims_v):
        raise InvalidParameterError(f"PhiMS must be finite, got: {phims_v}")
    if not math.isfinite(qox_c_m2):
        raise InvalidParameterError(f"Qox must be finite, got: {qox_c_m2}")
    if not math.isfinite(cox_f_m2) or cox_f_m2 <= 0.0:
        raise InvalidParameterError(f"Cox must be positive finite, got: {cox_f_m2}")

    return phims_v - (qox_c_m2 / cox_f_m2)


def calculate_level1_threshold_voltage(
    vfb_v: float,
    phi_f_v: float,
    na_cm3: float,
    cox_f_m2: float,
    epsilon_si: float = EPS_SI,
) -> float:
    """Calculate Level-1 MOSFET threshold voltage VT in Volts [V] for p-type substrate.

    threshold_voltage = VFB + 2*phi_F + sqrt(4 * q * epsilon_si * NA * phi_F) / Cox

    Args:
        vfb_v: Flat-band voltage in Volts [V].
        phi_f_v: Bulk Fermi potential in Volts [V] (> 0).
        na_cm3: Acceptor doping in cm^-3.
        cox_f_m2: Oxide capacitance in F/m^2.
        epsilon_si: Silicon permittivity in F/m.

    Returns:
        Threshold voltage in Volts.
    """
    if not math.isfinite(vfb_v):
        raise InvalidParameterError(f"VFB must be finite, got: {vfb_v}")
    if phi_f_v <= 0.0 or not math.isfinite(phi_f_v):
        raise InvalidParameterError(f"phi_F must be positive finite, got: {phi_f_v}")
    validate_doping(na_cm3, doping_name="Substrate acceptor (NA)")
    if cox_f_m2 <= 0.0 or not math.isfinite(cox_f_m2):
        raise InvalidParameterError(f"Cox must be positive finite, got: {cox_f_m2}")

    na_m3 = cm3_to_m3(na_cm3)
    q_dep_max = math.sqrt(4.0 * Q * epsilon_si * na_m3 * phi_f_v)
    v_th = vfb_v + 2.0 * phi_f_v + (q_dep_max / cox_f_m2)
    return v_th


def calculate_depletion_width(psi_s_v: float, na_cm3: float, epsilon_si: float = EPS_SI) -> float:
    """Calculate depletion width Wd in meters [m] under the depletion approximation.

    Wd = sqrt(2 * epsilon_si * psi_s / (q * NA))

    Valid only for psi_s >= 0 (depletion / inversion).

    Args:
        psi_s_v: Surface potential in Volts [V] (must be >= 0).
        na_cm3: Acceptor doping in cm^-3.
        epsilon_si: Silicon permittivity in F/m.

    Returns:
        Depletion width Wd in meters.

    Raises:
        ModelValidityError: If psi_s_v < 0 (accumulation has no depletion layer).
    """
    validate_doping(na_cm3, doping_name="Substrate acceptor (NA)")
    if psi_s_v < 0.0:
        raise ModelValidityError(
            f"Surface potential psi_s = {psi_s_v:.4e} V < 0 (accumulation). "
            "Depletion width formula is physically invalid for negative psi_s."
        )
    na_m3 = cm3_to_m3(na_cm3)
    return math.sqrt(2.0 * epsilon_si * psi_s_v / (Q * na_m3))


def classify_mos_operating_state(
    vg_v: float,
    vfb_v: float,
    threshold_voltage_v: float,
    tol_v: float = 1e-3,
) -> str:
    """Classify MOS operating state based on applied gate voltage VG.

    Args:
        vg_v: Applied gate voltage in Volts [V].
        vfb_v: Flat-band voltage in Volts [V].
        threshold_voltage_v: Threshold voltage in Volts [V].
        tol_v: Voltage tolerance deadband around boundary points [V].

    Returns:
        One of: 'Accumulation', 'Flat band', 'Depletion', 'Threshold / transition', 'Strong inversion'.
    """
    if not (math.isfinite(vg_v) and math.isfinite(vfb_v) and math.isfinite(threshold_voltage_v)):
        raise InvalidParameterError("Voltages must be finite numbers.")

    if abs(vg_v - vfb_v) <= tol_v:
        return "Flat band"
    elif vg_v < vfb_v - tol_v:
        return "Accumulation"
    elif abs(vg_v - threshold_voltage_v) <= tol_v:
        return "Threshold / transition"
    elif vg_v > threshold_voltage_v + tol_v:
        return "Strong inversion"
    else:
        return "Depletion"


def calculate_semiconductor_charge_level2(
    psi_s_v: float,
    na_cm3: float,
    phi_f_v: float,
    temp_k: float = 300.0,
    epsilon_si: float = EPS_SI,
) -> float:
    """Calculate exact total semiconductor space-charge per unit area Qs(psi_s) in C/m^2.

    Follows S. M. Sze & K. K. Ng, Physics of Semiconductor Devices, 3rd ed., Chapter 4:
        Qs = - sgn(psi_s) * sqrt(2 * eps_si * k * T * NA) * F(psi_s, phi_F, VT)
    where:
        F(u, uF) = [ (e^{-u} + u - 1) + (ni/NA)^2 * (e^{u} - u - 1) ]^(1/2)
        u = psi_s / VT
        (ni/NA)^2 = exp(-2*phi_F / VT)

    Sign convention:
        Qs > 0 for psi_s < 0 (accumulation of holes)
        Qs < 0 for psi_s > 0 (depletion of holes & inversion of electrons)
        Qs = 0 for psi_s = 0 (flat-band)

    Args:
        psi_s_v: Surface potential in Volts [V].
        na_cm3: Acceptor concentration in cm^-3.
        phi_f_v: Bulk Fermi potential in Volts [V].
        temp_k: Temperature in Kelvin [K].
        epsilon_si: Silicon permittivity in F/m.

    Returns:
        Qs in Coulombs per square meter [C/m^2].
    """
    if math.isclose(psi_s_v, 0.0, abs_tol=1e-12):
        return 0.0

    vt = thermal_voltage(temp_k)
    u = psi_s_v / vt
    u_f = phi_f_v / vt
    na_m3 = cm3_to_m3(na_cm3)

    # (ni/NA)^2 = exp(-2 * u_f)
    exp_minus_2uf = math.exp(-2.0 * u_f)

    # To avoid catastrophic floating-point cancellation for small u,
    # note that exp(-u) + u - 1 = u^2/2 + O(u^3)
    # math.expm1(x) = exp(x) - 1.
    # exp(-u) + u - 1 = expm1(-u) + u
    # exp(u) - u - 1 = expm1(u) - u
    if abs(u) < 1e-4:
        term_holes = 0.5 * u**2 - (u**3) / 6.0 + (u**4) / 24.0
        term_electrons = 0.5 * u**2 + (u**3) / 6.0 + (u**4) / 24.0
    else:
        # Cap large positive u to avoid overflow in exp(u)
        if u > 80.0:
            term_electrons = math.exp(u)
            term_holes = u - 1.0
        elif u < -80.0:
            term_holes = math.exp(-u)
            term_electrons = -u - 1.0
        else:
            term_holes = math.expm1(-u) + u
            term_electrons = math.expm1(u) - u

    inside_sqrt = term_holes + exp_minus_2uf * term_electrons
    if inside_sqrt < 0.0:
        inside_sqrt = 0.0

    f_val = math.sqrt(inside_sqrt)
    prefactor = math.sqrt(2.0 * epsilon_si * Q * na_m3 * vt)

    sign = -1.0 if psi_s_v > 0.0 else 1.0
    return sign * prefactor * f_val


def solve_surface_potential_level2(
    vg_v: float,
    vfb_v: float,
    cox_f_m2: float,
    na_cm3: float,
    phi_f_v: float,
    temp_k: float = 300.0,
    epsilon_si: float = EPS_SI,
    bracket_range: Optional[Tuple[float, float]] = None,
) -> Tuple[float, Dict[str, Any]]:
    """Numerically solve for surface potential psi_s in Volts [V] for a given gate voltage VG.

    Equation:
        f(psi_s) = VFB + psi_s - Qs(psi_s) / Cox - VG == 0

    Args:
        vg_v: Gate voltage in Volts [V].
        vfb_v: Flat-band voltage in Volts [V].
        cox_f_m2: Oxide capacitance in F/m^2.
        na_cm3: Acceptor doping in cm^-3.
        phi_f_v: Fermi potential in Volts [V].
        temp_k: Temperature in Kelvin [K].
        epsilon_si: Permittivity in F/m.
        bracket_range: Optional explicit (psi_min, psi_max) search bracket.

    Returns:
        (psi_s, solver_diagnostics)

    Raises:
        SolverConvergenceError: If bracketing or root finding fails.
    """
    def obj_func(psi_s: float) -> float:
        qs = calculate_semiconductor_charge_level2(
            psi_s_v=psi_s,
            na_cm3=na_cm3,
            phi_f_v=phi_f_v,
            temp_k=temp_k,
            epsilon_si=epsilon_si,
        )
        return (vfb_v + psi_s - (qs / cox_f_m2)) - vg_v

    # Determine bracketing interval
    # In accumulation (VG < VFB): psi_s is negative.
    # In depletion/inversion (VG > VFB): psi_s is positive.
    # Lower bound cannot be much less than ~ -0.8 V (strong accumulation).
    # Upper bound in strong inversion rarely exceeds ~ 2*phi_F + 0.5 V or bandgap.
    eg_ev = silicon_bandgap_ev(temp_k)

    if bracket_range is not None:
        bracket = bracket_range
    else:
        # Dynamic search for bracket containing opposite signs
        # psi_s must lie inside [-1.5, eg_ev + 0.5]
        low = -1.5
        high = max(2.5 * phi_f_v + 0.5, eg_ev + 0.5)

        # Ensure low and high bracket the root
        # If vg is very large positive or negative, widen gracefully
        f_low = obj_func(low)
        f_high = obj_func(high)

        # Check if 0 is between them
        f_zero = obj_func(0.0)

        if f_low * f_zero <= 0:
            bracket = (low, 0.0)
        elif f_zero * f_high <= 0:
            bracket = (0.0, high)
        elif f_low * f_high <= 0:
            bracket = (low, high)
        else:
            # Try expanded search
            if vg_v < vfb_v:
                bracket = (-3.0, 0.0)
            else:
                bracket = (0.0, high + 1.0)

    psi_s_sol, diag = solve_bracketed_root(
        obj_func,
        bracket=bracket,
        xtol=1e-11,
        rtol=1e-9,
        maxiter=120,
        name="Surface potential psi_s",
    )

    return psi_s_sol, diag


def simulate_mos(
    na_cm3: float = 1.0e16,
    tox_m: float = 10.0e-9,
    vg_v: float = 0.0,
    phims_v: float = -0.9,
    qox_c_m2: float = 0.0,
    temp_k: float = 300.0,
    level: int = 1,
    si_thickness_m: Optional[float] = None,
    num_points: int = 400,
) -> SimulationResult:
    """Simulate 1-D vertical MOS electrostatics and band diagram.

    Metal | SiO2 (x < 0) | p-type Silicon (x >= 0)
    x = 0 is the Si/SiO2 interface.

    Args:
        na_cm3: Substrate acceptor doping in cm^-3.
        tox_m: Oxide thickness in meters [m].
        vg_v: Applied gate voltage in Volts [V].
        phims_v: Metal-semiconductor work function difference in Volts [V].
        qox_c_m2: Equivalent oxide fixed charge in C/m^2.
        temp_k: Temperature in Kelvin [K].
        level: Simulation level (1 for analytical depletion approx, 2 for numerical solver).
        si_thickness_m: Total silicon domain depth in meters.
        num_points: Number of spatial grid points in silicon.

    Returns:
        SimulationResult containing MOS band diagram and parameters.
    """
    validate_temperature(temp_k)
    validate_doping(na_cm3, doping_name="Substrate acceptor (NA)")
    if not math.isfinite(vg_v):
        raise InvalidParameterError(f"Gate voltage VG must be finite: {vg_v}")

    cox = calculate_oxide_capacitance(tox_m)
    phi_f = calculate_fermi_potential(na_cm3, temp_k=temp_k)
    vfb = calculate_flat_band_voltage(phims_v=phims_v, qox_c_m2=qox_c_m2, cox_f_m2=cox)
    threshold_voltage = calculate_level1_threshold_voltage(
        vfb_v=vfb,
        phi_f_v=phi_f,
        na_cm3=na_cm3,
        cox_f_m2=cox,
    )
    eg_ev = silicon_bandgap_ev(temp_k)
    op_state = classify_mos_operating_state(vg_v, vfb, threshold_voltage)

    warnings = []
    solver_diagnostics: Dict[str, Any] = {}

    if level == 1:
        # Level 1 Analytical Depletion Approximation:
        # In accumulation (VG < VFB): psi_s ~ 0 (flat-band approx for bands, surface charge Qs > 0)
        # In depletion (VFB <= VG < VT):
        #   VG = VFB + psi_s + sqrt(2*eps_si*q*NA*psi_s) / Cox
        #   Let u = sqrt(psi_s). Then u^2 + (sqrt(2*eps*q*NA)/Cox)*u - (VG - VFB) = 0
        #   u = [ -gamma + sqrt(gamma^2 + 4*(VG - VFB)) ] / 2
        # In strong inversion (VG >= VT):
        #   psi_s is pinned at 2*phi_F
        gamma = math.sqrt(2.0 * EPS_SI * Q * cm3_to_m3(na_cm3)) / cox
        if vg_v < vfb:
            psi_s = 0.0
            qs = cox * (vfb - vg_v)  # Qs > 0 in accumulation
            wd = 0.0
        elif vg_v >= threshold_voltage:
            psi_s = 2.0 * phi_f
            wd = calculate_depletion_width(psi_s, na_cm3)
            # Total charge in inversion includes depletion + inversion electrons
            qs = -cox * (vg_v - vfb - psi_s)  # Qs < 0
        else:
            # Depletion regime: solve quadratic for sqrt(psi_s)
            v_dep = vg_v - vfb
            sqrt_psi = (-gamma + math.sqrt(gamma**2 + 4.0 * v_dep)) / 2.0
            psi_s = sqrt_psi**2
            wd = calculate_depletion_width(psi_s, na_cm3)
            qs = - math.sqrt(2.0 * EPS_SI * Q * cm3_to_m3(na_cm3) * psi_s)

    elif level == 2:
        # Level 2 Numerical Root Solving
        psi_s, solver_diagnostics = solve_surface_potential_level2(
            vg_v=vg_v,
            vfb_v=vfb,
            cox_f_m2=cox,
            na_cm3=na_cm3,
            phi_f_v=phi_f,
            temp_k=temp_k,
        )
        qs = calculate_semiconductor_charge_level2(
            psi_s_v=psi_s,
            na_cm3=na_cm3,
            phi_f_v=phi_f,
            temp_k=temp_k,
        )
        wd = calculate_depletion_width(psi_s, na_cm3) if psi_s > 0 else 0.0
    else:
        raise InvalidParameterError(f"MOS simulation level must be 1 or 2, got: {level}")

    # Spatial domain in silicon: x = 0 to x = x_bulk
    # Wd_max is at 2*phi_F
    wd_max = calculate_depletion_width(2.0 * phi_f, na_cm3)
    domain_depth = si_thickness_m if si_thickness_m is not None else max(3.0 * wd_max, 1.0e-6)
    if domain_depth <= 0.0:
        raise InvalidParameterError(f"Silicon thickness must be positive, got: {domain_depth}")

    # Spatial grid in silicon: x in [0, domain_depth]
    # More points near surface where bending occurs
    if wd > 0:
        x_dep = np.linspace(0.0, wd, max(num_points // 2, 20), endpoint=False)
        x_bulk = np.linspace(wd, domain_depth, max(num_points // 2, 20))
        x = np.concatenate([x_dep, x_bulk])
    else:
        x = np.linspace(0.0, domain_depth, num_points)
    x = np.unique(x)

    # Electrostatic potential profile phi(x) in silicon:
    # phi_bulk = 0 V
    # phi(0) = psi_s
    # Under the depletion approximation:
    # phi(x) = psi_s * (1 - x / Wd)^2 for 0 <= x <= Wd, and 0 for x > Wd
    potential = np.zeros_like(x)
    electric_field = np.zeros_like(x)
    if wd > 0:
        dep_mask = x <= wd
        potential[dep_mask] = psi_s * (1.0 - x[dep_mask] / wd)**2
        # E = -dphi/dx = 2*psi_s/Wd * (1 - x/Wd)
        electric_field[dep_mask] = (2.0 * psi_s / wd) * (1.0 - x[dep_mask] / wd)

    # Energy bands in Silicon:
    # EC(x) = EC_ref - phi(x)
    # EV(x) = EC(x) - Eg
    # Ei(x) = EC(x) - Eg/2
    # In bulk (phi = 0), p-type EF is flat at 0 eV (arbitrary ground reference).
    # Since bulk is p-type: EF_bulk = Ei_bulk - phi_F = 0 => Ei_bulk = phi_F
    # Since Ei = EC - Eg/2, EC_ref = phi_F + 0.5 * eg_ev
    ec_ref_ev = phi_f + 0.5 * eg_ev
    ec, ev, ei = bands_from_potential(potential, eg_ev=eg_ev, ec_ref_ev=ec_ref_ev)

    # Flat Fermi level EF = 0.0 eV in thermal equilibrium MOS structure
    ef = np.zeros_like(x)

    # Region classification array
    region = np.empty(x.shape, dtype=object)
    if wd > 0:
        region[x <= wd] = "semiconductor-depletion"
        region[x > wd] = "semiconductor-bulk"
    else:
        region[:] = "semiconductor-bulk"

    calculated_parameters = {
        "Cox_F_per_m2": cox,
        "phi_F_V": phi_f,
        "VFB_V": vfb,
        "threshold_voltage_V": threshold_voltage,
        "psi_s_V": psi_s,
        "Qs_C_per_m2": qs,
        "Wd_m": wd,
        "Wd_max_m": wd_max,
        "VG_V": vg_v,
        "operating_state": op_state,
        "temperature_K": temp_k,
        "NA_cm3": na_cm3,
        "tox_m": tox_m,
        "bandgap_eV": eg_ev,
        "level": level,
    }
    if solver_diagnostics:
        calculated_parameters["solver_diagnostics"] = solver_diagnostics

    assumptions = [
        "1-D vertical MOS electrostatics model",
        "p-type silicon substrate with zero bulk potential reference (phi_bulk = 0 V)",
        "Semiconductor energy bands plotted strictly in silicon (x >= 0, no bands in oxide)",
        "Zero gate leakage / ideal oxide dielectric",
        "Room temperature Boltzmann carrier statistics in silicon bulk",
    ]

    return SimulationResult(
        device="NMOS",
        model=f"1-D MOS electrostatic model (Level {level})",
        position=x,
        EC=ec,
        EV=ev,
        Ei=ei,
        EF=ef,
        EFn=None,
        EFp=None,
        potential=potential,
        electric_field=electric_field,
        region=region,
        calculated_parameters=calculated_parameters,
        operating_condition=f"{op_state} (VG = {vg_v:.3f} V, VFB = {vfb:.3f} V, VT = {threshold_voltage:.3f} V)",
        warnings=warnings,
        assumptions=assumptions,
        model_valid=True,
    )
