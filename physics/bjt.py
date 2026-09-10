"""NPN BJT two-junction analytical electrostatic model.

Structure:
    N+ emitter | P base | N collector

Voltage Definitions:
    VBE = VB - VE
    VBC = VB - VC

Junction Barrier Voltages:
    Vbarrier_EB = Vbi_EB - VBE
    Vbarrier_CB = Vbi_CB - VBC

Where built-in potentials are:
    Vbi_EB = V_T * ln(NE * NB / ni^2)
    Vbi_CB = V_T * ln(NC * NB / ni^2)

Depletion Widths (abrupt approximation):
    W_EB = sqrt( (2 * eps_si / q) * (1/NE + 1/NB) * Vbarrier_EB )
    x_e = (NB / (NE + NB)) * W_EB
    x_be = (NE / (NE + NB)) * W_EB

    W_CB = sqrt( (2 * eps_si / q) * (1/NC + 1/NB) * Vbarrier_CB )
    x_c = (NB / (NC + NB)) * W_CB
    x_bc = (NC / (NC + NB)) * W_CB

Operating Regions:
    Tolerance: tol = 1e-3 V (1 mV) around 0 V for near-zero transition.
    If |VBE| <= tol or |VBC| <= tol:
        Operating condition is classified as "Transition / near-zero bias" (or equilibrium if both == 0).
    Otherwise:
        EB reverse (VBE < -tol) + BC reverse (VBC < -tol) -> "Cutoff"
        EB forward (VBE > tol)  + BC reverse (VBC < -tol) -> "Forward active"
        EB forward (VBE > tol)  + BC forward (VBC > tol)  -> "Saturation"
        EB reverse (VBE < -tol) + BC forward (VBC > tol)  -> "Reverse active"

Base Depletion Overlap:
    If x_be + x_bc >= W_base:
        Depletion regions overlap (punch-through condition).
        Raises ModelValidityError (or invalid flag according to architecture; never silently clipped).

Energy Bands & Fermi Levels:
    At equilibrium (VBE == 0, VBC == 0):
        Flat EF across the entire device.
    Under bias:
        No global EF.
        Quasi-Fermi levels defined only in neutral regions:
            EFn in neutral emitter
            EFp in neutral base
            EFn in neutral collector
        No fabricated quasi-Fermi profiles through depletion or base transport regions.
        EC - EV = Eg strictly maintained.
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
from physics.pn_junction import calculate_pn_built_in_potential, calculate_pn_depletion_widths
from physics.results import SimulationResult


def classify_bjt_operating_region(vbe_v: float, vbc_v: float, tol_v: float = 1e-3) -> str:
    """Classify NPN BJT operating region based on junction voltages VBE and VBC.

    Args:
        vbe_v: Base-emitter voltage VB - VE [V].
        vbc_v: Base-collector voltage VB - VC [V].
        tol_v: Voltage deadband threshold around zero for transition state [V].

    Returns:
        Region string: 'Equilibrium', 'Transition / near-zero bias', 'Cutoff',
        'Forward active', 'Saturation', or 'Reverse active'.
    """
    if not (math.isfinite(vbe_v) and math.isfinite(vbc_v)):
        raise InvalidParameterError(f"Junction voltages must be finite: VBE={vbe_v}, VBC={vbc_v}")

    if abs(vbe_v) < 1e-9 and abs(vbc_v) < 1e-9:
        return "Equilibrium"

    if abs(vbe_v) <= tol_v or abs(vbc_v) <= tol_v:
        return "Transition / near-zero bias"

    eb_forward = vbe_v > tol_v
    bc_forward = vbc_v > tol_v

    if not eb_forward and not bc_forward:
        return "Cutoff"
    elif eb_forward and not bc_forward:
        return "Forward active"
    elif eb_forward and bc_forward:
        return "Saturation"
    else:
        return "Reverse active"


def simulate_npn_bjt(
    ne_cm3: float = 1.0e18,
    nb_cm3: float = 1.0e16,
    nc_cm3: float = 1.0e15,
    vbe_v: float = 0.0,
    vbc_v: float = 0.0,
    temp_k: float = 300.0,
    w_base_m: float = 1.0e-6,
    w_emitter_m: Optional[float] = None,
    w_collector_m: Optional[float] = None,
    num_points: int = 600,
) -> SimulationResult:
    """Simulate 1-D NPN BJT using coupled two-junction analytical electrostatic model.

    Spatial structure:
        N+ emitter (x < 0)
        EB metallurgical junction at x = 0
        P base (0 <= x <= W_base)
        BC metallurgical junction at x = W_base
        N collector (x > W_base)

    Args:
        ne_cm3: Emitter donor doping concentration in cm^-3 (N+).
        nb_cm3: Base acceptor doping concentration in cm^-3 (P).
        nc_cm3: Collector donor doping concentration in cm^-3 (N).
        vbe_v: Base-emitter voltage VB - VE in Volts [V].
        vbc_v: Base-collector voltage VB - VC in Volts [V].
        temp_k: Temperature in Kelvin [K].
        w_base_m: Metallurgical base width in meters [m].
        w_emitter_m: Neutral emitter thickness in meters [m]. Defaults to 2 * x_e.
        w_collector_m: Neutral collector thickness in meters [m]. Defaults to 2 * x_c.
        num_points: Total number of spatial grid points.

    Returns:
        SimulationResult containing BJT band diagram and electrostatics.

    Raises:
        InvalidParameterError: For invalid doping, temperature, dimensions.
        ModelValidityError: For barrier <= 0 or base punch-through (depletion overlap).
    """
    validate_temperature(temp_k)
    validate_doping(ne_cm3, doping_name="Emitter doping (NE)")
    validate_doping(nb_cm3, doping_name="Base doping (NB)")
    validate_doping(nc_cm3, doping_name="Collector doping (NC)")

    if not (math.isfinite(vbe_v) and math.isfinite(vbc_v)):
        raise InvalidParameterError(f"Voltages must be finite: VBE={vbe_v}, VBC={vbc_v}")

    if not math.isfinite(w_base_m) or w_base_m <= 0.0:
        raise InvalidParameterError(f"Base width must be positive finite: {w_base_m}")

    eg_ev = silicon_bandgap_ev(temp_k)
    vt_v = thermal_voltage(temp_k)
    ni_cm3 = intrinsic_carrier_concentration(temp_k)

    # Junction built-in potentials
    vbi_eb = calculate_pn_built_in_potential(na_cm3=nb_cm3, nd_cm3=ne_cm3, temp_k=temp_k)
    vbi_cb = calculate_pn_built_in_potential(na_cm3=nb_cm3, nd_cm3=nc_cm3, temp_k=temp_k)

    # Barrier voltages: Vbarrier = Vbi - V_applied_junction
    # For EB junction: P-side is base, N-side is emitter => VA_EB = VB - VE = VBE
    # For CB junction: P-side is base, N-side is collector => VA_CB = VB - VC = VBC
    vbarrier_eb = vbi_eb - vbe_v
    vbarrier_cb = vbi_cb - vbc_v

    # Model validity checks: barriers must be > 0 for abrupt depletion approximation
    if vbarrier_eb <= 0.0:
        raise ModelValidityError(
            f"EB junction barrier Vbarrier_EB = {vbarrier_eb:.4e} V <= 0. "
            "Abrupt depletion approximation is invalid under heavy forward bias."
        )
    if vbarrier_cb <= 0.0:
        raise ModelValidityError(
            f"CB junction barrier Vbarrier_CB = {vbarrier_cb:.4e} V <= 0. "
            "Abrupt depletion approximation is invalid under heavy forward bias."
        )

    # Calculate depletion widths
    w_eb, x_be, x_e = calculate_pn_depletion_widths(na_cm3=nb_cm3, nd_cm3=ne_cm3, vd_v=vbarrier_eb)
    w_cb, x_bc, x_c = calculate_pn_depletion_widths(na_cm3=nb_cm3, nd_cm3=nc_cm3, vd_v=vbarrier_cb)

    # Base depletion overlap detection (punch-through)
    total_base_depletion = x_be + x_bc
    if total_base_depletion >= w_base_m:
        raise ModelValidityError(
            f"Base depletion overlap detected: x_be ({x_be*1e9:.1f} nm) + x_bc ({x_bc*1e9:.1f} nm) = "
            f"{total_base_depletion*1e9:.1f} nm >= W_base ({w_base_m*1e9:.1f} nm). "
            "Independent-junction electrostatic approximation is invalid (punch-through)."
        )

    # Region dimensions
    le_neutral = w_emitter_m if w_emitter_m is not None else 2.5 * x_e
    lc_neutral = w_collector_m if w_collector_m is not None else 2.5 * x_c
    if le_neutral <= 0.0 or lc_neutral <= 0.0:
        raise InvalidParameterError("Emitter and collector neutral thicknesses must be positive.")

    # Spatial coordinates:
    # x = 0 is EB metallurgical junction
    # Emitter depletion: [-x_e, 0]
    # Base EB depletion: [0, x_be]
    # Base neutral: [x_be, w_base_m - x_bc]
    # Base CB depletion: [w_base_m - x_bc, w_base_m]
    # Collector depletion: [w_base_m, w_base_m + x_c]
    x_min = - (le_neutral + x_e)
    x_max = w_base_m + x_c + lc_neutral

    pts = max(num_points // 6, 20)
    g_e_neutral = np.linspace(x_min, -x_e, pts, endpoint=False)
    g_eb_dep_n = np.linspace(-x_e, 0.0, pts, endpoint=False)
    g_eb_dep_p = np.linspace(0.0, x_be, pts, endpoint=False)
    g_b_neutral = np.linspace(x_be, w_base_m - x_bc, pts, endpoint=False)
    g_cb_dep_p = np.linspace(w_base_m - x_bc, w_base_m, pts, endpoint=False)
    g_cb_dep_n = np.linspace(w_base_m, w_base_m + x_c, pts, endpoint=True)
    g_c_neutral = np.linspace(w_base_m + x_c, x_max, pts + 1)[1:]

    x = np.concatenate([
        g_e_neutral, g_eb_dep_n, g_eb_dep_p, g_b_neutral, g_cb_dep_p, g_cb_dep_n, g_c_neutral
    ])
    x = np.unique(x)

    # Electrostatic potential phi(x)
    # Define reference: phi = 0 V in neutral base (x_be <= x <= w_base_m - x_bc)
    # At EB junction (x <= x_be):
    # Potential drops by Vbarrier_EB from base (p) down to emitter (n):
    # In base EB depletion (0 <= x <= x_be):
    #   phi(x) = - (q * NB / (2 * eps_si)) * (x_be - x)^2
    # In emitter EB depletion (-x_e <= x <= 0):
    #   phi(x) = - Vbarrier_EB + (q * NE / (2 * eps_si)) * (x + x_e)^2
    # In neutral emitter (x < -x_e):
    #   phi(x) = - Vbarrier_EB
    #
    # At CB junction (x >= w_base_m - x_bc):
    # Potential drops by Vbarrier_CB from base (p) down to collector (n):
    # In base CB depletion (w_base_m - x_bc <= x <= w_base_m):
    #   phi(x) = - (q * NB / (2 * eps_si)) * (x - (w_base_m - x_bc))^2
    # In collector CB depletion (w_base_m <= x <= w_base_m + x_c):
    #   phi(x) = - Vbarrier_CB + (q * NC / (2 * eps_si)) * (w_base_m + x_c - x)^2
    # In neutral collector (x > w_base_m + x_c):
    #   phi(x) = - Vbarrier_CB

    ne_m3 = cm3_to_m3(ne_cm3)
    nb_m3 = cm3_to_m3(nb_cm3)
    nc_m3 = cm3_to_m3(nc_cm3)

    potential = np.zeros_like(x)
    electric_field = np.zeros_like(x)
    region = np.empty(x.shape, dtype=object)

    m_e_neut = x < -x_e
    m_eb_n = (x >= -x_e) & (x < 0.0)
    m_eb_p = (x >= 0.0) & (x < x_be)
    m_b_neut = (x >= x_be) & (x <= w_base_m - x_bc)
    m_cb_p = (x > w_base_m - x_bc) & (x <= w_base_m)
    m_cb_n = (x > w_base_m) & (x <= w_base_m + x_c)
    m_c_neut = x > w_base_m + x_c

    # Potential
    potential[m_e_neut] = -vbarrier_eb
    potential[m_eb_n] = -vbarrier_eb + (Q * ne_m3 / (2.0 * EPS_SI)) * (x[m_eb_n] + x_e)**2
    potential[m_eb_p] = - (Q * nb_m3 / (2.0 * EPS_SI)) * (x_be - x[m_eb_p])**2
    potential[m_b_neut] = 0.0
    potential[m_cb_p] = - (Q * nb_m3 / (2.0 * EPS_SI)) * (x[m_cb_p] - (w_base_m - x_bc))**2
    potential[m_cb_n] = -vbarrier_cb + (Q * nc_m3 / (2.0 * EPS_SI)) * ((w_base_m + x_c) - x[m_cb_n])**2
    potential[m_c_neut] = -vbarrier_cb

    # Electric field E = -dphi/dx
    # E(-x_e) = 0, E(0) = - q*NE*x_e/eps, E(x_be) = 0
    # E points in -x direction across EB junction (from n to p? No: E points from + to -, so + donor in n to - acceptor in p, which is +x direction!)
    # Let's verify carefully:
    # In EB: N+ on left (x < 0), P on right (x > 0).
    # Donors on left (+ charge), acceptors on right (- charge).
    # Field points from + to -, so E points in +x direction!
    # dphi/dx: phi is -Vbarrier on left, 0 in base. So dphi/dx > 0, which means -dphi/dx < 0 if phi increases?
    # Wait: phi(-x_e) = -Vbarrier_EB, phi(x_be) = 0. Potential INCREASES from N+ emitter to P base.
    # Therefore dphi/dx > 0, so E = -dphi/dx < 0.
    # Wait, in electrostatics:
    # Poisson: d^2 phi / dx^2 = - rho / eps.
    # In N+ region (x < 0): rho = +q*ND > 0.
    # d^2 phi / dx^2 = - (positive) < 0, so phi is concave down.
    # But earlier we had phi(p) = 0 and phi(n) = VD in PN junction!
    # In standard PN junction with P on left and N on right:
    # Left is P (negative space charge), right is N (positive space charge).
    # Field points from + (right) to - (left): points in -x direction.
    # Potential at left (P) was 0, potential at right (N) was +VD.
    # In NPN BJT:
    # Emitter is N+ on LEFT (x < 0). Base is P in MIDDLE. Collector is N on RIGHT.
    # Left (Emitter) has + space charge (donors).
    # Middle (Base) has - space charge (acceptors).
    # Right (Collector) has + space charge (donors).
    # Field across EB junction: points from Emitter (+) to Base (-), so points in +x direction!
    # Potential: E = -dphi/dx => dphi/dx = -E < 0 => phi decreases from Emitter to Base!
    # Let's verify: In N+ silicon, conduction band is CLOSE to Fermi level.
    # Since EC = EC_ref - phi, a lower EC in N+ means phi is HIGHER in N+ than in P!
    # In N region: EC is lower than in P region.
    # Since EC = EC_ref - phi, phi_N > phi_P!
    # Yes! Electrostatic potential of N-type silicon is HIGHER than P-type silicon!
    # Let's re-verify:
    # In PN junction (P left, N right): phi_P = 0, phi_N = +VD > 0. (N is HIGHER than P).
    # In NPN BJT: Emitter (N) is left, Base (P) is middle, Collector (N) is right.
    # Therefore, phi_Emitter > phi_Base, and phi_Collector > phi_Base!
    # In the neutral base (P), let phi_Base = 0 V.
    # Then phi_Emitter = + Vbarrier_EB!
    # And phi_Collector = + Vbarrier_CB!
    # Let's check EC:
    # In Base (P): phi = 0 => EC_base = EC_ref.
    # In Emitter (N): phi = +Vbarrier_EB => EC_emitter = EC_ref - Vbarrier_EB (LOWER than Base, as it must be for N+!).
    # In Collector (N): phi = +Vbarrier_CB => EC_collector = EC_ref - Vbarrier_CB (LOWER than Base, as it must be for N!).
    # This is physically and mathematically consistent!

    potential[m_e_neut] = vbarrier_eb
    potential[m_eb_n] = vbarrier_eb - (Q * ne_m3 / (2.0 * EPS_SI)) * (x[m_eb_n] + x_e)**2
    potential[m_eb_p] = (Q * nb_m3 / (2.0 * EPS_SI)) * (x_be - x[m_eb_p])**2
    potential[m_b_neut] = 0.0
    potential[m_cb_p] = (Q * nb_m3 / (2.0 * EPS_SI)) * (x[m_cb_p] - (w_base_m - x_bc))**2
    potential[m_cb_n] = vbarrier_cb - (Q * nc_m3 / (2.0 * EPS_SI)) * ((w_base_m + x_c) - x[m_cb_n])**2
    potential[m_c_neut] = vbarrier_cb

    # Electric field E = -dphi/dx:
    # For EB junction:
    # In EB emitter depletion (-x_e <= x <= 0):
    # phi = vbarrier_eb - (Q*ne / 2eps) * (x + x_e)^2
    # dphi/dx = - (Q*ne / eps) * (x + x_e) <= 0
    # E = - dphi/dx = + (Q*ne / eps) * (x + x_e) >= 0 (points +x, from N+ to P)
    # In EB base depletion (0 <= x <= x_be):
    # phi = (Q*nb / 2eps) * (x_be - x)^2
    # dphi/dx = - (Q*nb / eps) * (x_be - x) <= 0
    # E = - dphi/dx = + (Q*nb / eps) * (x_be - x) >= 0
    # Peak at x = 0: E = + Q*ne*x_e / eps = + Q*nb*x_be / eps.
    electric_field[m_eb_n] = (Q * ne_m3 / EPS_SI) * (x[m_eb_n] + x_e)
    electric_field[m_eb_p] = (Q * nb_m3 / EPS_SI) * (x_be - x[m_eb_p])

    # For CB junction:
    # Base is left (P, - charge), Collector is right (N, + charge).
    # Field points from + (Collector) to - (Base), which is -x direction!
    # In CB base depletion (w_base_m - x_bc <= x <= w_base_m):
    # phi = (Q*nb / 2eps) * (x - (w_base_m - x_bc))^2
    # dphi/dx = + (Q*nb / eps) * (x - (w_base_m - x_bc)) >= 0
    # E = - dphi/dx = - (Q*nb / eps) * (x - (w_base_m - x_bc)) <= 0
    # In CB collector depletion (w_base_m <= x <= w_base_m + x_c):
    # phi = vbarrier_cb - (Q*nc / 2eps) * ((w_base_m + x_c) - x)^2
    # dphi/dx = + (Q*nc / eps) * ((w_base_m + x_c) - x) >= 0
    # E = - dphi/dx = - (Q*nc / eps) * ((w_base_m + x_c) - x) <= 0
    electric_field[m_cb_p] = -(Q * nb_m3 / EPS_SI) * (x[m_cb_p] - (w_base_m - x_bc))
    electric_field[m_cb_n] = -(Q * nc_m3 / EPS_SI) * ((w_base_m + x_c) - x[m_cb_n])

    region[m_e_neut] = "emitter-neutral"
    region[m_eb_n] = "eb-depletion-emitter"
    region[m_eb_p] = "eb-depletion-base"
    region[m_b_neut] = "base-neutral"
    region[m_cb_p] = "cb-depletion-base"
    region[m_cb_n] = "cb-depletion-collector"
    region[m_c_neut] = "collector-neutral"

    # Energy bands [eV]
    # EC = EC_ref - phi
    # EV = EC - Eg
    # Ei = EC - Eg / 2
    # Reference: In base (phi = 0), EFp is at delta_ef_b below Ei.
    # Let Base EFp = 0 eV (arbitrary ground reference).
    # Then Ei_base = delta_ef_b.
    # EC_base = Ei_base + Eg/2 = delta_ef_b + Eg/2.
    # Since phi_base = 0, EC_ref = delta_ef_b + 0.5 * eg_ev.
    delta_ef_b = (K_B_EV * temp_k) * math.log(nb_cm3 / ni_cm3)
    ec_ref_ev = delta_ef_b + 0.5 * eg_ev
    ec, ev, ei = bands_from_potential(potential, eg_ev=eg_ev, ec_ref_ev=ec_ref_ev)

    is_equilibrium = (abs(vbe_v) < 1e-9) and (abs(vbc_v) < 1e-9)
    op_region = classify_bjt_operating_region(vbe_v, vbc_v)

    ef: Optional[np.ndarray] = None
    efn: Optional[np.ndarray] = None
    efp: Optional[np.ndarray] = None

    if is_equilibrium:
        # Single flat EF = 0.0 eV throughout
        ef = np.zeros_like(x)
    else:
        # Under bias: do not invent a global EF or spatial transport quasi-Fermi curve
        ef = None
        efn = np.full_like(x, fill_value=np.nan)
        efp = np.full_like(x, fill_value=np.nan)

        # In neutral base: EFp = 0.0 eV
        efp[m_b_neut] = 0.0

        # Since VBE = VB - VE, VE = VB - VBE = -VBE relative to base.
        # In neutral emitter: EFn = -VBE
        efn[m_e_neut] = -vbe_v

        # Since VBC = VB - VC, VC = VB - VBC = -VBC relative to base.
        # In neutral collector: EFn = -VBC
        efn[m_c_neut] = -vbc_v

    calculated_parameters = {
        "VBE_V": vbe_v,
        "VBC_V": vbc_v,
        "Vbi_EB_V": vbi_eb,
        "Vbi_CB_V": vbi_cb,
        "Vbarrier_EB_V": vbarrier_eb,
        "Vbarrier_CB_V": vbarrier_cb,
        "W_EB_m": w_eb,
        "W_CB_m": w_cb,
        "x_e_m": x_e,
        "x_be_m": x_be,
        "x_bc_m": x_bc,
        "x_c_m": x_c,
        "W_base_m": w_base_m,
        "neutral_base_width_m": float(w_base_m - total_base_depletion),
        "total_base_depletion_m": float(total_base_depletion),
        "operating_region": op_region,
        "temperature_K": temp_k,
        "bandgap_eV": eg_ev,
        "NE_cm3": ne_cm3,
        "NB_cm3": nb_cm3,
        "NC_cm3": nc_cm3,
    }

    assumptions = [
        "1-D NPN abrupt two-junction analytical electrostatic model",
        "Independent abrupt depletion approximations for EB and CB junctions",
        "Uniform doping within emitter (N+), base (P), and collector (N)",
        "Zero electric field in neutral regions",
        "No carrier transport, recombination-generation, or current/gain calculations",
        "Quasi-Fermi levels shown only as reference values in neutral regions",
    ]

    warnings = []
    if (w_base_m - total_base_depletion) < 0.2 * w_base_m:
        warnings.append(
            "Base width is heavily depleted (>80% of metallurgical base width). "
            "Base-width modulation (Early effect) is severe."
        )

    return SimulationResult(
        device="NPN BJT",
        model="Two-junction analytical electrostatic approximation",
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
        operating_condition=f"{op_region} (VBE = {vbe_v:.3f} V, VBC = {vbc_v:.3f} V)",
        warnings=warnings,
        assumptions=assumptions,
        model_valid=True,
    )
