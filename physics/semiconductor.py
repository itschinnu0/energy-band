"""Shared semiconductor physics calculations.

Reference: S. M. Sze and Kwok K. Ng, 'Physics of Semiconductor Devices', 3rd ed.
Calculations follow Section 7 of IMPLEMENTATION_PLAN_SOURCE.md.
"""

from typing import Union, Tuple, Optional
import math
import numpy as np

from physics.constants import (
    Q,
    K_B,
    K_B_EV,
    EG_0_EV,
    VARSHNI_ALPHA_EV,
    VARSHNI_BETA_K,
    NC_300_CM3,
    NV_300_CM3,
    NC_300_M3,
    NV_300_M3,
    T_REF_K,
)
from physics.units import cm3_to_m3, m3_to_cm3, ev_to_joules, joules_to_ev
from physics.validation import validate_temperature, validate_doping
from physics.results import SemiconductorProperties


def thermal_voltage(temp_k: float) -> float:
    """Calculate thermal voltage V_T = k_B * T / q in Volts [V]."""
    validate_temperature(temp_k)
    return (K_B * temp_k) / Q


def silicon_bandgap_ev(temp_k: float) -> float:
    """Calculate silicon bandgap Eg(T) in eV using Varshni equation.

    Eg(T) = Eg0 - alpha * T^2 / (T + beta)
    Eg0 = 1.166 eV, alpha = 4.73e-4 eV/K, beta = 636 K.
    At 300 K, this yields approximately 1.125 eV.
    """
    validate_temperature(temp_k)
    return EG_0_EV - (VARSHNI_ALPHA_EV * temp_k**2) / (temp_k + VARSHNI_BETA_K)


def silicon_bandgap_joules(temp_k: float) -> float:
    """Calculate silicon bandgap in Joules [J]."""
    return ev_to_joules(silicon_bandgap_ev(temp_k))


def effective_density_of_states(temp_k: float) -> Tuple[float, float]:
    """Calculate effective density of states (Nc, Nv) in cm^-3 at temperature temp_k.

    Nc(T) = Nc(300) * (T / 300)^(3/2)
    Nv(T) = Nv(300) * (T / 300)^(3/2)
    """
    validate_temperature(temp_k)
    factor = (temp_k / T_REF_K) ** 1.5
    nc = NC_300_CM3 * factor
    nv = NV_300_CM3 * factor
    return nc, nv


def intrinsic_carrier_concentration(temp_k: float) -> float:
    """Calculate intrinsic carrier concentration n_i(T) in cm^-3.

    n_i(T) = sqrt(Nc(T) * Nv(T)) * exp(-Eg(T) / (2 * k_B * T))
    """
    validate_temperature(temp_k)
    nc, nv = effective_density_of_states(temp_k)
    eg = silicon_bandgap_ev(temp_k)
    vt_ev = K_B_EV * temp_k
    ni = math.sqrt(nc * nv) * math.exp(-eg / (2.0 * vt_ev))
    return ni


def get_semiconductor_properties(temp_k: float = 300.0) -> SemiconductorProperties:
    """Return common semiconductor properties for silicon at temperature temp_k."""
    validate_temperature(temp_k)
    eg_ev = silicon_bandgap_ev(temp_k)
    vt_v = thermal_voltage(temp_k)
    nc_cm3, nv_cm3 = effective_density_of_states(temp_k)
    ni_cm3 = intrinsic_carrier_concentration(temp_k)
    return SemiconductorProperties(
        temperature_k=temp_k,
        bandgap_ev=eg_ev,
        thermal_voltage_v=vt_v,
        nc_cm3=nc_cm3,
        nv_cm3=nv_cm3,
        ni_cm3=ni_cm3,
        nc_m3=cm3_to_m3(nc_cm3),
        nv_m3=cm3_to_m3(nv_cm3),
        ni_m3=cm3_to_m3(ni_cm3),
    )


def equilibrium_fermi_level(
    doping_cm3: float,
    doping_type: str,
    temp_k: float = 300.0,
    ei_ev: float = 0.0,
) -> float:
    """Calculate equilibrium Fermi level EF in eV referenced to Ei.

    For non-degenerate p-type:
        Ei - EF = kT * ln(NA / ni) => EF = Ei - kT * ln(NA / ni)
    For non-degenerate n-type:
        EF - Ei = kT * ln(ND / ni) => EF = Ei + kT * ln(ND / ni)
        
    Args:
        doping_cm3: Net acceptor NA or donor ND concentration in cm^-3.
        doping_type: 'p' (acceptor) or 'n' (donor).
        temp_k: Temperature in Kelvin.
        ei_ev: Intrinsic Fermi level energy in eV (default: 0.0 eV reference).
        
    Returns:
        EF in eV.
    """
    validate_temperature(temp_k)
    validate_doping(doping_cm3)
    doping_type_lower = doping_type.lower().strip()
    if doping_type_lower not in ('p', 'n'):
        raise ValueError(f"doping_type must be 'p' or 'n', got '{doping_type}'")

    ni = intrinsic_carrier_concentration(temp_k)
    vt_ev = K_B_EV * temp_k  # kT in eV

    # Fermi potential / separation from Ei in eV
    delta_ef = vt_ev * math.log(doping_cm3 / ni)

    if doping_type_lower == 'n':
        return ei_ev + delta_ef
    else:
        return ei_ev - delta_ef


def carrier_densities_from_quasi_fermi(
    efn_ev: Union[float, np.ndarray],
    efp_ev: Union[float, np.ndarray],
    ei_ev: Union[float, np.ndarray],
    temp_k: float = 300.0,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """Calculate electron and hole concentrations from quasi-Fermi levels under Boltzmann approximation.

    n = n_i * exp((EFn - Ei) / (k_B * T))
    p = n_i * exp((Ei - EFp) / (k_B * T))
    
    np = n_i^2 * exp((EFn - EFp) / (k_B * T))

    Args:
        efn_ev: Electron quasi-Fermi level in eV.
        efp_ev: Hole quasi-Fermi level in eV.
        ei_ev: Intrinsic Fermi level in eV.
        temp_k: Temperature in Kelvin.

    Returns:
        (n, p) concentrations in cm^-3.
    """
    validate_temperature(temp_k)
    ni = intrinsic_carrier_concentration(temp_k)
    vt_ev = K_B_EV * temp_k

    n = ni * np.exp((efn_ev - ei_ev) / vt_ev)
    p = ni * np.exp((ei_ev - efp_ev) / vt_ev)
    return n, p


def bands_from_potential(
    potential_v: np.ndarray,
    eg_ev: float,
    ec_ref_ev: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate EC, EV, and Ei band profiles in eV from electrostatic potential phi(x).

    Conventions (locked in Section 2.2 and Section 46):
        EC(x) = EC_ref - phi(x)
        EV(x) = EC(x) - Eg
        Ei(x) = (EC(x) + EV(x)) / 2  (mid-gap approximation for non-degenerate Boltzmann)
        EC(x) - EV(x) = Eg is strictly maintained.

    Args:
        potential_v: Array of electrostatic potential phi(x) in Volts.
        eg_ev: Bandgap in eV.
        ec_ref_ev: Reference energy for EC at phi = 0 in eV (default: 0.0).

    Returns:
        (EC_eV, EV_eV, Ei_eV) arrays in eV.
    """
    ec = ec_ref_ev - potential_v
    ev = ec - eg_ev
    ei = ec - (eg_ev / 2.0)
    return ec, ev, ei
