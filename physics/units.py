"""Unit conversions and handling.

Explicit conversions between engineering units and SI units.
Maintains clear distinction between:
- potential: Volts [V]
- energy: Joules [J]
- energy: electron-volts [eV]
"""

from typing import Union
import numpy as np
from physics.constants import Q


def ev_to_joules(energy_ev: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert energy from electron-volts (eV) to Joules (J)."""
    return energy_ev * Q


def joules_to_ev(energy_j: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert energy from Joules (J) to electron-volts (eV)."""
    return energy_j / Q


def potential_to_energy_ev(
    potential_v: Union[float, np.ndarray],
    ref_energy_ev: float = 0.0,
    charge_multiplier: float = -1.0,
) -> Union[float, np.ndarray]:
    """Convert electrostatic potential (V) to electron energy (eV).

    By definition for electrons (charge q_e = -q):
    E(x) = E_ref - q * phi(x)
    In eV units where 1 eV = q * (1 V):
    E_eV(x) = E_ref_eV - phi(x)
    
    Here charge_multiplier is -1.0 for electrons.
    """
    return ref_energy_ev + charge_multiplier * potential_v


def cm3_to_m3(conc_cm3: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert concentration from cm^-3 to m^-3."""
    return conc_cm3 * 1.0e6


def m3_to_cm3(conc_m3: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert concentration from m^-3 to cm^-3."""
    return conc_m3 * 1.0e-6


def um_to_m(length_um: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from micrometers (um) to meters (m)."""
    return length_um * 1.0e-6


def m_to_um(length_m: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from meters (m) to micrometers (um)."""
    return length_m * 1.0e6


def nm_to_m(length_nm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from nanometers (nm) to meters (m)."""
    return length_nm * 1.0e-9


def m_to_nm(length_m: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from meters (m) to nanometers (nm)."""
    return length_m * 1.0e9


def cm_to_m(length_cm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from centimeters (cm) to meters (m)."""
    return length_cm * 1.0e-2


def m_to_cm(length_m: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert length from meters (m) to centimeters (cm)."""
    return length_m * 1.0e2
