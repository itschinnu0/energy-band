"""Common result object definitions and contracts."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
import numpy as np


@dataclass
class SemiconductorProperties:
    """Equilibrium and material properties of a semiconductor region."""
    temperature_k: float
    bandgap_ev: float
    thermal_voltage_v: float
    nc_cm3: float
    nv_cm3: float
    ni_cm3: float
    nc_m3: float
    nv_m3: float
    ni_m3: float


@dataclass
class SimulationResult:
    """Common result contract for device simulation.

    Follows the specification in AGENTS.md and Section 48 of IMPLEMENTATION_PLAN_SOURCE.md.
    Quantities not applicable or unavailable should be None. Never use fabricated zero curves.
    """
    device: str
    model: str
    position: np.ndarray  # [m]
    EC: np.ndarray  # [eV]
    EV: np.ndarray  # [eV]
    Ei: np.ndarray  # [eV]
    EF: Optional[np.ndarray] = None  # [eV] (flat in equilibrium; None if biased)
    EFn: Optional[np.ndarray] = None  # [eV] (electron quasi-Fermi level)
    EFp: Optional[np.ndarray] = None  # [eV] (hole quasi-Fermi level)
    potential: Optional[np.ndarray] = None  # [V]
    electric_field: Optional[np.ndarray] = None  # [V/m]
    region: Optional[np.ndarray] = None  # Region labels (e.g., 'p', 'n', 'depletion', 'oxide', 'bulk')
    calculated_parameters: Dict[str, Any] = field(default_factory=dict)
    operating_condition: str = ""
    warnings: List[str] = field(default_factory=list)
    assumptions: List[str] = field(default_factory=list)
    model_valid: bool = True
