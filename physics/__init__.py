"""Physics package for Energy Band Diagram Simulator."""

from physics.exceptions import (
    EnergyBandError,
    InvalidParameterError,
    ModelValidityError,
    SolverConvergenceError,
    UnitConversionError,
)
from physics.results import SemiconductorProperties, SimulationResult
from physics.electrostatics import (
    electric_field_from_potential,
    potential_from_electric_field,
    electric_field_from_charge_density,
    potential_from_charge_density,
)

from physics.pn_junction import (
    calculate_pn_built_in_potential,
    calculate_pn_depletion_widths,
    simulate_pn_junction,
)
from physics.bjt import (
    classify_bjt_operating_region,
    simulate_npn_bjt,
)

__all__ = [
    "EnergyBandError",
    "InvalidParameterError",
    "ModelValidityError",
    "SolverConvergenceError",
    "UnitConversionError",
    "SemiconductorProperties",
    "SimulationResult",
    "electric_field_from_potential",
    "potential_from_electric_field",
    "electric_field_from_charge_density",
    "potential_from_charge_density",
    "calculate_pn_built_in_potential",
    "calculate_pn_depletion_widths",
    "simulate_pn_junction",
    "classify_bjt_operating_region",
    "simulate_npn_bjt",
]
