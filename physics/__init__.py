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
]
