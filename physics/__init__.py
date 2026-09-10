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

from physics.mos import (
    calculate_oxide_capacitance,
    calculate_fermi_potential,
    calculate_flat_band_voltage,
    calculate_level1_threshold_voltage,
    calculate_depletion_width,
    classify_mos_operating_state,
    calculate_semiconductor_charge_level2,
    solve_surface_potential_level2,
    simulate_mos,
)
from physics.solver import solve_bracketed_root

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
    "calculate_oxide_capacitance",
    "calculate_fermi_potential",
    "calculate_flat_band_voltage",
    "calculate_level1_threshold_voltage",
    "calculate_depletion_width",
    "classify_mos_operating_state",
    "calculate_semiconductor_charge_level2",
    "solve_surface_potential_level2",
    "simulate_mos",
    "solve_bracketed_root",
]
