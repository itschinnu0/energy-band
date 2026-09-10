"""Tests for result dataclass structure and invariants."""

import pytest
import numpy as np
from physics.results import SimulationResult, SemiconductorProperties


def test_simulation_result_contract():
    """Verify SimulationResult fields and defaults."""
    x = np.linspace(0, 1e-6, 100)
    ec = np.zeros(100)
    ev = ec - 1.12
    ei = ec - 0.56

    result = SimulationResult(
        device="PN Junction",
        model="Abrupt 1-D depletion approximation",
        position=x,
        EC=ec,
        EV=ev,
        Ei=ei,
    )
    assert result.device == "PN Junction"
    assert result.model == "Abrupt 1-D depletion approximation"
    assert len(result.position) == 100
    assert result.EF is None
    assert result.EFn is None
    assert result.EFp is None
    assert result.model_valid is True
    assert isinstance(result.warnings, list)
    assert isinstance(result.assumptions, list)


def test_semiconductor_properties_contract():
    """Verify SemiconductorProperties fields."""
    props = SemiconductorProperties(
        temperature_k=300.0,
        bandgap_ev=1.125,
        thermal_voltage_v=0.02586,
        nc_cm3=2.86e19,
        nv_cm3=2.66e19,
        ni_cm3=1.0e10,
        nc_m3=2.86e25,
        nv_m3=2.66e25,
        ni_m3=1.0e16,
    )
    assert props.temperature_k == 300.0
    assert props.bandgap_ev == 1.125
    assert props.nc_m3 == 2.86e25

