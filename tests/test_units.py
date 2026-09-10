"""Tests for physical unit conversions and energy/potential relations."""

import pytest
import numpy as np
from physics.constants import Q
from physics.units import (
    ev_to_joules,
    joules_to_ev,
    potential_to_energy_ev,
    cm3_to_m3,
    m3_to_cm3,
    um_to_m,
    m_to_um,
    nm_to_m,
    m_to_nm,
    cm_to_m,
    m_to_cm,
)


def test_ev_joules_roundtrip():
    """Verify roundtrip conversion between eV and Joules."""
    energies_ev = np.array([0.0, 1.0, 1.12, 5.5, -2.3])
    joules = ev_to_joules(energies_ev)
    assert np.allclose(joules, energies_ev * Q)
    recovered_ev = joules_to_ev(joules)
    assert np.allclose(recovered_ev, energies_ev)


def test_potential_to_energy_ev_definition():
    """Verify EC(x) = EC_ref - q*phi(x) produces correct energy in eV.
    
    In eV units: E_eV = E_ref_eV - phi(x)[V].
    A positive electrostatic potential lowers electron energy.
    """
    phi_v = np.array([0.0, 0.5, 1.0, -0.5])
    ec_ev = potential_to_energy_ev(phi_v, ref_energy_ev=0.0)
    assert np.allclose(ec_ev, [0.0, -0.5, -1.0, 0.5])

    # Test with custom reference level
    ec_ev_shifted = potential_to_energy_ev(phi_v, ref_energy_ev=1.5)
    assert np.allclose(ec_ev_shifted, [1.5, 1.0, 0.5, 2.0])


def test_concentration_conversions():
    """Verify 1 cm^-3 = 10^6 m^-3."""
    c_cm3 = 1.0e16
    c_m3 = cm3_to_m3(c_cm3)
    assert c_m3 == pytest.approx(1.0e22)
    assert m3_to_cm3(c_m3) == pytest.approx(c_cm3)

    # Array support
    arr_cm3 = np.array([1e15, 1e18, 5e19])
    assert np.allclose(m3_to_cm3(cm3_to_m3(arr_cm3)), arr_cm3)


def test_length_conversions():
    """Verify length unit conversions."""
    # um to m
    assert um_to_m(1.0) == pytest.approx(1e-6)
    assert m_to_um(1e-6) == pytest.approx(1.0)

    # nm to m
    assert nm_to_m(10.0) == pytest.approx(1e-8)
    assert m_to_nm(1e-8) == pytest.approx(10.0)

    # cm to m
    assert cm_to_m(1.0) == pytest.approx(0.01)
    assert m_to_cm(0.01) == pytest.approx(1.0)

