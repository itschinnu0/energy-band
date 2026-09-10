"""Tests for semiconductor physics foundation calculations."""

import pytest
import numpy as np
from physics.constants import EG_0_EV, T_REF_K
from physics.exceptions import InvalidParameterError
from physics.semiconductor import (
    thermal_voltage,
    silicon_bandgap_ev,
    silicon_bandgap_joules,
    effective_density_of_states,
    intrinsic_carrier_concentration,
    get_semiconductor_properties,
    equilibrium_fermi_level,
    carrier_densities_from_quasi_fermi,
    bands_from_potential,
)
from physics.units import joules_to_ev


def test_thermal_voltage():
    """Verify thermal voltage at T = 300 K is ~25.85 mV."""
    vt = thermal_voltage(300.0)
    assert vt == pytest.approx(0.025852, rel=1e-3)

    # Scaling with temperature: Vt(600) = 2 * Vt(300)
    assert thermal_voltage(600.0) == pytest.approx(2.0 * vt, rel=1e-5)


def test_silicon_bandgap_temperature_dependence():
    """Verify bandgap Varshni relation: Eg(0 K) = 1.166 eV, Eg(300 K) ~ 1.125 eV."""
    # Eg at T = 300 K
    eg_300 = silicon_bandgap_ev(300.0)
    assert 1.11 < eg_300 < 1.14
    assert eg_300 == pytest.approx(1.1246, abs=0.005)

    # Eg in Joules roundtrip
    eg_j = silicon_bandgap_joules(300.0)
    assert joules_to_ev(eg_j) == pytest.approx(eg_300)

    # Monotonically decreasing with temperature
    assert silicon_bandgap_ev(100.0) > silicon_bandgap_ev(300.0)
    assert silicon_bandgap_ev(300.0) > silicon_bandgap_ev(400.0)


def test_effective_density_of_states():
    """Verify Nc and Nv at 300 K and temperature scaling T^(3/2)."""
    nc_300, nv_300 = effective_density_of_states(300.0)
    assert nc_300 == pytest.approx(2.86e19)
    assert nv_300 == pytest.approx(2.66e19)

    nc_600, nv_600 = effective_density_of_states(600.0)
    factor = (600.0 / 300.0) ** 1.5
    assert nc_600 == pytest.approx(nc_300 * factor)
    assert nv_600 == pytest.approx(nv_300 * factor)


def test_intrinsic_carrier_concentration():
    """Verify intrinsic carrier concentration ni in silicon at 300 K (~1e10 cm^-3)."""
    ni = intrinsic_carrier_concentration(300.0)
    # Typically ~ 1.0e10 cm^-3 using Sze parameters
    assert 5e9 < ni < 2e10

    # Temperature sensitivity: ni increases sharply with T
    ni_400 = intrinsic_carrier_concentration(400.0)
    assert ni_400 > ni * 100


def test_equilibrium_fermi_level():
    """Verify equilibrium Fermi level position for n-type and p-type silicon."""
    temp = 300.0
    ni = intrinsic_carrier_concentration(temp)
    vt = thermal_voltage(temp)

    # For n-type with ND = 1e16 cm^-3: EF > Ei
    ef_n = equilibrium_fermi_level(1e16, 'n', temp_k=temp, ei_ev=0.0)
    expected_ef_n = vt * np.log(1e16 / ni)
    assert ef_n == pytest.approx(expected_ef_n, rel=1e-4)
    assert ef_n > 0.0

    # For p-type with NA = 1e16 cm^-3: EF < Ei
    ef_p = equilibrium_fermi_level(1e16, 'p', temp_k=temp, ei_ev=0.0)
    expected_ef_p = -vt * np.log(1e16 / ni)
    assert ef_p == pytest.approx(expected_ef_p, rel=1e-4)
    assert ef_p < 0.0

    # Symmetry for equal doping
    assert ef_n == pytest.approx(-ef_p)

    # Invalid doping type
    with pytest.raises(ValueError):
        equilibrium_fermi_level(1e16, 'invalid')


def test_carrier_densities_from_quasi_fermi():
    """Verify Boltzmann quasi-Fermi carrier relationships:
    n = ni * exp((EFn - Ei) / kT)
    p = ni * exp((Ei - EFp) / kT)
    np = ni^2 * exp((EFn - EFp) / kT)
    """
    temp = 300.0
    ni = intrinsic_carrier_concentration(temp)
    vt = thermal_voltage(temp)

    # Equilibrium case: EFn = EFp = EF
    ef = 0.3  # eV above Ei
    n_eq, p_eq = carrier_densities_from_quasi_fermi(ef, ef, 0.0, temp_k=temp)
    assert n_eq * p_eq == pytest.approx(ni**2, rel=1e-4)
    assert n_eq == pytest.approx(ni * np.exp(0.3 / vt), rel=1e-4)

    # Non-equilibrium case: EFn - EFp = 0.2 eV
    efn = 0.3
    efp = 0.1
    n, p = carrier_densities_from_quasi_fermi(efn, efp, 0.0, temp_k=temp)
    assert n * p == pytest.approx(ni**2 * np.exp(0.2 / vt), rel=1e-4)


def test_bands_from_potential_and_ec_ev_identity():
    """Verify EC - EV = Eg identity and EC = EC_ref - phi(x) relation."""
    eg = 1.12
    phi = np.array([0.0, 0.2, 0.5, 0.8, -0.4])
    ec, ev, ei = bands_from_potential(phi, eg_ev=eg, ec_ref_ev=0.0)

    # Identity 1: EC - EV == Eg everywhere
    assert np.allclose(ec - ev, eg)

    # Identity 2: Ei is exactly mid-gap
    assert np.allclose(ei, (ec + ev) / 2.0)

    # Identity 3: EC = -phi for ec_ref = 0
    assert np.allclose(ec, -phi)

    # Downward band bending for positive potential
    # If phi increases, EC and EV decrease (bend downwards)
    assert ec[1] < ec[0]
    assert ev[1] < ev[0]


def test_get_semiconductor_properties():
    """Verify package semiconductor property helper."""
    props = get_semiconductor_properties(300.0)
    assert props.temperature_k == 300.0
    assert props.bandgap_ev == pytest.approx(silicon_bandgap_ev(300.0))
    assert props.ni_m3 == pytest.approx(props.ni_cm3 * 1e6)


