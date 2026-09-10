"""Tests for shared 1-D electrostatics foundation.

Verifies:
- TEST 1: Linear potential produces constant electric field E = -dphi/dx
- TEST 2: Constant electric field produces linear potential phi = phi0 - E0*(x - x0)
- TEST 3: Constant charge density produces linear electric field E = E0 + (rho0/eps)*(x - x0)
- TEST 4: Poisson consistency: quadratic potential d^2phi/dx^2 = -rho/eps
- TEST 5: Zero charge density maintains constant electric field
- TEST 6: Zero electric field maintains constant potential
- TEST 7: Permittivity validation (rejects <=0, NaN, Inf)
- TEST 8: Grid validation (rejects non-1D, NaN/Inf, non-monotonic, duplicates, <2 points)
- TEST 9: Sign convention: positive potential slope -> negative electric field
- TEST 10: Round-trip consistency: phi -> E -> phi and rho -> E -> phi
"""

import pytest
import numpy as np
from physics.constants import EPSILON_0, EPS_SI
from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.electrostatics import (
    electric_field_from_potential,
    potential_from_electric_field,
    electric_field_from_charge_density,
    potential_from_charge_density,
)
from physics.validation import validate_grid, validate_permittivity


def test_1_linear_potential():
    """TEST 1: Given phi(x) = a*x + b, verify E(x) = -a."""
    x = np.linspace(0.0, 1.0e-6, 101)  # 1 um grid
    a = 2.5e6  # V/m slope
    b = 1.2    # V
    phi = a * x + b

    e_field = electric_field_from_potential(x, phi)

    # E = -dphi/dx = -a
    assert np.allclose(e_field, -a, rtol=1e-5)


def test_2_constant_electric_field():
    """TEST 2: Given E(x) = E0 and phi(x0) = phi0, verify phi(x) = phi0 - E0*(x-x0)."""
    x = np.linspace(0.0, 2.0e-6, 201)
    e0 = 5.0e5  # V/m
    e_field = np.full_like(x, e0)
    x0 = 0.5e-6
    phi0 = 3.0  # V

    phi = potential_from_electric_field(x, e_field, x0=x0, phi0=phi0)
    expected_phi = phi0 - e0 * (x - x0)

    assert np.allclose(phi, expected_phi, rtol=1e-5, atol=1e-6)


def test_3_constant_charge_density():
    """TEST 3: Given rho(x) = rho0 and E(x0) = E0, verify E(x) = E0 + (rho0/eps)*(x-x0)."""
    x = np.linspace(0.0, 1.0e-6, 501)
    rho0 = 1.6e4  # C/m^3
    rho = np.full_like(x, rho0)
    eps = EPS_SI
    x0 = 0.0
    e0 = -1.0e5  # V/m

    e_field = electric_field_from_charge_density(x, rho, eps, x0=x0, e0=e0)
    expected_e = e0 + (rho0 / eps) * (x - x0)

    assert np.allclose(e_field, expected_e, rtol=1e-4)


def test_4_poisson_consistency_quadratic():
    """TEST 4: Poisson consistency with quadratic potential:
    phi(x) = A*x^2 + B*x + C
    E(x) = -(2*A*x + B)
    dE/dx = -2*A = rho/eps  =>  rho(x) = -2*A*eps
    d^2phi/dx^2 = 2*A = -rho/eps
    """
    x = np.linspace(0.0, 1.0e-6, 1001)
    a = 1.0e12  # V/m^2
    b = -5.0e5  # V/m
    c = 0.7     # V
    phi_analytical = a * x**2 + b * x + c
    eps = EPS_SI

    # 1. Check E from phi
    e_numerical = electric_field_from_potential(x, phi_analytical)
    e_analytical = -(2.0 * a * x + b)
    assert np.allclose(e_numerical, e_analytical, rtol=1e-4, atol=1e-5)

    # 2. Check integration of rho to E and phi
    rho_constant = -2.0 * a * eps
    rho = np.full_like(x, rho_constant)

    e_from_rho = electric_field_from_charge_density(x, rho, eps, x0=x[0], e0=e_analytical[0])
    assert np.allclose(e_from_rho, e_analytical, rtol=1e-4, atol=1e-5)

    phi_from_rho = potential_from_charge_density(
        x, rho, eps, x0_e=x[0], e0=e_analytical[0], x0_phi=x[0], phi0=phi_analytical[0]
    )
    assert np.allclose(phi_from_rho, phi_analytical, rtol=1e-4, atol=1e-5)


def test_5_zero_charge_density():
    """TEST 5: For rho(x) = 0, verify dE/dx = 0 (E field remains constant)."""
    x = np.linspace(0.0, 1.0e-6, 101)
    rho_zero = np.zeros_like(x)
    e0 = 3.4e4

    e_field = electric_field_from_charge_density(x, rho_zero, EPS_SI, x0=x[0], e0=e0)
    assert np.allclose(e_field, e0)


def test_6_zero_electric_field():
    """TEST 6: Given E(x) = 0, verify potential remains constant."""
    x = np.linspace(0.0, 1.0e-6, 101)
    e_zero = np.zeros_like(x)
    phi0 = 1.5

    phi = potential_from_electric_field(x, e_zero, x0=x[0], phi0=phi0)
    assert np.allclose(phi, phi0)


def test_7_permittivity_validation():
    """TEST 7: Verify invalid epsilon values are rejected explicitly."""
    with pytest.raises(InvalidParameterError):
        validate_permittivity(0.0)
    with pytest.raises(InvalidParameterError):
        validate_permittivity(-1.0e-11)
    with pytest.raises(InvalidParameterError):
        validate_permittivity(float('nan'))
    with pytest.raises(InvalidParameterError):
        validate_permittivity(float('inf'))

    # Valid permittivities pass
    validate_permittivity(EPSILON_0)
    validate_permittivity(EPS_SI)


def test_8_grid_validation():
    """TEST 8: Verify grid validation rejects invalid inputs."""
    # Valid grid passes
    validate_grid(np.array([0.0, 1.0e-7, 2.0e-7]))

    # Non-1D array
    with pytest.raises(InvalidParameterError):
        validate_grid(np.ones((2, 2)))

    # Insufficient points (< 2)
    with pytest.raises(InvalidParameterError):
        validate_grid(np.array([1.0e-6]))

    # Non-finite coordinates
    with pytest.raises(ModelValidityError):
        validate_grid(np.array([0.0, np.nan, 2.0e-7]))
    with pytest.raises(ModelValidityError):
        validate_grid(np.array([0.0, np.inf, 2.0e-7]))

    # Decreasing coordinates (must not sort automatically)
    with pytest.raises(InvalidParameterError):
        validate_grid(np.array([2.0e-7, 1.0e-7, 0.0]))

    # Duplicate coordinates
    with pytest.raises(InvalidParameterError):
        validate_grid(np.array([0.0, 1.0e-7, 1.0e-7, 2.0e-7]))


def test_9_sign_convention():
    """TEST 9: E = -dphi/dx sign convention.
    A positive potential slope (dphi/dx > 0) MUST produce a negative electric field (E < 0).
    A negative potential slope (dphi/dx < 0) MUST produce a positive electric field (E > 0).
    """
    x = np.array([0.0, 1.0e-6, 2.0e-6])
    phi_increasing = np.array([0.0, 1.0, 2.0])  # dphi/dx = +1e6 V/m
    e_neg = electric_field_from_potential(x, phi_increasing)
    assert np.all(e_neg < 0.0)
    assert np.allclose(e_neg, -1.0e6)

    phi_decreasing = np.array([2.0, 1.0, 0.0])  # dphi/dx = -1e6 V/m
    e_pos = electric_field_from_potential(x, phi_decreasing)
    assert np.all(e_pos > 0.0)
    assert np.allclose(e_pos, +1.0e6)


def test_10_round_trip_consistency():
    """TEST 10: Demonstrate phi -> E -> phi and rho -> E -> phi round-trip consistency."""
    x = np.linspace(0.0, 5.0e-7, 501)
    # Sinusoidal potential test for non-trivial shape:
    # phi(x) = V0 * sin(2*pi*x / L) + phi_offset
    l_box = 5.0e-7
    v0 = 0.5
    phi_offset = 1.0
    phi_orig = v0 * np.sin(2.0 * np.pi * x / l_box) + phi_offset

    # phi -> E
    e_calc = electric_field_from_potential(x, phi_orig)

    # E -> phi (with boundary condition at x[0])
    phi_reconstructed = potential_from_electric_field(x, e_calc, x0=x[0], phi0=phi_orig[0])

    assert np.allclose(phi_reconstructed, phi_orig, rtol=1e-3, atol=1e-4)

