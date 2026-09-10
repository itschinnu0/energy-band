"""Shared 1-D electrostatics foundation.

Implements core Maxwell/Poisson electrostatics relations in 1-D:
    dE/dx = rho / epsilon
    E = -d(phi)/dx
    d^2(phi)/dx^2 = -rho / epsilon

All internal variables must strictly adhere to SI units:
    position x: [m]
    potential phi: [V]
    electric field E: [V/m]
    charge density rho: [C/m^3]
    permittivity epsilon: [F/m]
"""

from typing import Union, Optional
import numpy as np
from scipy import integrate

from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.validation import validate_grid, validate_finite_array, validate_permittivity
from physics.constants import EPSILON_0, EPS_SI


def electric_field_from_potential(
    x: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """Calculate 1-D electric field E(x) from electrostatic potential phi(x).

    Relation:
        E(x) = -d(phi)/dx

    Uses NumPy gradient with second-order central differences in interior points
    and first-order one-sided differences at boundaries, taking arbitrary grid spacing into account.

    Args:
        x: 1-D spatial coordinate array in meters [m].
        phi: 1-D electrostatic potential array in Volts [V].

    Returns:
        Electric field array E(x) in Volts per meter [V/m].

    Raises:
        InvalidParameterError: If inputs are invalid or grids do not match.
        ModelValidityError: If inputs contain non-finite numbers.
    """
    validate_grid(x, min_points=2)
    if not isinstance(phi, np.ndarray):
        phi = np.asarray(phi, dtype=float)

    if phi.ndim != 1 or phi.shape != x.shape:
        raise InvalidParameterError(
            f"Potential array shape {phi.shape} must match grid array shape {x.shape} as 1-D"
        )
    validate_finite_array(phi, name="Potential phi")

    # E = -d(phi)/dx
    # edge_order=2 uses second-order accurate differences at boundaries
    dphi_dx = np.gradient(phi, x, edge_order=2)
    return -dphi_dx


def potential_from_electric_field(
    x: np.ndarray,
    electric_field: np.ndarray,
    x0: float,
    phi0: float,
) -> np.ndarray:
    """Calculate electrostatic potential phi(x) from electric field E(x).

    Relation:
        d(phi)/dx = -E(x)
        phi(x) = phi(x0) - integral_{x0}^x E(x') dx'

    Args:
        x: 1-D spatial coordinate array in meters [m].
        electric_field: 1-D electric field array in Volts per meter [V/m].
        x0: Reference coordinate in meters [m], must lie within or at the boundary of x.
        phi0: Reference boundary potential phi(x0) in Volts [V].

    Returns:
        Electrostatic potential array phi(x) in Volts [V].

    Raises:
        InvalidParameterError: If boundary reference x0 is outside the grid, or dimensions mismatch.
        ModelValidityError: If array values are non-finite.
    """
    validate_grid(x, min_points=2)
    if not isinstance(electric_field, np.ndarray):
        electric_field = np.asarray(electric_field, dtype=float)

    if electric_field.ndim != 1 or electric_field.shape != x.shape:
        raise InvalidParameterError(
            f"Electric field shape {electric_field.shape} must match grid array shape {x.shape} as 1-D"
        )
    validate_finite_array(electric_field, name="Electric field")

    if not np.isfinite(x0) or not np.isfinite(phi0):
        raise InvalidParameterError(f"Boundary condition (x0={x0}, phi0={phi0}) must be finite real numbers.")

    if x0 < x[0] - 1e-14 or x0 > x[-1] + 1e-14:
        raise InvalidParameterError(
            f"Reference coordinate x0={x0} m is outside the grid range [{x[0]}, {x[-1]}] m"
        )

    # Compute cumulative integral from grid start x[0]:
    # int_E[i] = integral_{x[0]}^{x[i]} E(x') dx'
    int_e_from_start = integrate.cumulative_trapezoid(electric_field, x, initial=0.0)

    # Value of integral at x0
    int_e_at_x0 = float(np.interp(x0, x, int_e_from_start))

    # phi(x) = phi0 - (int_{x[0]}^x E dx - int_{x[0]}^x0 E dx)
    phi = phi0 - (int_e_from_start - int_e_at_x0)
    return phi


def electric_field_from_charge_density(
    x: np.ndarray,
    rho: np.ndarray,
    epsilon: float,
    x0: float,
    e0: float,
) -> np.ndarray:
    """Calculate electric field E(x) from charge density rho(x).

    Relation:
        dE/dx = rho(x) / epsilon
        E(x) = E(x0) + integral_{x0}^x (rho(x') / epsilon) dx'

    Args:
        x: 1-D spatial coordinate array in meters [m].
        rho: 1-D charge density array in Coulombs per cubic meter [C/m^3].
        epsilon: Dielectric permittivity in Farads per meter [F/m].
        x0: Reference coordinate in meters [m], must lie within or at the boundary of x.
        e0: Reference boundary electric field E(x0) in Volts per meter [V/m].

    Returns:
        Electric field array E(x) in Volts per meter [V/m].

    Raises:
        InvalidParameterError: If inputs/parameters are invalid or out of bounds.
        ModelValidityError: If array values are non-finite.
    """
    validate_grid(x, min_points=2)
    validate_permittivity(epsilon)

    if not isinstance(rho, np.ndarray):
        rho = np.asarray(rho, dtype=float)

    if rho.ndim != 1 or rho.shape != x.shape:
        raise InvalidParameterError(
            f"Charge density shape {rho.shape} must match grid array shape {x.shape} as 1-D"
        )
    validate_finite_array(rho, name="Charge density")

    if not np.isfinite(x0) or not np.isfinite(e0):
        raise InvalidParameterError(f"Boundary condition (x0={x0}, e0={e0}) must be finite real numbers.")

    if x0 < x[0] - 1e-14 or x0 > x[-1] + 1e-14:
        raise InvalidParameterError(
            f"Reference coordinate x0={x0} m is outside the grid range [{x[0]}, {x[-1]}] m"
        )

    # Integrand = rho / epsilon
    source = rho / epsilon
    int_source_from_start = integrate.cumulative_trapezoid(source, x, initial=0.0)
    int_source_at_x0 = float(np.interp(x0, x, int_source_from_start))

    # E(x) = E0 + integral_{x0}^x (rho/epsilon) dx'
    e_field = e0 + (int_source_from_start - int_source_at_x0)
    return e_field


def potential_from_charge_density(
    x: np.ndarray,
    rho: np.ndarray,
    epsilon: float,
    x0_e: float,
    e0: float,
    x0_phi: float,
    phi0: float,
) -> np.ndarray:
    """Calculate electrostatic potential phi(x) from charge density rho(x) via double integration.

    Relations:
        d^2(phi)/dx^2 = -rho / epsilon
        E(x) = e0 + integral_{x0_e}^x (rho/epsilon) dx'
        phi(x) = phi0 - integral_{x0_phi}^x E(x') dx'

    Args:
        x: 1-D spatial coordinate array in meters [m].
        rho: 1-D charge density in [C/m^3].
        epsilon: Permittivity in [F/m].
        x0_e: Reference coordinate for electric field boundary condition [m].
        e0: Boundary electric field E(x0_e) in [V/m].
        x0_phi: Reference coordinate for potential boundary condition [m].
        phi0: Boundary potential phi(x0_phi) in [V].

    Returns:
        Electrostatic potential array phi(x) in Volts [V].
    """
    e_field = electric_field_from_charge_density(x, rho, epsilon, x0=x0_e, e0=e0)
    phi = potential_from_electric_field(x, e_field, x0=x0_phi, phi0=phi0)
    return phi

