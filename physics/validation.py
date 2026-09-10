"""Parameter validation for semiconductor models.

All checks reject invalid physical states explicitly without silent clipping.
Raises InvalidParameterError or ModelValidityError.
"""

from typing import Optional
import math
import numpy as np
from physics.exceptions import InvalidParameterError, ModelValidityError


def validate_temperature(temp_k: float, min_temp_k: float = 50.0, max_temp_k: float = 600.0) -> None:
    """Validate operating temperature in Kelvin.
    
    Silicon semiconductor models with standard Boltzmann statistics and ionization
    are typically physically valid between ~50 K and ~600 K.
    """
    if not math.isfinite(temp_k):
        raise InvalidParameterError(f"Temperature must be a finite real number, got: {temp_k}")
    if temp_k <= 0:
        raise InvalidParameterError(f"Temperature must be strictly positive (T > 0 K), got: {temp_k} K")
    if temp_k < min_temp_k or temp_k > max_temp_k:
        raise InvalidParameterError(
            f"Temperature {temp_k} K is outside the supported physical range [{min_temp_k}, {max_temp_k}] K"
        )


def validate_doping(
    doping_cm3: float,
    doping_name: str = "Doping",
    min_doping_cm3: float = 1.0e12,
    max_doping_cm3: float = 1.0e20,
) -> None:
    """Validate doping concentration in cm^-3.
    
    Rejects negative, non-finite, or unphysical doping levels.
    """
    if not math.isfinite(doping_cm3):
        raise InvalidParameterError(f"{doping_name} concentration must be finite, got: {doping_cm3}")
    if doping_cm3 <= 0:
        raise InvalidParameterError(f"{doping_name} concentration must be strictly positive, got: {doping_cm3} cm^-3")
    if doping_cm3 < min_doping_cm3 or doping_cm3 > max_doping_cm3:
        raise InvalidParameterError(
            f"{doping_name} concentration {doping_cm3:e} cm^-3 is outside physical range "
            f"[{min_doping_cm3:e}, {max_doping_cm3:e}] cm^-3"
        )


def validate_finite_array(arr: np.ndarray, name: str = "Array") -> None:
    """Ensure array contains only finite floating point numbers (no NaN, no Inf)."""
    if not np.all(np.isfinite(arr)):
        raise ModelValidityError(f"{name} contains non-finite values (NaN or Inf)")


def validate_grid(x: np.ndarray, min_points: int = 2) -> None:
    """Validate 1-D spatial grid coordinate array.

    Requirements:
    - 1-D numpy array
    - Finite values (no NaN, no Inf)
    - At least min_points
    - Strictly increasing (dx > 0, no duplicate coordinates, no decreasing steps)
    """
    if not isinstance(x, np.ndarray):
        try:
            x = np.asarray(x, dtype=float)
        except Exception as e:
            raise InvalidParameterError(f"Grid must be convertible to a NumPy float array: {e}")

    if x.ndim != 1:
        raise InvalidParameterError(f"Spatial grid must be a 1-D array, got shape {x.shape}")

    if len(x) < min_points:
        raise InvalidParameterError(f"Spatial grid must have at least {min_points} points, got {len(x)}")

    validate_finite_array(x, name="Spatial grid")

    dx = np.diff(x)
    if np.any(dx <= 0):
        if np.any(dx == 0):
            raise InvalidParameterError("Spatial grid contains duplicate/repeated coordinates.")
        raise InvalidParameterError("Spatial grid must be strictly monotonically increasing.")


def validate_permittivity(epsilon: float, name: str = "Permittivity") -> None:
    """Validate dielectric permittivity [F/m].

    Must be finite and strictly positive (> 0).
    """
    if not math.isfinite(epsilon):
        raise InvalidParameterError(f"{name} must be a finite real number, got: {epsilon}")
    if epsilon <= 0:
        raise InvalidParameterError(f"{name} must be strictly positive (> 0 F/m), got: {epsilon}")

