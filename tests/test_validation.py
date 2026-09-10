"""Tests for parameter validation and error raising."""

import pytest
import numpy as np
from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.validation import validate_temperature, validate_doping, validate_finite_array


def test_validate_temperature_valid():
    """Verify valid temperatures pass validation."""
    validate_temperature(300.0)
    validate_temperature(77.0)
    validate_temperature(500.0)


def test_validate_temperature_invalid():
    """Verify invalid temperatures raise InvalidParameterError."""
    with pytest.raises(InvalidParameterError):
        validate_temperature(0.0)  # Absolute zero not allowed in Boltzmann
    with pytest.raises(InvalidParameterError):
        validate_temperature(-10.0)
    with pytest.raises(InvalidParameterError):
        validate_temperature(float('nan'))
    with pytest.raises(InvalidParameterError):
        validate_temperature(float('inf'))
    with pytest.raises(InvalidParameterError):
        validate_temperature(10.0)  # Below min 50 K
    with pytest.raises(InvalidParameterError):
        validate_temperature(1000.0)  # Above max 600 K


def test_validate_doping_valid():
    """Verify valid doping concentrations pass."""
    validate_doping(1.0e15)
    validate_doping(1.0e17)
    validate_doping(1.0e19)


def test_validate_doping_invalid():
    """Verify invalid doping raises InvalidParameterError."""
    with pytest.raises(InvalidParameterError):
        validate_doping(0.0)
    with pytest.raises(InvalidParameterError):
        validate_doping(-1.0e15)
    with pytest.raises(InvalidParameterError):
        validate_doping(float('nan'))
    with pytest.raises(InvalidParameterError):
        validate_doping(1.0e10)  # Below min 1e12
    with pytest.raises(InvalidParameterError):
        validate_doping(1.0e22)  # Above max 1e20


def test_validate_finite_array():
    """Verify finite array validation."""
    valid_arr = np.array([1.0, 2.0, 3.5, -4.0])
    validate_finite_array(valid_arr)

    with pytest.raises(ModelValidityError):
        validate_finite_array(np.array([1.0, np.nan, 3.0]))

    with pytest.raises(ModelValidityError):
        validate_finite_array(np.array([1.0, np.inf, 3.0]))

