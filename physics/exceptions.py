"""Domain-specific exceptions for energy band simulator."""


class EnergyBandError(Exception):
    """Base exception for all physics and simulation errors."""
    pass


class InvalidParameterError(EnergyBandError):
    """Raised when an input parameter fails validation or is outside physical limits."""
    pass


class ModelValidityError(EnergyBandError):
    """Raised when a physical condition violates the domain of validity for the model."""
    pass


class SolverConvergenceError(EnergyBandError):
    """Raised when a numerical solver fails to converge within tolerance."""
    pass


class UnitConversionError(EnergyBandError):
    """Raised when an invalid unit or unit conversion is requested."""
    pass
