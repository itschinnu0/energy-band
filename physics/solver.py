"""Numerical methods and bounded root solving.

Provides bounded, bracketed root solving with explicit convergence checks,
residual reporting, and diagnostic metadata.
"""

from typing import Callable, Tuple, Dict, Any, Optional
import math
from scipy import optimize

from physics.exceptions import SolverConvergenceError, InvalidParameterError


def solve_bracketed_root(
    f: Callable[[float], float],
    bracket: Tuple[float, float],
    xtol: float = 1e-12,
    rtol: float = 1e-10,
    maxiter: int = 100,
    name: str = "Variable",
) -> Tuple[float, Dict[str, Any]]:
    """Solve for a scalar root f(x) == 0 within a bounded interval [a, b].

    Uses SciPy's Brent method (`brentq`) with full convergence checks.
    Rejects cases where the root cannot be bracketed or solving does not converge.
    Never silently clips or invents solutions.

    Args:
        f: Continuous real-valued function of a single float variable.
        bracket: 2-tuple of (a, b) defining the bounding interval.
        xtol: Absolute tolerance for termination.
        rtol: Relative tolerance for termination.
        maxiter: Maximum number of iterations allowed.
        name: Name of variable being solved (for error diagnostic messages).

    Returns:
        (root, diagnostics_dict) where diagnostics_dict contains:
            - "converged": bool
            - "iterations": int
            - "function_calls": int
            - "residual": float (|f(root)|)
            - "bracket": Tuple[float, float]

    Raises:
        InvalidParameterError: If bracket bounds are invalid, non-finite, or a >= b.
        SolverConvergenceError: If f(a) and f(b) do not bracket a zero crossing,
                                or if Brent's method fails to converge.
    """
    a, b = bracket
    if not (math.isfinite(a) and math.isfinite(b)):
        raise InvalidParameterError(f"Root bracket bounds must be finite numbers: [{a}, {b}]")

    if a >= b:
        raise InvalidParameterError(f"Lower bracket bound must be strictly less than upper bound: [{a}, {b}]")

    try:
        fa = f(a)
        fb = f(b)
    except Exception as e:
        raise SolverConvergenceError(f"Failed to evaluate objective function at bracket edges [{a}, {b}]: {e}")

    if not (math.isfinite(fa) and math.isfinite(fb)):
        raise SolverConvergenceError(f"Objective function returned non-finite value at bracket edges: f({a})={fa}, f({b})={fb}")

    # Check exact boundary matches
    if abs(fa) <= xtol:
        return a, {
            "converged": True,
            "iterations": 0,
            "function_calls": 1,
            "residual": abs(fa),
            "bracket": (a, b),
        }
    if abs(fb) <= xtol:
        return b, {
            "converged": True,
            "iterations": 0,
            "function_calls": 2,
            "residual": abs(fb),
            "bracket": (a, b),
        }

    # Verify signs opposite
    if fa * fb > 0.0:
        raise SolverConvergenceError(
            f"Root for {name} is not bracketed on [{a:.4e}, {b:.4e}]: "
            f"f(a) = {fa:.4e} and f(b) = {fb:.4e} have the same sign."
        )

    try:
        root, r = optimize.brentq(
            f,
            a,
            b,
            xtol=xtol,
            rtol=rtol,
            maxiter=maxiter,
            full_output=True,
            disp=False,
        )
    except Exception as e:
        raise SolverConvergenceError(f"Root solving for {name} failed on [{a}, {b}]: {e}")

    if not r.converged:
        raise SolverConvergenceError(
            f"Root solving for {name} did not converge within {maxiter} iterations on [{a}, {b}]."
        )

    f_val = abs(f(root))
    diagnostics = {
        "converged": bool(r.converged),
        "iterations": int(r.iterations),
        "function_calls": int(r.function_calls),
        "residual": float(f_val),
        "bracket": (a, b),
    }

    return float(root), diagnostics
