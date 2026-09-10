"""Tests for 1-D abrupt PN junction physics model."""

import math
import pytest
import numpy as np

from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.constants import EPS_SI, Q
from physics.pn_junction import (
    calculate_pn_built_in_potential,
    calculate_pn_depletion_widths,
    simulate_pn_junction,
)


def test_vbi_positive_for_physical_doping():
    """Test 1: Vbi is positive for physical doping values."""
    vbi = calculate_pn_built_in_potential(na_cm3=1e16, nd_cm3=1e16, temp_k=300.0)
    assert vbi > 0.0
    # For Si at 300K with NA=ND=1e16, Vbi ~ 0.65 - 0.75 V
    assert 0.6 < vbi < 0.8


def test_va_zero_gives_vd_equals_vbi():
    """Test 2: VA=0 gives VD=Vbi."""
    res = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0, temp_k=300.0)
    vbi = res.calculated_parameters["Vbi_V"]
    vd = res.calculated_parameters["VD_V"]
    assert math.isclose(vd, vbi, rel_tol=1e-9)


def test_forward_bias_reduces_vd():
    """Test 3: Forward bias (VA > 0) reduces VD."""
    res_eq = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0, temp_k=300.0)
    res_fwd = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.3, temp_k=300.0)
    assert res_fwd.calculated_parameters["VD_V"] < res_eq.calculated_parameters["VD_V"]
    assert math.isclose(
        res_fwd.calculated_parameters["VD_V"],
        res_eq.calculated_parameters["Vbi_V"] - 0.3,
        rel_tol=1e-9,
    )


def test_reverse_bias_increases_vd():
    """Test 4: Reverse bias (VA < 0) increases VD."""
    res_eq = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0, temp_k=300.0)
    res_rev = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=-2.0, temp_k=300.0)
    assert res_rev.calculated_parameters["VD_V"] > res_eq.calculated_parameters["VD_V"]
    assert math.isclose(
        res_rev.calculated_parameters["VD_V"],
        res_eq.calculated_parameters["Vbi_V"] + 2.0,
        rel_tol=1e-9,
    )


def test_depletion_width_changes_with_vd():
    """Test 5: W changes correctly with VD (proportional to sqrt(VD))."""
    res_eq = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0, temp_k=300.0)
    res_rev = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=-2.0, temp_k=300.0)
    res_fwd = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.3, temp_k=300.0)

    w_eq = res_eq.calculated_parameters["W_m"]
    w_rev = res_rev.calculated_parameters["W_m"]
    w_fwd = res_fwd.calculated_parameters["W_m"]

    assert w_rev > w_eq > w_fwd

    # Exact scaling: W_rev / W_eq == sqrt(VD_rev / VD_eq)
    vd_eq = res_eq.calculated_parameters["VD_V"]
    vd_rev = res_rev.calculated_parameters["VD_V"]
    expected_ratio = math.sqrt(vd_rev / vd_eq)
    actual_ratio = w_rev / w_eq
    assert math.isclose(actual_ratio, expected_ratio, rel_tol=1e-6)


def test_charge_neutrality():
    """Test 6: Charge neutrality NA*xp = ND*xn."""
    # Asymmetric junction
    na = 5e16
    nd = 2e15
    res = simulate_pn_junction(na_cm3=na, nd_cm3=nd, va_v=0.0, temp_k=300.0)
    xp = res.calculated_parameters["xp_m"]
    xn = res.calculated_parameters["xn_m"]

    # More lightly doped side (n) has larger depletion width
    assert xn > xp
    assert math.isclose(na * xp, nd * xn, rel_tol=1e-9)


def test_potential_drop_across_depletion():
    """Test 7: Potential drop across depletion = VD."""
    va = -1.5
    res = simulate_pn_junction(na_cm3=2e16, nd_cm3=1e16, va_v=va, temp_k=300.0)
    vd = res.calculated_parameters["VD_V"]
    phi = res.potential
    # Drop between n-neutral and p-neutral
    delta_phi = phi[-1] - phi[0]
    assert math.isclose(delta_phi, vd, rel_tol=1e-6)


def test_correct_electric_field_sign():
    """Test 8: Correct electric field sign and peak magnitude."""
    res = simulate_pn_junction(na_cm3=1e16, nd_cm3=2e16, va_v=-1.0, temp_k=300.0)
    e_field = res.electric_field
    # Electric field points in -x direction (from + charge in n to - charge in p)
    # So E(x) <= 0 everywhere
    assert np.all(e_field <= 1e-12)

    # Maximum magnitude at x=0
    peak_e = res.calculated_parameters["peak_electric_field_V_per_m"]
    assert peak_e < 0.0
    # Peak analytical formula: q * NA * xp / eps_si
    na_m3 = 1e16 * 1e6
    xp = res.calculated_parameters["xp_m"]
    expected_peak = - (Q * na_m3 * xp) / EPS_SI
    assert math.isclose(peak_e, expected_peak, rel_tol=1e-5)


def test_bandgap_consistency():
    """Test 9: EC - EV = Eg strictly everywhere."""
    res = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.2, temp_k=300.0)
    eg = res.calculated_parameters["bandgap_eV"]
    diff = res.EC - res.EV
    np.testing.assert_allclose(diff, eg, atol=1e-9)

    ei_expected = res.EC - 0.5 * eg
    np.testing.assert_allclose(res.Ei, ei_expected, atol=1e-9)


def test_equilibrium_ef_is_flat():
    """Test 10: Equilibrium EF is flat and non-None."""
    res = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0, temp_k=300.0)
    assert res.EF is not None
    assert np.all(np.isfinite(res.EF))
    assert np.max(res.EF) - np.min(res.EF) < 1e-12
    # In equilibrium, EFn and EFp should be None
    assert res.EFn is None
    assert res.EFp is None


def test_biased_model_quasi_fermi_rules():
    """Test 11: Biased model does not fabricate a global flat EF and enforces EFn_n - EFp_p = VA."""
    va = 0.4
    res = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=va, temp_k=300.0)
    # Must NOT have global EF
    assert res.EF is None

    # Must have EFn and EFp defined in neutral regions
    assert res.EFn is not None
    assert res.EFp is not None

    xp = res.calculated_parameters["xp_m"]
    xn = res.calculated_parameters["xn_m"]

    # In p-neutral (x < -xp), EFp = -VA
    p_neutral_mask = res.position < -xp
    np.testing.assert_allclose(res.EFp[p_neutral_mask], -va, atol=1e-9)

    # In n-neutral (x > xn), EFn = 0
    n_neutral_mask = res.position > xn
    np.testing.assert_allclose(res.EFn[n_neutral_mask], 0.0, atol=1e-9)

    # Enforce EFn_n - EFp_p = VA
    efn_n = float(np.nanmean(res.EFn[n_neutral_mask]))
    efp_p = float(np.nanmean(res.EFp[p_neutral_mask]))
    assert math.isclose(efn_n - efp_p, va, abs_tol=1e-9)

    # Depletion region must NOT have fabricated values
    dep_mask = (res.position >= -xp) & (res.position <= xn)
    assert np.all(np.isnan(res.EFn[dep_mask]))
    assert np.all(np.isnan(res.EFp[dep_mask]))


def test_vd_non_positive_raises_model_validity_error():
    """Test 12 & 13: VD <= 0 raises ModelValidityError without silent clipping."""
    # Built in potential for NA=ND=1e16 is ~0.7V. Setting VA = 1.0 V forces VD < 0
    with pytest.raises(ModelValidityError, match="outside its physical validity range"):
        simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=1.5, temp_k=300.0)

    # Also test calculate_pn_depletion_widths directly with VD <= 0
    with pytest.raises(ModelValidityError, match="VD"):
        calculate_pn_depletion_widths(na_cm3=1e16, nd_cm3=1e16, vd_v=0.0)

    with pytest.raises(ModelValidityError, match="VD"):
        calculate_pn_depletion_widths(na_cm3=1e16, nd_cm3=1e16, vd_v=-0.5)


def test_invalid_parameters_rejected():
    """Test 14: Invalid doping, temperature, and voltages are rejected."""
    with pytest.raises(InvalidParameterError):
        simulate_pn_junction(na_cm3=-1e16, nd_cm3=1e16)

    with pytest.raises(InvalidParameterError):
        simulate_pn_junction(na_cm3=1e16, nd_cm3=0.0)

    with pytest.raises(InvalidParameterError):
        simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, temp_k=20.0)  # Below 50 K

    with pytest.raises(InvalidParameterError):
        simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, temp_k=700.0)  # Above 600 K

    with pytest.raises(InvalidParameterError):
        simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=float("nan"))
