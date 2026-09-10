"""Tests for 1-D vertical MOS electrostatics and NMOS band model."""

import math
import pytest
import numpy as np

from physics.constants import EPS_OX, EPS_SI, Q
from physics.exceptions import InvalidParameterError, ModelValidityError, SolverConvergenceError
from physics.mos import (
    calculate_oxide_capacitance,
    calculate_fermi_potential,
    calculate_flat_band_voltage,
    calculate_level1_threshold_voltage,
    calculate_depletion_width,
    classify_mos_operating_state,
    calculate_semiconductor_charge_level2,
    solve_surface_potential_level2,
    simulate_mos,
)


def test_oxide_capacitance_equation():
    """Test 1: Cox = epsilon_ox / tox."""
    tox = 10e-9  # 10 nm
    expected_cox = EPS_OX / tox
    assert math.isclose(calculate_oxide_capacitance(tox), expected_cox, rel_tol=1e-9)


def test_invalid_tox_rejected():
    """Test 2: Invalid tox (<= 0 or NaN) is rejected without silent clipping."""
    with pytest.raises(InvalidParameterError):
        calculate_oxide_capacitance(0.0)
    with pytest.raises(InvalidParameterError):
        calculate_oxide_capacitance(-5e-9)
    with pytest.raises(InvalidParameterError):
        calculate_oxide_capacitance(float("nan"))


def test_fermi_potential_sign_for_ptype():
    """Test 3: phi_F has correct positive sign for p-type substrate."""
    phi_f = calculate_fermi_potential(na_cm3=1e17, temp_k=300.0)
    # For NA = 1e17 cm^-3 in Si at 300 K, phi_F is around 0.41 - 0.43 V
    assert phi_f > 0.0
    assert 0.35 < phi_f < 0.50


def test_flat_band_voltage_equation():
    """Test 4: VFB = PhiMS - Qox/Cox."""
    phims = -0.95
    qox = 1e-4  # C/m^2 positive oxide charge
    cox = 3.45e-3  # F/m^2
    expected_vfb = phims - (qox / cox)
    assert math.isclose(calculate_flat_band_voltage(phims, qox, cox), expected_vfb, rel_tol=1e-9)


def test_level1_threshold_equation():
    """Test 5: Level-1 threshold voltage VT = VFB + 2*phi_F + sqrt(4*q*eps_si*NA*phi_F)/Cox."""
    vfb = -0.9
    phi_f = 0.35
    na = 1e16
    cox = 3.45e-3
    v_th = calculate_level1_threshold_voltage(vfb, phi_f, na, cox)
    expected = vfb + 2.0 * phi_f + math.sqrt(4.0 * Q * EPS_SI * (na * 1e6) * phi_f) / cox
    assert math.isclose(v_th, expected, rel_tol=1e-9)
    assert v_th > vfb + 2.0 * phi_f


def test_gate_voltage_relation_sign_convention():
    """Test 6: VG = VFB + psi_s - Qs/Cox sign convention."""
    res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=1.5, level=1)
    p = res.calculated_parameters
    vg = p["VG_V"]
    vfb = p["VFB_V"]
    psi_s = p["psi_s_V"]
    qs = p["Qs_C_per_m2"]
    cox = p["Cox_F_per_m2"]

    # In depletion/inversion, Qs < 0, so -Qs/Cox > 0, which contributes positively to VG
    assert qs < 0.0
    vg_recalc = vfb + psi_s - (qs / cox)
    assert math.isclose(vg, vg_recalc, rel_tol=1e-5)


def test_operating_state_classification():
    """Tests 7, 8, 9, 10: Classification of accumulation, flat-band, depletion, strong-inversion."""
    vfb = -1.0
    vth = 0.8

    assert classify_mos_operating_state(vg_v=-2.0, vfb_v=vfb, threshold_voltage_v=vth) == "Accumulation"
    assert classify_mos_operating_state(vg_v=-1.0, vfb_v=vfb, threshold_voltage_v=vth) == "Flat band"
    assert classify_mos_operating_state(vg_v=-0.9995, vfb_v=vfb, threshold_voltage_v=vth) == "Flat band"  # within tol
    assert classify_mos_operating_state(vg_v=0.0, vfb_v=vfb, threshold_voltage_v=vth) == "Depletion"
    assert classify_mos_operating_state(vg_v=0.7995, vfb_v=vfb, threshold_voltage_v=vth) == "Threshold / transition"
    assert classify_mos_operating_state(vg_v=2.0, vfb_v=vfb, threshold_voltage_v=vth) == "Strong inversion"


def test_depletion_width_equation():
    """Test 11: Wd = sqrt(2*eps_si*psi_s / (q*NA))."""
    psi_s = 0.4
    na = 1e16
    wd = calculate_depletion_width(psi_s_v=psi_s, na_cm3=na)
    expected = math.sqrt(2.0 * EPS_SI * psi_s / (Q * (na * 1e6)))
    assert math.isclose(wd, expected, rel_tol=1e-9)

    # Negative psi_s in accumulation raises ModelValidityError
    with pytest.raises(ModelValidityError, match="negative psi_s"):
        calculate_depletion_width(psi_s_v=-0.1, na_cm3=na)


def test_qs_sign_convention():
    """Test 12: Qs > 0 for accumulation (psi_s < 0), Qs < 0 for depletion/inversion (psi_s > 0)."""
    phi_f = 0.35
    na = 1e16
    qs_acc = calculate_semiconductor_charge_level2(psi_s_v=-0.2, na_cm3=na, phi_f_v=phi_f)
    qs_dep = calculate_semiconductor_charge_level2(psi_s_v=0.3, na_cm3=na, phi_f_v=phi_f)
    qs_inv = calculate_semiconductor_charge_level2(psi_s_v=0.8, na_cm3=na, phi_f_v=phi_f)

    assert qs_acc > 0.0
    assert qs_dep < 0.0
    assert qs_inv < qs_dep  # stronger negative inversion charge


def test_psi_s_zero_gives_flat_bands():
    """Test 13: psi_s = 0 gives flat potential and flat bands in silicon."""
    res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=-0.9, phims_v=-0.9, qox_c_m2=0.0, level=1)
    assert math.isclose(res.calculated_parameters["psi_s_V"], 0.0, abs_tol=1e-9)
    assert np.allclose(res.potential, 0.0, atol=1e-12)
    assert np.allclose(np.diff(res.EC), 0.0, atol=1e-12)
    assert np.allclose(np.diff(res.EV), 0.0, atol=1e-12)


def test_bandgap_consistency_mos():
    """Test 14: EC - EV = Eg everywhere in silicon."""
    res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=2.0, level=1)
    eg = res.calculated_parameters["bandgap_eV"]
    np.testing.assert_allclose(res.EC - res.EV, eg, atol=1e-9)
    np.testing.assert_allclose(res.Ei, res.EC - 0.5 * eg, atol=1e-9)


def test_bulk_potential_reference_is_zero():
    """Test 15: Bulk potential reference is zero (phi_bulk = 0 V)."""
    res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=1.5, level=1)
    # The last point in deep silicon bulk must be 0 V
    assert math.isclose(res.potential[-1], 0.0, abs_tol=1e-9)


def test_band_bending_direction():
    """Test 16: psi_s > 0 produces downward band bending from bulk to surface."""
    res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=1.5, level=1)
    psi_s = res.calculated_parameters["psi_s_V"]
    assert psi_s > 0.0
    # Surface is at index 0, bulk is at index -1
    # phi(0) = psi_s > 0, phi(bulk) = 0
    # EC(0) = EC_bulk - psi_s < EC(bulk) -> downward band bending toward surface
    assert res.EC[0] < res.EC[-1]
    assert math.isclose(res.EC[-1] - res.EC[0], psi_s, rel_tol=1e-6)


def test_level2_numerical_solver_converges():
    """Test 17 & 18: Level-2 numerical solver converges for normal cases and residual is within tolerance."""
    for vg in [-1.5, -0.5, 0.0, 1.0, 2.5]:
        res = simulate_mos(na_cm3=1e16, tox_m=10e-9, vg_v=vg, level=2)
        assert res.model_valid is True
        diag = res.calculated_parameters["solver_diagnostics"]
        assert diag["converged"] is True
        assert diag["residual"] < 1e-8


def test_solver_failure_raises_solver_convergence_error():
    """Test 19: Unbracketable or impossible solver conditions raise SolverConvergenceError."""
    # Test solve_surface_potential_level2 with an impossible bracket range where f(a)*f(b) > 0
    with pytest.raises(SolverConvergenceError):
        solve_surface_potential_level2(
            vg_v=5.0,
            vfb_v=0.0,
            cox_f_m2=1e-3,
            na_cm3=1e16,
            phi_f_v=0.35,
            bracket_range=(1.5, 2.0),  # Both produce values far from 0
        )


def test_no_silent_clipping_and_input_rejection():
    """Test 20: Invalid parameters raise exceptions rather than being silently clipped."""
    with pytest.raises(InvalidParameterError):
        simulate_mos(na_cm3=-1e16)

    with pytest.raises(InvalidParameterError):
        simulate_mos(temp_k=10.0)  # Below min 50K

    with pytest.raises(InvalidParameterError):
        simulate_mos(tox_m=-10e-9)

    with pytest.raises(InvalidParameterError):
        simulate_mos(level=3)
