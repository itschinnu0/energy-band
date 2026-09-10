"""Tests for NPN BJT two-junction analytical model."""

import math
import pytest
import numpy as np

from physics.exceptions import InvalidParameterError, ModelValidityError
from physics.bjt import classify_bjt_operating_region, simulate_npn_bjt


def test_npn_structure_and_regions():
    """Test 1: NPN structure contains emitter, base, and collector regions."""
    res = simulate_npn_bjt(ne_cm3=1e18, nb_cm3=1e16, nc_cm3=1e15, vbe_v=0.0, vbc_v=0.0)
    assert res.device == "NPN BJT"
    assert res.model == "Two-junction analytical electrostatic approximation"

    # Verify regions present
    regions = set(res.region)
    expected_regions = {
        "emitter-neutral",
        "eb-depletion-emitter",
        "eb-depletion-base",
        "base-neutral",
        "cb-depletion-base",
        "cb-depletion-collector",
        "collector-neutral",
    }
    assert expected_regions.issubset(regions)


def test_vbe_vbc_sign_conventions():
    """Test 2: VBE = VB - VE and VBC = VB - VC conventions."""
    res = simulate_npn_bjt(vbe_v=0.6, vbc_v=-2.0)
    params = res.calculated_parameters
    assert params["VBE_V"] == 0.6
    assert params["VBC_V"] == -2.0


def test_correct_junction_barrier_calculation():
    """Test 3: Correct junction barrier voltages Vbarrier = Vbi - V_applied."""
    vbe = 0.5
    vbc = -1.0
    res = simulate_npn_bjt(ne_cm3=1e18, nb_cm3=1e16, nc_cm3=1e15, vbe_v=vbe, vbc_v=vbc)
    params = res.calculated_parameters

    vbi_eb = params["Vbi_EB_V"]
    vbi_cb = params["Vbi_CB_V"]

    assert math.isclose(params["Vbarrier_EB_V"], vbi_eb - vbe, rel_tol=1e-9)
    assert math.isclose(params["Vbarrier_CB_V"], vbi_cb - vbc, rel_tol=1e-9)


def test_operating_region_classification():
    """Test 4, 5, 6, 7, 8, 9: Operating region classification for all modes."""
    # Cutoff: EB reverse, BC reverse
    assert classify_bjt_operating_region(vbe_v=-0.5, vbc_v=-1.0) == "Cutoff"

    # Forward active: EB forward, BC reverse
    assert classify_bjt_operating_region(vbe_v=0.65, vbc_v=-2.0) == "Forward active"

    # Saturation: EB forward, BC forward
    assert classify_bjt_operating_region(vbe_v=0.65, vbc_v=0.4) == "Saturation"

    # Reverse active: EB reverse, BC forward
    assert classify_bjt_operating_region(vbe_v=-0.5, vbc_v=0.5) == "Reverse active"

    # Equilibrium: both 0
    assert classify_bjt_operating_region(vbe_v=0.0, vbc_v=0.0) == "Equilibrium"

    # Near-zero transition: |V| <= 1 mV
    assert classify_bjt_operating_region(vbe_v=0.0005, vbc_v=-1.0) == "Transition / near-zero bias"
    assert classify_bjt_operating_region(vbe_v=0.6, vbc_v=-0.0005) == "Transition / near-zero bias"


def test_depletion_width_validity():
    """Test 10: Depletion widths are calculated correctly and scale with barrier."""
    res_fa = simulate_npn_bjt(ne_cm3=1e18, nb_cm3=1e16, nc_cm3=1e15, vbe_v=0.5, vbc_v=-2.0)
    p = res_fa.calculated_parameters

    # Emitter is N+ heavily doped, Base is P moderately doped -> x_be >> x_e
    assert p["x_be_m"] > p["x_e_m"]
    # Base is P moderately doped, Collector is N lightly doped -> x_c >> x_bc
    assert p["x_c_m"] > p["x_bc_m"]


def test_base_depletion_overlap_detection():
    """Test 11: Base-width overlap (punch-through) raises ModelValidityError without silent clipping."""
    # Small base width with high reverse bias on CB junction
    with pytest.raises(ModelValidityError, match="Base depletion overlap detected"):
        simulate_npn_bjt(
            ne_cm3=1e18,
            nb_cm3=1e15,  # light base doping gives large depletion width
            nc_cm3=1e15,
            vbe_v=-1.0,
            vbc_v=-10.0,  # high reverse bias expands depletion into base
            w_base_m=0.2e-6,  # 200 nm thin base
        )


def test_equilibrium_flat_ef():
    """Test 12: Equilibrium BJT has a single flat EF."""
    res = simulate_npn_bjt(vbe_v=0.0, vbc_v=0.0)
    assert res.EF is not None
    assert np.all(np.isfinite(res.EF))
    assert np.max(res.EF) - np.min(res.EF) < 1e-12
    assert res.EFn is None
    assert res.EFp is None


def test_bandgap_consistency_bjt():
    """Test 13: EC - EV = Eg strictly throughout BJT."""
    res = simulate_npn_bjt(vbe_v=0.5, vbc_v=-1.5)
    eg = res.calculated_parameters["bandgap_eV"]
    np.testing.assert_allclose(res.EC - res.EV, eg, atol=1e-9)

    ei_expected = res.EC - 0.5 * eg
    np.testing.assert_allclose(res.Ei, ei_expected, atol=1e-9)


def test_no_fabricated_quasi_fermi_transport():
    """Test 14: Biased BJT has quasi-Fermi levels ONLY in neutral regions and no fake transport profile."""
    vbe = 0.6
    vbc = -1.0
    res = simulate_npn_bjt(vbe_v=vbe, vbc_v=vbc)

    # Must NOT have global EF
    assert res.EF is None
    assert res.EFn is not None
    assert res.EFp is not None

    m_e_neut = res.region == "emitter-neutral"
    m_b_neut = res.region == "base-neutral"
    m_c_neut = res.region == "collector-neutral"

    # Neutral base EFp = 0
    np.testing.assert_allclose(res.EFp[m_b_neut], 0.0, atol=1e-9)
    # Neutral emitter EFn = -VBE
    np.testing.assert_allclose(res.EFn[m_e_neut], -vbe, atol=1e-9)
    # Neutral collector EFn = -VBC
    np.testing.assert_allclose(res.EFn[m_c_neut], -vbc, atol=1e-9)

    # In depletion regions, no fake profile
    dep_mask = np.isin(
        res.region,
        ["eb-depletion-emitter", "eb-depletion-base", "cb-depletion-base", "cb-depletion-collector"],
    )
    assert np.all(np.isnan(res.EFn[dep_mask]))
    assert np.all(np.isnan(res.EFp[dep_mask]))


def test_no_silent_clipping_and_invalid_inputs():
    """Test 15: Barrier <= 0 raises ModelValidityError; invalid inputs raise InvalidParameterError."""
    # Excessive forward bias on EB junction reduces barrier to <= 0
    with pytest.raises(ModelValidityError, match="EB junction barrier"):
        simulate_npn_bjt(vbe_v=2.0, vbc_v=0.0)

    # Excessive forward bias on CB junction reduces barrier to <= 0
    with pytest.raises(ModelValidityError, match="CB junction barrier"):
        simulate_npn_bjt(vbe_v=0.0, vbc_v=2.0)

    # Invalid doping
    with pytest.raises(InvalidParameterError):
        simulate_npn_bjt(ne_cm3=-1e18)

    # Invalid base width
    with pytest.raises(InvalidParameterError):
        simulate_npn_bjt(w_base_m=-1e-6)
