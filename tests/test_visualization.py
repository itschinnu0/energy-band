"""Tests for band diagram visualization, exports, and application integration."""

import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from physics.pn_junction import simulate_pn_junction
from physics.bjt import simulate_npn_bjt
from physics.mos import simulate_mos
from visualization.band_plot import (
    plot_band_diagram,
    export_result_to_csv_dataframe,
    export_figure_to_png_bytes,
)


def test_plot_pn_equilibrium_and_biased():
    """Verify plotting a valid PN result in equilibrium and bias does not fail."""
    # Equilibrium
    res_eq = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.0)
    fig_eq = plot_band_diagram(res_eq)
    assert fig_eq is not None
    plt.close(fig_eq)

    # Biased
    res_fwd = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.3)
    fig_fwd = plot_band_diagram(res_fwd)
    assert fig_fwd is not None
    plt.close(fig_fwd)


def test_plot_bjt_all_modes():
    """Verify plotting BJT results across different operating modes does not fail."""
    modes = [
        (0.0, 0.0),    # Equilibrium
        (0.65, -1.5),  # Forward active
        (-1.0, -2.0),  # Cutoff
        (0.6, 0.4),    # Saturation
    ]
    for vbe, vbc in modes:
        res = simulate_npn_bjt(vbe_v=vbe, vbc_v=vbc)
        fig = plot_band_diagram(res)
        assert fig is not None
        plt.close(fig)


def test_plot_mos_levels():
    """Verify plotting MOS results for both Level 1 and Level 2 does not fail."""
    # Level 1
    res_l1 = simulate_mos(vg_v=1.5, level=1)
    fig_l1 = plot_band_diagram(res_l1)
    assert fig_l1 is not None
    plt.close(fig_l1)

    # Level 2
    res_l2 = simulate_mos(vg_v=-1.5, level=2)
    fig_l2 = plot_band_diagram(res_l2)
    assert fig_l2 is not None
    plt.close(fig_l2)


def test_unavailable_fermi_levels_not_fabricated():
    """Verify that when EF is None, no line is drawn for it, and NaN regions in EFn/EFp are respected."""
    res_biased = simulate_pn_junction(na_cm3=1e16, nd_cm3=1e16, va_v=0.4)
    assert res_biased.EF is None
    assert res_biased.EFn is not None
    assert np.any(np.isnan(res_biased.EFn))

    fig = plot_band_diagram(res_biased)
    ax = fig.axes[0]
    labels = [line.get_label() for line in ax.get_lines()]

    # Should have EC, EV, Ei, EFn, EFp, but NOT EF
    assert r"$E_C$ (Conduction band)" in labels
    assert r"$E_V$ (Valence band)" in labels
    assert r"$E_i$ (Intrinsic level)" in labels
    assert r"$E_{Fn}$ (Electron quasi-Fermi)" in labels
    assert r"$E_{Fp}$ (Hole quasi-Fermi)" in labels
    assert r"$E_F$ (Fermi level)" not in labels
    plt.close(fig)


def test_export_png_bytes():
    """Verify figure export to PNG produces non-empty PNG bytes."""
    res = simulate_pn_junction()
    fig = plot_band_diagram(res)
    png_bytes = export_figure_to_png_bytes(fig)
    assert isinstance(png_bytes, bytes)
    assert len(png_bytes) > 1000
    # Check PNG magic bytes: \x89PNG\r\n\x1a\n
    assert png_bytes[:8] == b"\x89PNG\r\n\x1a\n"
    plt.close(fig)


def test_export_csv_dataframe():
    """Verify spatial profile export to DataFrame conforms to Section 48 specification."""
    res = simulate_pn_junction(va_v=0.2)
    df = export_result_to_csv_dataframe(res)
    assert isinstance(df, pd.DataFrame)
    expected_cols = [
        "position_m",
        "EC_eV",
        "EV_eV",
        "Ei_eV",
        "EF_eV",
        "EFn_eV",
        "EFp_eV",
        "potential_V",
        "electric_field_V_per_m",
        "region",
    ]
    assert list(df.columns) == expected_cols
    assert len(df) == len(res.position)
    # Unavailable EF should be NaN
    assert df["EF_eV"].isna().all()
    # EC - EV should equal Eg
    np.testing.assert_allclose(df["EC_eV"] - df["EV_eV"], res.calculated_parameters["bandgap_eV"], atol=1e-9)


def test_app_model_invocations():
    """Verify all three models can be invoked directly from application layer parameters."""
    pn = simulate_pn_junction(na_cm3=1e17, nd_cm3=1e16, va_v=-1.0)
    assert pn.device == "PN Junction"
    assert pn.model_valid is True

    bjt = simulate_npn_bjt(vbe_v=0.6, vbc_v=-1.0)
    assert bjt.device == "NPN BJT"
    assert bjt.model_valid is True

    mos = simulate_mos(vg_v=0.5, level=1)
    assert mos.device == "NMOS"
    assert mos.model_valid is True
