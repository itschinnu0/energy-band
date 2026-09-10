"""Matplotlib-specific energy band visualization.

Strictly consumes SimulationResult. Does NOT calculate device physics.
Complies with all AGENTS.md rules:
- Displays EC, EV, Ei.
- Displays only physically available Fermi levels (EF, EFn, EFp).
- Never fabricates fake flat EF when EF is None.
- Never connects EFn/EFp across depletion regions when they are NaN.
- Annotates region boundaries (depletion, neutral, oxide interface).
- Supports export to PNG bytes and CSV data.
"""

from typing import Optional, Tuple, Dict, Any, List
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from physics.results import SimulationResult


def determine_spatial_scale(position_m: np.ndarray) -> Tuple[np.ndarray, str, float]:
    """Determine best spatial scale (nm or um) for x-axis.

    Args:
        position_m: Spatial coordinate array in meters.

    Returns:
        (scaled_position, unit_str, scale_factor)
    """
    total_span = float(np.ptp(position_m))
    if total_span < 1.0e-6:
        # Nanometers
        return position_m * 1.0e9, "nm", 1.0e9
    else:
        # Micrometers
        return position_m * 1.0e6, "µm", 1.0e6


def plot_band_diagram(
    result: SimulationResult,
    figsize: Tuple[float, float] = (9.0, 5.2),
    dpi: int = 100,
    show_grid: bool = True,
    show_regions: bool = True,
) -> Figure:
    """Plot energy band diagram from a SimulationResult.

    Args:
        result: Validated SimulationResult object.
        figsize: Figure size in inches (width, height).
        dpi: Dots per inch.
        show_grid: Whether to display a light grid.
        show_regions: Whether to shade / mark device region boundaries.

    Returns:
        Matplotlib Figure object.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    x_m = result.position
    x_scaled, x_unit, scale_factor = determine_spatial_scale(x_m)

    # Plot fundamental band edges
    ax.plot(x_scaled, result.EC, label=r"$E_C$ (Conduction band)", color="#1f77b4", lw=2.2)
    ax.plot(x_scaled, result.EV, label=r"$E_V$ (Valence band)", color="#d62728", lw=2.2)
    ax.plot(x_scaled, result.Ei, label=r"$E_i$ (Intrinsic level)", color="#7f7f7f", ls="--", lw=1.2, alpha=0.85)

    # Plot Fermi / Quasi-Fermi levels ONLY where physically present
    if result.EF is not None:
        ax.plot(x_scaled, result.EF, label=r"$E_F$ (Fermi level)", color="#2ca02c", ls="-.", lw=1.8)

    if result.EFn is not None:
        # Use masked array or plot with nan to break lines across depletion
        ax.plot(x_scaled, result.EFn, label=r"$E_{Fn}$ (Electron quasi-Fermi)", color="#ff7f0e", ls="-.", lw=1.8)

    if result.EFp is not None:
        ax.plot(x_scaled, result.EFp, label=r"$E_{Fp}$ (Hole quasi-Fermi)", color="#9467bd", ls=":", lw=2.0)

    # Annotate region boundaries and shade depletion/interfaces if region array exists
    if show_regions and result.region is not None and len(result.region) == len(x_scaled):
        _annotate_regions(ax, x_scaled, result.region, x_unit)

    ax.set_xlabel(f"Position $x$ [{x_unit}]", fontsize=11, fontweight="medium")
    ax.set_ylabel("Energy [eV]", fontsize=11, fontweight="medium")

    title_str = f"{result.device}: {result.operating_condition}" if result.operating_condition else result.device
    ax.set_title(title_str, fontsize=12, fontweight="bold", pad=12)

    if show_grid:
        ax.grid(True, linestyle=":", alpha=0.6)

    ax.legend(loc="best", framealpha=0.9, fontsize=9.5)

    # Tight layout and clean margins
    fig.tight_layout()
    return fig


def _annotate_regions(ax: Axes, x_scaled: np.ndarray, region: np.ndarray, x_unit: str) -> None:
    """Helper to add subtle background shading and boundary dividers for device regions."""
    # Find boundary transition indices where region label changes
    changes = np.where(region[:-1] != region[1:])[0]

    unique_regions = []
    reg_starts = [0] + [c + 1 for c in changes]
    reg_ends = [c for c in changes] + [len(region) - 1]

    colors = {
        "p-depletion": "#fff2cc",
        "n-depletion": "#fff2cc",
        "depletion": "#fff2cc",
        "eb-depletion-emitter": "#fee8d6",
        "eb-depletion-base": "#fee8d6",
        "cb-depletion-base": "#e5f5e0",
        "cb-depletion-collector": "#e5f5e0",
        "semiconductor-depletion": "#f2f0f7",
    }

    for start_idx, end_idx in zip(reg_starts, reg_ends):
        lbl = str(region[start_idx])
        x0 = x_scaled[start_idx]
        x1 = x_scaled[end_idx]

        # Depletion shading
        for key, bg_col in colors.items():
            if key in lbl.lower():
                ax.axvspan(x0, x1, color=bg_col, alpha=0.45, zorder=0)
                break

    for c in changes:
        x_trans = (x_scaled[c] + x_scaled[c + 1]) / 2.0
        ax.axvline(x_trans, color="#888888", linestyle="--", linewidth=0.8, alpha=0.6, zorder=1)


def export_result_to_csv_dataframe(result: SimulationResult) -> pd.DataFrame:
    """Convert SimulationResult spatial profile into a clean pandas DataFrame for CSV export.

    Follows Section 48 export specification:
        position_m
        EC_eV
        EV_eV
        Ei_eV
        EF_eV
        EFn_eV
        EFp_eV
        potential_V
        electric_field_V_per_m
        region
    Unavailable quantities are represented as NaN / blank.
    """
    data = {
        "position_m": result.position,
        "EC_eV": result.EC,
        "EV_eV": result.EV,
        "Ei_eV": result.Ei,
    }

    data["EF_eV"] = result.EF if result.EF is not None else np.nan
    data["EFn_eV"] = result.EFn if result.EFn is not None else np.nan
    data["EFp_eV"] = result.EFp if result.EFp is not None else np.nan

    if result.potential is not None:
        data["potential_V"] = result.potential
    else:
        data["potential_V"] = np.nan

    if result.electric_field is not None:
        data["electric_field_V_per_m"] = result.electric_field
    else:
        data["electric_field_V_per_m"] = np.nan

    if result.region is not None:
        data["region"] = result.region
    else:
        data["region"] = ""

    df = pd.DataFrame(data)
    return df


def export_figure_to_png_bytes(fig: Figure, dpi: int = 150) -> bytes:
    """Export Matplotlib figure to PNG byte array for download."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()
