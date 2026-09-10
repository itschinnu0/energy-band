"""Streamlit GUI for Energy Band Diagram Simulator.

Integrates:
- 1-D Abrupt PN Junction
- NPN BJT (Two-Junction Electrostatic Model)
- NMOS / MOS Structure (Level-1 & Level-2)
- Band diagram visualization
- Voltage sweeps and sweep plots
- CSV & PNG data export
- Clean physics error handling without silent clipping
"""

from typing import Dict, Any, List
import io
import math
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from physics.exceptions import (
    EnergyBandError,
    InvalidParameterError,
    ModelValidityError,
    SolverConvergenceError,
    UnitConversionError,
)
from physics.pn_junction import simulate_pn_junction
from physics.bjt import simulate_npn_bjt
from physics.mos import simulate_mos
from visualization.band_plot import (
    plot_band_diagram,
    export_result_to_csv_dataframe,
    export_figure_to_png_bytes,
)


st.set_page_config(
    page_title="Energy Band Diagram Simulator",
    page_icon="⚡",
    layout="wide",
)


def format_param_table(params: Dict[str, Any]) -> pd.DataFrame:
    """Format dictionary of calculated parameters into a clean display DataFrame with units."""
    formatted = []
    unit_map = {
        "Vbi_V": ("Built-in potential (Vbi)", "V"),
        "VD_V": ("Depletion barrier (VD)", "V"),
        "VA_V": ("Applied junction voltage (VA)", "V"),
        "W_m": ("Total depletion width (W)", "µm", 1e6),
        "xp_m": ("p-side depletion width (xp)", "µm", 1e6),
        "xn_m": ("n-side depletion width (xn)", "µm", 1e6),
        "peak_electric_field_V_per_m": ("Peak electric field", "kV/cm", 1e-5),
        "potential_drop_V": ("Potential drop across depletion", "V"),
        "temperature_K": ("Temperature", "K"),
        "bandgap_eV": ("Bandgap (Eg)", "eV"),
        "NA_cm3": ("Acceptor doping (NA)", "cm⁻³"),
        "ND_cm3": ("Donor doping (ND)", "cm⁻³"),
        "NE_cm3": ("Emitter doping (NE)", "cm⁻³"),
        "NB_cm3": ("Base doping (NB)", "cm⁻³"),
        "NC_cm3": ("Collector doping (NC)", "cm⁻³"),
        "VBE_V": ("Base-Emitter voltage (VBE)", "V"),
        "VBC_V": ("Base-Collector voltage (VBC)", "V"),
        "Vbi_EB_V": ("EB built-in potential", "V"),
        "Vbi_CB_V": ("CB built-in potential", "V"),
        "Vbarrier_EB_V": ("EB junction barrier", "V"),
        "Vbarrier_CB_V": ("CB junction barrier", "V"),
        "W_EB_m": ("EB depletion width", "µm", 1e6),
        "W_CB_m": ("CB depletion width", "µm", 1e6),
        "x_e_m": ("Emitter depletion width (xe)", "nm", 1e9),
        "x_be_m": ("Base EB depletion width (xbe)", "nm", 1e9),
        "x_bc_m": ("Base CB depletion width (xbc)", "nm", 1e9),
        "x_c_m": ("Collector depletion width (xc)", "µm", 1e6),
        "W_base_m": ("Metallurgical base width (WB)", "µm", 1e6),
        "neutral_base_width_m": ("Neutral base width", "µm", 1e6),
        "operating_region": ("Operating region", ""),
        "Cox_F_per_m2": ("Oxide capacitance (Cox)", "µF/cm²", 1e2),
        "phi_F_V": ("Fermi potential (phi_F)", "V"),
        "VFB_V": ("Flat-band voltage (VFB)", "V"),
        "threshold_voltage_V": ("Threshold voltage (VT)", "V"),
        "psi_s_V": ("Surface potential (psi_s)", "V"),
        "Qs_C_per_m2": ("Semiconductor charge (Qs)", "µC/cm²", 1e2),
        "Wd_m": ("Depletion width (Wd)", "nm", 1e9),
        "Wd_max_m": ("Max depletion width (Wd,max)", "nm", 1e9),
        "VG_V": ("Gate voltage (VG)", "V"),
        "operating_state": ("Operating state", ""),
        "tox_m": ("Oxide thickness (tox)", "nm", 1e9),
    }

    for k, v in params.items():
        if k in ("solver_diagnostics", "level"):
            continue
        if k in unit_map:
            info = unit_map[k]
            label = info[0]
            unit = info[1]
            if len(info) == 3 and isinstance(v, (int, float)):
                scaled_val = v * info[2]
                val_str = f"{scaled_val:.4g}"
            elif isinstance(v, float):
                val_str = f"{v:.4g}"
            else:
                val_str = str(v)
            formatted.append({"Parameter": label, "Value": val_str, "Unit": unit})
        else:
            formatted.append({"Parameter": k, "Value": str(v), "Unit": ""})

    return pd.DataFrame(formatted)


def main():
    st.title("Semiconductor Energy Band Diagram Simulator")
    st.markdown(
        "A physics-first analytical 1-D energy-band simulator for **PN Junctions**, **NPN BJTs**, "
        "and **MOS Structures**."
    )

    device = st.sidebar.selectbox(
        "Select Device",
        ["1-D PN Junction", "NPN BJT", "NMOS / MOS Capacitor"],
    )

    st.sidebar.markdown("---")
    st.sidebar.header("Device Parameters")

    if device == "1-D PN Junction":
        render_pn_gui()
    elif device == "NPN BJT":
        render_bjt_gui()
    else:
        render_mos_gui()


def render_pn_gui():
    col1, col2 = st.sidebar.columns(2)
    na_exp = col1.number_input("Acceptor NA [10^x cm⁻³]", min_value=14.0, max_value=19.0, value=16.0, step=0.5)
    nd_exp = col2.number_input("Donor ND [10^x cm⁻³]", min_value=14.0, max_value=19.0, value=16.0, step=0.5)
    na = 10.0**na_exp
    nd = 10.0**nd_exp

    temp_k = st.sidebar.slider("Temperature [K]", min_value=100.0, max_value=500.0, value=300.0, step=10.0)
    va = st.sidebar.slider("Applied Voltage VA [V]", min_value=-5.0, max_value=0.75, value=0.0, step=0.05)

    st.sidebar.caption("Positive VA = Forward bias; Negative VA = Reverse bias.")

    st.sidebar.markdown("---")
    st.sidebar.header("Voltage Sweep (Optional)")
    enable_sweep = st.sidebar.checkbox("Enable VA Sweep", value=False)
    if enable_sweep:
        sw_c1, sw_c2 = st.sidebar.columns(2)
        va_start = sw_c1.number_input("VA Start [V]", value=-4.0, step=0.5)
        va_stop = sw_c2.number_input("VA Stop [V]", value=0.5, step=0.1)
        va_steps = st.sidebar.slider("Points", min_value=10, max_value=100, value=30)

    try:
        res = simulate_pn_junction(na_cm3=na, nd_cm3=nd, va_v=va, temp_k=temp_k)

        # Header Info
        st.subheader(f"Condition: {res.operating_condition}")
        if res.warnings:
            for w in res.warnings:
                st.warning(w)

        col_plot, col_data = st.columns([1.6, 1.0])

        with col_plot:
            fig = plot_band_diagram(res)
            st.pyplot(fig)

            # PNG Export button
            png_bytes = export_figure_to_png_bytes(fig)
            st.download_button(
                label="📥 Download Band Diagram (PNG)",
                data=png_bytes,
                file_name=f"pn_band_VA_{va:+.2f}V.png",
                mime="image/png",
            )

        with col_data:
            st.markdown("#### Calculated Parameters")
            df_params = format_param_table(res.calculated_parameters)
            st.dataframe(df_params, use_container_width=True, hide_index=True)

            # CSV Export button
            df_profile = export_result_to_csv_dataframe(res)
            csv_bytes = df_profile.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 Download Profile Data (CSV)",
                data=csv_bytes,
                file_name=f"pn_profile_VA_{va:+.2f}V.csv",
                mime="text/csv",
            )

        if enable_sweep:
            st.markdown("---")
            st.subheader("Voltage Sweep Analysis")
            va_vals = np.linspace(va_start, va_stop, va_steps)
            sw_data = []

            for v in va_vals:
                try:
                    r_sw = simulate_pn_junction(na_cm3=na, nd_cm3=nd, va_v=v, temp_k=temp_k)
                    p = r_sw.calculated_parameters
                    sw_data.append({
                        "VA [V]": v,
                        "VD [V]": p["VD_V"],
                        "W [µm]": p["W_m"] * 1e6,
                        "xp [µm]": p["xp_m"] * 1e6,
                        "xn [µm]": p["xn_m"] * 1e6,
                        "Peak E-field [kV/cm]": p["peak_electric_field_V_per_m"] * 1e-5,
                    })
                except ModelValidityError:
                    continue

            if sw_data:
                df_sw = pd.DataFrame(sw_data)
                fig_sw, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.0))

                ax1.plot(df_sw["VA [V]"], df_sw["W [µm]"], color="#1f77b4", lw=2, marker=".")
                ax1.set_xlabel("Applied Voltage VA [V]")
                ax1.set_ylabel("Depletion Width W [µm]")
                ax1.set_title("Depletion Width vs Applied Voltage")
                ax1.grid(True, linestyle=":", alpha=0.6)

                ax2.plot(df_sw["VA [V]"], df_sw["Peak E-field [kV/cm]"], color="#d62728", lw=2, marker=".")
                ax2.set_xlabel("Applied Voltage VA [V]")
                ax2.set_ylabel("Peak E-Field [kV/cm]")
                ax2.set_title("Peak Electric Field vs Applied Voltage")
                ax2.grid(True, linestyle=":", alpha=0.6)

                fig_sw.tight_layout()
                st.pyplot(fig_sw)

                st.download_button(
                    label="📥 Download Sweep Data (CSV)",
                    data=df_sw.to_csv(index=False).encode("utf-8"),
                    file_name="pn_voltage_sweep.csv",
                    mime="text/csv",
                )

    except (InvalidParameterError, ModelValidityError, UnitConversionError) as e:
        st.error(f"❌ Physical Validation Error: {e}")


def render_bjt_gui():
    col1, col2, col3 = st.sidebar.columns(3)
    ne_exp = col1.number_input("NE [10^x]", min_value=16.0, max_value=20.0, value=18.0, step=0.5)
    nb_exp = col2.number_input("NB [10^x]", min_value=14.0, max_value=18.0, value=16.0, step=0.5)
    nc_exp = col3.number_input("NC [10^x]", min_value=14.0, max_value=17.0, value=15.0, step=0.5)
    ne = 10.0**ne_exp
    nb = 10.0**nb_exp
    nc = 10.0**nc_exp

    w_base_um = st.sidebar.slider("Base Width [µm]", min_value=0.2, max_value=3.0, value=1.0, step=0.1)
    temp_k = st.sidebar.slider("Temperature [K]", min_value=100.0, max_value=500.0, value=300.0, step=10.0)

    vbe = st.sidebar.slider("Base-Emitter VBE [V]", min_value=-3.0, max_value=0.75, value=0.6, step=0.05)
    vbc = st.sidebar.slider("Base-Collector VBC [V]", min_value=-10.0, max_value=0.75, value=-1.5, step=0.25)

    st.sidebar.caption("VBE = VB - VE; VBC = VB - VC.")

    st.sidebar.markdown("---")
    st.sidebar.header("VBE Sweep (Optional)")
    enable_sweep = st.sidebar.checkbox("Enable VBE Sweep", value=False)
    if enable_sweep:
        sw_c1, sw_c2 = st.sidebar.columns(2)
        vbe_start = sw_c1.number_input("VBE Start [V]", value=-2.0, step=0.2)
        vbe_stop = sw_c2.number_input("VBE Stop [V]", value=0.65, step=0.05)
        vbe_steps = st.sidebar.slider("Points", min_value=10, max_value=60, value=25)

    try:
        res = simulate_npn_bjt(
            ne_cm3=ne,
            nb_cm3=nb,
            nc_cm3=nc,
            vbe_v=vbe,
            vbc_v=vbc,
            temp_k=temp_k,
            w_base_m=w_base_um * 1e-6,
        )

        st.subheader(f"Operating Mode: {res.operating_condition}")
        if res.warnings:
            for w in res.warnings:
                st.warning(w)

        col_plot, col_data = st.columns([1.6, 1.0])

        with col_plot:
            fig = plot_band_diagram(res)
            st.pyplot(fig)

            png_bytes = export_figure_to_png_bytes(fig)
            st.download_button(
                label="📥 Download Band Diagram (PNG)",
                data=png_bytes,
                file_name=f"bjt_band_VBE_{vbe:+.2f}V_VBC_{vbc:+.2f}V.png",
                mime="image/png",
            )

        with col_data:
            st.markdown("#### Calculated Parameters")
            df_params = format_param_table(res.calculated_parameters)
            st.dataframe(df_params, use_container_width=True, hide_index=True)

            df_profile = export_result_to_csv_dataframe(res)
            csv_bytes = df_profile.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 Download Profile Data (CSV)",
                data=csv_bytes,
                file_name=f"bjt_profile_VBE_{vbe:+.2f}V_VBC_{vbc:+.2f}V.csv",
                mime="text/csv",
            )

        if enable_sweep:
            st.markdown("---")
            st.subheader("BJT Sweep Analysis")
            vbe_vals = np.linspace(vbe_start, vbe_stop, vbe_steps)
            sw_data = []

            for vb_pt in vbe_vals:
                try:
                    r_sw = simulate_npn_bjt(
                        ne_cm3=ne,
                        nb_cm3=nb,
                        nc_cm3=nc,
                        vbe_v=vb_pt,
                        vbc_v=vbc,
                        temp_k=temp_k,
                        w_base_m=w_base_um * 1e-6,
                    )
                    p = r_sw.calculated_parameters
                    sw_data.append({
                        "VBE [V]": vb_pt,
                        "Region": p["operating_region"],
                        "EB Barrier [V]": p["Vbarrier_EB_V"],
                        "Neutral Base Width [µm]": p["neutral_base_width_m"] * 1e6,
                        "W_EB [nm]": p["W_EB_m"] * 1e9,
                    })
                except ModelValidityError:
                    continue

            if sw_data:
                df_sw = pd.DataFrame(sw_data)
                fig_sw, ax = plt.subplots(figsize=(8, 3.8))
                ax.plot(df_sw["VBE [V]"], df_sw["Neutral Base Width [µm]"], color="#2ca02c", lw=2, marker="o")
                ax.set_xlabel("VBE [V]")
                ax.set_ylabel("Neutral Base Width [µm]")
                ax.set_title(f"Neutral Base Width vs VBE (at VBC = {vbc:.2f} V)")
                ax.grid(True, linestyle=":", alpha=0.6)
                fig_sw.tight_layout()
                st.pyplot(fig_sw)

                st.download_button(
                    label="📥 Download BJT Sweep CSV",
                    data=df_sw.to_csv(index=False).encode("utf-8"),
                    file_name="bjt_vbe_sweep.csv",
                    mime="text/csv",
                )

    except (InvalidParameterError, ModelValidityError, UnitConversionError) as e:
        st.error(f"❌ Physical Validation Error: {e}")


def render_mos_gui():
    na_exp = st.sidebar.number_input("Substrate NA [10^x cm⁻³]", min_value=14.0, max_value=18.0, value=16.0, step=0.5)
    na = 10.0**na_exp
    tox_nm = st.sidebar.slider("Oxide Thickness tox [nm]", min_value=2.0, max_value=50.0, value=10.0, step=1.0)
    phims = st.sidebar.number_input("Work Function Diff PhiMS [V]", value=-0.95, step=0.05)
    qox_c_m2 = st.sidebar.number_input("Oxide Charge Qox [10⁻⁴ C/m²]", value=0.0, step=0.1) * 1e-4

    temp_k = st.sidebar.slider("Temperature [K]", min_value=100.0, max_value=500.0, value=300.0, step=10.0)
    vg = st.sidebar.slider("Gate Voltage VG [V]", min_value=-4.0, max_value=4.0, value=1.0, step=0.1)

    level_choice = st.sidebar.radio("Model Level", ["Level 1 (Analytical Depletion Approx)", "Level 2 (Numerical Poisson-Boltzmann)"])
    level = 1 if "Level 1" in level_choice else 2

    st.sidebar.markdown("---")
    st.sidebar.header("Gate Voltage Sweep (Optional)")
    enable_sweep = st.sidebar.checkbox("Enable VG Sweep", value=False)
    if enable_sweep:
        sw_c1, sw_c2 = st.sidebar.columns(2)
        vg_start = sw_c1.number_input("VG Start [V]", value=-3.0, step=0.5)
        vg_stop = sw_c2.number_input("VG Stop [V]", value=3.0, step=0.5)
        vg_steps = st.sidebar.slider("Points", min_value=10, max_value=80, value=30)

    try:
        res = simulate_mos(
            na_cm3=na,
            tox_m=tox_nm * 1e-9,
            vg_v=vg,
            phims_v=phims,
            qox_c_m2=qox_c_m2,
            temp_k=temp_k,
            level=level,
        )

        st.subheader(f"MOS Condition: {res.operating_condition}")
        if res.warnings:
            for w in res.warnings:
                st.warning(w)

        col_plot, col_data = st.columns([1.6, 1.0])

        with col_plot:
            fig = plot_band_diagram(res)
            st.pyplot(fig)

            png_bytes = export_figure_to_png_bytes(fig)
            st.download_button(
                label="📥 Download Band Diagram (PNG)",
                data=png_bytes,
                file_name=f"mos_band_VG_{vg:+.2f}V_L{level}.png",
                mime="image/png",
            )

        with col_data:
            st.markdown("#### Calculated Parameters")
            df_params = format_param_table(res.calculated_parameters)
            st.dataframe(df_params, use_container_width=True, hide_index=True)

            if "solver_diagnostics" in res.calculated_parameters:
                diag = res.calculated_parameters["solver_diagnostics"]
                st.caption(
                    f"Solver Diagnostics: converged={diag['converged']}, "
                    f"residual={diag['residual']:.2e}, iterations={diag['iterations']}"
                )

            df_profile = export_result_to_csv_dataframe(res)
            csv_bytes = df_profile.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 Download Profile Data (CSV)",
                data=csv_bytes,
                file_name=f"mos_profile_VG_{vg:+.2f}V_L{level}.csv",
                mime="text/csv",
            )

        if enable_sweep:
            st.markdown("---")
            st.subheader("MOS Voltage Sweep Analysis")
            vg_vals = np.linspace(vg_start, vg_stop, vg_steps)
            sw_data = []

            for v in vg_vals:
                try:
                    r_sw = simulate_mos(
                        na_cm3=na,
                        tox_m=tox_nm * 1e-9,
                        vg_v=v,
                        phims_v=phims,
                        qox_c_m2=qox_c_m2,
                        temp_k=temp_k,
                        level=level,
                    )
                    p = r_sw.calculated_parameters
                    sw_data.append({
                        "VG [V]": v,
                        "Surface Potential psi_s [V]": p["psi_s_V"],
                        "Qs [µC/cm²]": p["Qs_C_per_m2"] * 1e2,
                        "Wd [nm]": p["Wd_m"] * 1e9,
                        "Operating State": p["operating_state"],
                    })
                except (ModelValidityError, SolverConvergenceError):
                    continue

            if sw_data:
                df_sw = pd.DataFrame(sw_data)
                fig_sw, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.0))

                ax1.plot(df_sw["VG [V]"], df_sw["Surface Potential psi_s [V]"], color="#1f77b4", lw=2, marker=".")
                ax1.set_xlabel("Gate Voltage VG [V]")
                ax1.set_ylabel("Surface Potential ψs [V]")
                ax1.set_title("ψs vs VG")
                ax1.grid(True, linestyle=":", alpha=0.6)

                ax2.plot(df_sw["VG [V]"], df_sw["Qs [µC/cm²]"], color="#e377c2", lw=2, marker=".")
                ax2.set_xlabel("Gate Voltage VG [V]")
                ax2.set_ylabel("Semiconductor Charge Qs [µC/cm²]")
                ax2.set_title("Qs vs VG")
                ax2.grid(True, linestyle=":", alpha=0.6)

                fig_sw.tight_layout()
                st.pyplot(fig_sw)

                st.download_button(
                    label="📥 Download MOS Sweep CSV",
                    data=df_sw.to_csv(index=False).encode("utf-8"),
                    file_name=f"mos_sweep_L{level}.csv",
                    mime="text/csv",
                )

    except (InvalidParameterError, ModelValidityError, SolverConvergenceError, UnitConversionError) as e:
        st.error(f"❌ Simulation Error: {e}")


if __name__ == "__main__":
    main()
