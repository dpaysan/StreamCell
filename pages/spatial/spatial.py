import streamlit as st
from utils.helpers import *

# Persist per-embedding default sizes (like your feature page)
if "figsize_map" not in st.session_state:
    st.session_state["figsize_map"] = {
        "spatial": (18.0, 9.0),
    }
else:
    st.session_state["figsize_map"]["spatial"] = (18.0, 9.0)




def spatdata_spatial_page():

    st.set_page_config(page_title="Spatial display", page_icon="🧭", layout="wide")

    # ---- Get cached AnnData loaded on the Load page ----
    spatdata = st.session_state.get("spatdata", None)
    if spatdata is None:
        st.error("Error: No dataset found.")
        st.stop()


    # =========================
    #       TOP CONTROLS
    # =========================
    cat_cols = get_categorical_obs(spatdata.obs)

    samples = spatdata.obs.sample_id.unique
    if not samples:
        st.error("No samples found.")
        st.stop()


    st.subheader("Configuration")

    c0, c1 = st.columns([1.2, 1.2])
    with c0:
        selected_sample = st.selectbox("Sample", list(spatdata.obs.sample_id.unique()))
    with c1:
        color_mode = st.selectbox("Color by", ["Metadata column", "Gene"])


    c2, c3 = st.columns([1.2, 1.2])
    with c2:
        # Color controls
        obs_color_col = "(none)"
        gene_text = ""
        if color_mode == "Metadata column":
            obs_color_col = st.selectbox("Metadata column (obs)", options=["(none)"] + list(spatdata.obs.columns), index=0)
        else:
            gene_text = st.text_input("Gene (var name)", placeholder="e.g. MALAT1")


    st.divider()
    st.markdown("**Graphic controls**")

    c4,c5,c6, c7 = st.columns([1.2, 1.2, 1.2, 1.2])
    with c4:
        pt_size = st.slider("Point size", 1.0, 10.0, 1.0, 1.0)
    with c5:
        alpha = st.slider("Opacity", 0.0, 1.0, 1.0, 0.1)
    with c6:
            vmin = st.text_input("Min. value", placeholder=0, value=0)
    with c7:
            vmax = st.text_input("Max. value", placeholder="p95", value="p95")

    # ---- Apply subsetting (only if categorical filter exists) ----
    mask = spatdata.obs.sample_id.isin([selected_sample]).values
    spatdata_view = spatdata[mask]

    # Optional: keep only spots with gene > 0 when coloring by a gene
    only_pos = False
    if color_mode == "Gene":
        only_pos = st.checkbox("Show only spots with selected gene expression > 0", value=True)
        if only_pos and gene_text.strip():
            gmask = gene_positive_mask(spatdata_view, gene_text.strip(), layer=None)  # or pass your `layer`
            if gmask.any():
                spatdata_view = spatdata_view[gmask]
            else:
                st.info("No spots with expression > 0 for this gene in the current sample.")


    # ---- Status strip ----
    st.divider()
    m1, m2, m3 = st.columns(3)
    m1.metric("n_obs (filtered)", f"{spatdata_view.n_obs:,}")
    m2.metric("n_vars", f"{spatdata_view.n_vars:,}")
    m3.metric("Sample", selected_sample)

    plot_type = "spatial"



    # ---- Plot trigger (single, large plot, no column layout around it) ----

    st.divider()
    curr_w, curr_h = st.session_state["figsize_map"][plot_type]
    color_arg = resolve_color(spatdata, color_mode, obs_color_col, gene_text)

    # helpful feedback if user selected Gene mode but the gene is missing
    if color_mode == "Gene" and gene_text.strip() and color_arg is None:
        st.info(
            f"Gene '{gene_text.strip()}' not found in `adata.var_names` (case-insensitive check). Showing uncolored points.")

    if spatdata_view.n_obs == 0:
        st.warning("No cells to display. Adjust your subset or remove the filter.")
    else:
        render_plot(spatdata_view, plot_type, color_arg=color_arg, size=pt_size, fig_w=curr_w, fig_h=curr_h,
                    alpha=alpha, vmin=vmin, vmax=vmax, tissue = selected_sample)

        # ---- Resize controls BELOW the plot ----
        st.subheader("Resize plot")
        # Use dedicated keys so different plot types don't step on each other
        w_key = f"fig_w_{plot_type}"
        h_key = f"fig_h_{plot_type}"

        # Initialize slider state if first time on this plot type
        if w_key not in st.session_state:
            st.session_state[w_key] = curr_w
        if h_key not in st.session_state:
            st.session_state[h_key] = curr_h

        def _apply_resize():
            # Update the per-plot-type map and re-run to draw with new size above
            st.session_state["figsize_map"][plot_type] = (
                st.session_state[w_key],
                st.session_state[h_key],
            )
            st.rerun()

        st.slider("Width (inches)", 4.0, 20.0, st.session_state[w_key], 0.5, key=w_key, on_change=_apply_resize)
        st.slider("Height (inches)", 3.0, 16.0, st.session_state[h_key], 0.5, key=h_key, on_change=_apply_resize)