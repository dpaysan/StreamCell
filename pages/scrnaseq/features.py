from utils.helpers import *

def scdata_feature_page():

    st.set_page_config(page_title="Feature display", page_icon="🎯", layout="wide")

    # ---- dataset ----
    scdata = st.session_state.get("scdata", None)
    if scdata is None:
        st.error("Error: No dataset found.")
        st.stop()

    # ---- figsize defaults (persist across reruns) ----
    if "figsize_map" not in st.session_state:
        st.session_state["figsize_map"] = {"Violin": (12.0, 7.0), "Heatmap": (10.0, 6.0), "Dot plot": (12.0, 7.0)}
    for k in ("Violin", "Heatmap", "Dot plot"):
        st.session_state["figsize_map"].setdefault(k, (12.0, 7.0))

    # =========================
    #       TOP CONTROLS
    # =========================
    cat_cols = get_categorical_obs(scdata.obs)
    num_cols = get_numeric_obs(scdata.obs)

    st.subheader("Configuration")
    c1, c2, c3, c4 = st.columns([1.1, 1.1, 1.1, 1.1])
    with c1:
        plot_type = st.radio("Plot type", ["Violin", "Heatmap", "Dot plot"], horizontal=True)
    with c2:
        feature_source = st.selectbox("Feature source", ["Genes", "Metadata"])
    with c3:
        if cat_cols:
            groupby_sel = st.selectbox("Group by (obs column)", options=["(none)"] + cat_cols, index=0)
            groupby = None if groupby_sel == "(none)" else groupby_sel
        else:
            groupby = None
            st.write("**Group by**: (no categorical columns found — using a single group)")
    with c4:
        layer_sel = st.selectbox("Data layer", options=["(none)"] + list(scdata.layers), index=0)
        layer = None if layer_sel == "(none)" else layer_sel

    # Feature input
    features = []
    if feature_source == "Genes":
        features_text = st.text_input(
            "Genes (comma-separated var names), e.g. SERPINA7, KRT7, MARCO",
            placeholder="GeneA, GeneB, GeneC, ...",
        )
        features, missing = parse_features(features_text, scdata.var_names)
        if features_text and missing:
            st.info("Ignoring not-found genes: " + ", ".join(missing[:10]) + ("…" if len(missing) > 10 else ""))
    else:
        # Only allow obs numerics for Violin & Dot plot
        if plot_type == "Heatmap":
            st.warning("Heatmap supports gene features only. Switch Feature source to 'Gene (var)'.")
            features = []
        else:
            features = st.multiselect("Obs numeric columns", options=num_cols, placeholder="Select numeric obs columns")

    # Subset controls (optional)
    subset_available = bool(cat_cols)
    if subset_available:
        st.markdown("**Subset (optional)**")
        sub_col = st.selectbox("Filter column", options=["(none)"] + cat_cols, index=0)
        subset_values = []
        if sub_col != "(none)":
            cats = pd.Series(scdata.obs[sub_col].astype("category")).cat.categories
            all_values = list(cats) if len(cats) > 0 else sorted(scdata.obs[sub_col].dropna().unique().tolist())
            subset_values = st.multiselect("Keep only these values", options=all_values, default=all_values)
    else:
        sub_col, subset_values = "(none)", []

    # ---- Apply subsetting ----
    scdata_view = scdata
    if subset_available and sub_col != "(none)":
        if len(subset_values) == 0:
            scdata_view = scdata[[]]
        else:
            mask = scdata.obs[sub_col].isin(subset_values)
            scdata_view = scdata[mask]

    # ---- Status strip ----
    st.divider()
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("n_obs (filtered)", f"{scdata_view.n_obs:,}")
    m2.metric("n_vars", f"{scdata_view.n_vars:,}")
    m3.metric("Group by", groupby if groupby else "all")
    m4.metric("Subset", sub_col if (subset_available and sub_col != "(none)") else "None")

    # ---- Plot trigger (single, large plot, no column layout around it) ----

    st.divider()
    curr_w, curr_h = st.session_state["figsize_map"][plot_type]

    generate = st.button("🎨 Generate plot", type="primary")
    if generate:
        if scdata_view.n_obs == 0:
            st.warning("No cells to display. Adjust your subset or remove the filter.")
        elif not features:
            st.warning("Enter at least one feature to plot.")
        else:
            render_plot(scdata_view, plot_type, groupby=groupby, features=features, fig_w=curr_w, fig_h=curr_h, layer=layer)

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

