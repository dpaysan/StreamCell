import streamlit as st

from appdata.datasets import get_adata
from utils.helpers import fmt_bytes, adata_memory_bytes


def scdata_data_page():

    st.set_page_config(page_title="Load single-cell data", page_icon="📥")

    st.subheader("📥 Load single-cell dataset")

    with st.status("Initializing dataset cache ...", expanded=True) as status:
        st.write("Loading single cell data...")
        try:
            ad = get_adata("scdata")  # first call triggers actual load; subsequent calls are instant
        except Exception as e:
            st.error(f"  Failed to load single-cell data: {e}")

    status.update(label="Dataset ready in shared cache ✅", state="complete")

    st.success("Single-cell dataset loaded and ready for use across the app.")

    # Keep convenient references (no copies) for current session if you like:
    st.session_state[f"scdata"] = ad

    # --- Summaries ---
    st.subheader("Data set")
    c1, c2, c3 = st.columns(3)
    c1.metric("Observations (n_obs)", f"{ad.n_obs:,}")
    c2.metric("Genes (n_vars)", f"{ad.n_vars:,}")
    c3.metric("Approx. memory", fmt_bytes(adata_memory_bytes(ad)))