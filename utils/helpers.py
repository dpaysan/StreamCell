import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import scanpy as sc
from typing import List, Tuple, Any

from scipy import sparse

if "figsize_map" not in st.session_state:
    st.session_state["figsize_map"] = {
        "Violin": (14.0, 9.0),
        "Heatmap": (8.0, 6.0),
        "Dot plot": (13.0, 8.0),
    }


def available_embeddings(adata) -> List[str]:
    # Canonical in AnnData: .obsm["X_umap"], .obsm["X_pca"]
    obsm = getattr(adata, "obsm", {})
    # Be forgiving if someone wrote "OpsM" by mistake
    if not hasattr(adata, "obsm") and hasattr(adata, "OpsM"):
        obsm = getattr(adata, "OpsM")
    opts = []
    if "X_umap" in obsm: opts.append("UMAP")
    if "X_pca" in obsm: opts.append("PCA")
    return opts


def resolve_color(adata, mode: str, obs_col: str, gene_text: str):
    """Return a valid Scanpy 'color' arg or None; handles gene name case-insensitively."""
    if mode == "Metadata column":
        return None if obs_col == "(none)" else obs_col
    # Gene mode
    gene = (gene_text or "").strip()
    if not gene:
        return None
    if gene in adata.var_names:
        return gene
    # case-insensitive match
    hits = [g for g in adata.var_names if g.lower() == gene.lower()]
    return hits[0] if hits else None


# ---- Helpers ----
def get_categorical_obs(df: pd.DataFrame, max_unique: int = 50) -> List[str]:
    cats = []
    for c in df.columns:
        if pd.api.types.is_categorical_dtype(df[c]) or df[c].dtype == "object":
            cats.append(c)
        else:
            try:
                if df[c].nunique(dropna=True) <= max_unique:
                    cats.append(c)
            except Exception:
                pass
    return sorted(cats)


def parse_features(text: str, universe: pd.Index) -> Tuple[List[str], List[str]]:
    feats = [t.strip() for t in text.split(",") if t.strip()] if text else []
    present = [f for f in feats if f in universe]
    missing = [f for f in feats if f not in universe]
    return present, missing


def ensure_groupby_column(adata_view, groupby: str) -> str:
    """Guarantee a valid groupby column; if none provided, create '__all__'."""
    if groupby:
        return groupby
    if "__all__" not in adata_view.obs.columns:
        adata_view.obs["__all__"] = "all"
    return "__all__"


def render_plot(adata_view, plot_type, groupby=None, features=None, fig_w=None, fig_h=None,
                color_arg=None, size=None,
                alpha=None, vmin=None, vmax=None, alpha_img=1.0,
                layer=None, tissue = None):
    with st.spinner("Plotting data — please wait…"):
        plt.close("all")
        gcol = ensure_groupby_column(adata_view, groupby)

        if plot_type == "Violin":
            sc.pl.violin(
                adata_view,
                keys=features,
                groupby=gcol,
                stripplot=False,
                jitter=0.3,
                rotation=90,
                show=False,
                layer=layer,
            )
            fig = plt.gcf()

        elif plot_type == "Heatmap":
            sc.pl.heatmap(
                adata_view,
                var_names=features,
                groupby=gcol,
                use_raw=None,
                standard_scale="var",
                cmap="rocket_r",
                figsize=(fig_w, fig_h),  # direct control
                show=False,
                layer=layer,
            )
            fig = plt.gcf()

        elif plot_type == "Dot plot":  # Dot plot
            dp = sc.pl.dotplot(
                adata_view,
                var_names=features,
                groupby=gcol,
                figsize=(fig_w, fig_h),  # direct control
                show=False,
                layer=layer,
            )
            fig = getattr(dp, "fig", None) or getattr(dp, "_fig", None) or plt.gcf()

        elif plot_type == "UMAP":
            sc.pl.umap(
                adata_view,
                color=color_arg,
                size=size,
                show=False,
                vmin=vmin,
                vmax=vmax,
                alpha=alpha,
                cmap="rocket_r",
                layer=layer,
            )
            fig = plt.gcf()
            fig.set_size_inches(fig_w, fig_h)
        elif plot_type == "PCA":
            sc.pl.pca(
                adata_view,
                color=color_arg,
                size=size,
                show=False,
                vmin=vmin,
                vmax=vmax,
                alpha=alpha,
                cmap="rocket_r",
                layer=layer,
            )
            fig = plt.gcf()
            fig.set_size_inches(fig_w, fig_h)
        elif plot_type == "spatial":
            fig = sc.pl.spatial(adata_view,
                                color=color_arg,
                                ncols=1,
                                na_color=None,
                                library_id=tissue,
                                vmin=vmin,
                                vmax=vmax,
                                return_fig=True,
                                size=size,
                                show=False,
                                cmap="rocket_r",
                                marker = "s",
                                alpha_img = alpha_img
                                )
            fig.set_size_inches(fig_w, fig_h)
        else:
            st.error("Unknown plot type: {}".format(plot_type))

        st.pyplot(fig, clear_figure=True, width="content")
        plt.close(fig)


def fmt_bytes(num: int) -> str:
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if num < 1024 or unit == "TB":
            return f"{num:.2f} {unit}"
        num /= 1024.0


def _array_nbytes(x: Any) -> int:
    """Safely get approximate memory usage (in bytes) of arrays/sparse arrays."""
    if x is None:
        return 0
    if sparse.issparse(x):
        # count data + index arrays
        total = 0
        for attr in ("data", "indices", "indptr"):
            a = getattr(x, attr, None)
            if a is not None and hasattr(a, "nbytes"):
                total += a.nbytes
        return total
    if hasattr(x, "nbytes"):
        return int(x.nbytes)
    return 0


def adata_memory_bytes(adata) -> int:
    """Approximate total bytes used by key AnnData components."""
    total = 0
    # core matrix
    total += _array_nbytes(adata.X)

    # obs/var dataframes
    if isinstance(adata.obs, pd.DataFrame):
        total += int(adata.obs.memory_usage(deep=True).sum())
    if isinstance(adata.var, pd.DataFrame):
        total += int(adata.var.memory_usage(deep=True).sum())

    # layers
    if isinstance(adata.layers, dict):
        for v in adata.layers.values():
            total += _array_nbytes(v)

    # obsm/varm (dict of arrays)
    for mapping in (getattr(adata, "obsm", {}), getattr(adata, "varm", {})):
        if isinstance(mapping, dict):
            for v in mapping.values():
                total += _array_nbytes(v)

    # obsp/varp (pairwise/sparse matrices)
    for mapping in (getattr(adata, "obsp", {}), getattr(adata, "varp", {})):
        if hasattr(mapping, "items"):
            for v in mapping.values():
                total += _array_nbytes(v)

    return total


def show_summary(adata):
    mem_bytes = adata_memory_bytes(adata)
    st.subheader("Dataset summary")
    col1, col2, col3 = st.columns(3)
    col1.metric("Observations (n_obs)", f"{adata.n_obs:,}")
    col2.metric("Genes (n_vars)", f"{adata.n_vars:,}")
    col3.metric("Approx. memory", fmt_bytes(mem_bytes))


import numpy as np
from scipy import sparse

def gene_positive_mask(adata, gene: str, layer: str | None = None, threshold: float = 0.0) -> np.ndarray:
    """Return boolean mask for cells with expression > threshold for `gene`."""
    if gene not in adata.var_names:
        # if the gene isn't present, don't drop anything
        return np.ones(adata.n_obs, dtype=bool)

    gi = int(np.where(adata.var_names == gene)[0][0])
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    col = X[:, gi]  # (n_obs, 1)

    if sparse.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col > float(threshold)

def resolve_gene_list(candidates: List[str], universe: pd.Index) -> Tuple[List[str], List[str]]:
    present = [g for g in candidates if g in universe]
    missing = [g for g in candidates if g not in universe]
    return present, missing

def get_numeric_obs(df: pd.DataFrame) -> List[str]:
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    return sorted(num)
