# appdata/datasets.py
import os

import numpy as np
import pandas as pd
import streamlit as st
import scanpy as sc
from scipy import sparse

# --- Hard-coded dataset paths (adjust to your server) ---
DATASETS = {
    "scdata": "data/snrnaseq/snrnaseq_example.h5ad",
    "spatdata": "data/visiumhd/visiumhd_example.h5ad",
}

def _path_ok(path: str) -> bool:
    try:
        return os.path.exists(path)
    except Exception:
        return False

def _make_readonly(x):
    """Mark ndarray or sparse matrix as read-only."""
    if x is None:
        return x
    if sparse.issparse(x):
        try:
            x.data.flags.writeable = False
        except Exception:
            pass
        return x
    if isinstance(x, np.ndarray):
        try:
            x.setflags(write=False)
        except Exception:
            pass
    return x

def _freeze_df(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame that raises if someone tries to modify it."""
    df = df.copy(deep=True)

    class FrozenDataFrame(pd.DataFrame):
        _metadata = []

        def __setitem__(self, key, value):
            raise AttributeError("This DataFrame is read-only")

        def drop(self, *args, **kwargs):
            raise AttributeError("This DataFrame is read-only")

        def assign(self, *args, **kwargs):
            raise AttributeError("This DataFrame is read-only")

    frozen = FrozenDataFrame(df)
    return frozen

@st.cache_resource(show_spinner=False)
def _load_dataset(key: str):
    path = DATASETS[key]
    adata = sc.read_h5ad(path)

    # Freeze numerical data
    adata.X = _make_readonly(adata.X)
    if hasattr(adata, "layers"):
        for k, v in list(adata.layers.items()):
            adata.layers[k] = _make_readonly(v)

    # Freeze obs/var dataframes
    adata.obs = _freeze_df(adata.obs)
    adata.var = _freeze_df(adata.var)

    # Lock attribute setting (e.g. adata.uns["foo"] = ...)
    _orig_setattr = adata.__setattr__

    def _locked_setattr(self, name, value):
        if name not in {"_locked", "__setattr__"} and getattr(self, "_locked", False):
            raise AttributeError(f"AnnData object is locked (read-only); cannot set attribute {name}")
        _orig_setattr(self, name, value)

    adata._locked = True
    adata.__setattr__ = _locked_setattr.__get__(adata, type(adata))
    return adata


def get_adata(key: str):
    """Return the shared AnnData for the given key (loads once per process)."""
    if key not in DATASETS:
        raise KeyError(f"Unknown dataset key: {key}. Known: {list(DATASETS)}")
    return _load_dataset(key)

def available_datasets():
    """List available dataset keys."""
    return list(DATASETS.keys())
