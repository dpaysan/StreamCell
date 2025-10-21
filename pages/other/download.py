import streamlit as st
import os

def show_file_download(label, path):
    if not os.path.exists(path):
        st.error(f"File not found: `{path}`")
        return

    size = os.path.getsize(path)
    units = ["B","KB","MB","GB"]
    for u in units:
        if size < 1024 or u == "GB": break
        size /= 1024
    fmt = f"{size:.2f} {u}"

    with open(path, "rb") as f:
        st.download_button(
            label=f"Download {label} ({fmt})",
            data=f,
            file_name=os.path.basename(path),
            mime="application/octet-stream",
            use_container_width=False,
        )

def downloads_page():
    import streamlit as st

    st.set_page_config(page_title="Download datasets", page_icon="⬇️", layout="wide")
    st.title("⬇️ Download datasets")

    # Store your file-share URLs in Streamlit secrets or env vars
    # e.g. in .streamlit/secrets.toml:
    # [files]
    # scdata_url = "https://files.example.org/data/scdata.h5ad"
    # spatdata_url = "https://files.example.org/data/spatdata.h5ad"
    FILES = st.secrets.get("files", {})

    sc_url = FILES.get("scdata_url", "")
    sp_url = FILES.get("spatdata_url", "")

    def link_block(title: str, url: str):
        st.subheader(title)
        if url:
            # Opens in a new tab; the file is fetched directly from your file server/CDN
            st.link_button(f"Download {title} (as .h5ad)", url)
        else:
            st.warning("No URL configured. Add it to `st.secrets['files']`.")

    c0, c1 = st.columns(2)
    with c0:
        link_block("Single-cell data", sc_url)
    with c1:
        link_block("Spatial data", sp_url)

    st.info("The above downloads the datasets as h5ad files from our fileshare")