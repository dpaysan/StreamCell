import streamlit as st

def home_page():


    st.set_page_config(page_title="Single-cell & Spatial Visualizer", page_icon="🧬", layout="wide")

    st.header("🧬 Single-cell & Spatial Visualizer")
    st.markdown("""
    Welcome to the **Single-cell & Spatial Visualizer**, a web-based interactive application designed to help you explore single-cell (snRNA-seq) and spatial transcriptomics datasets with ease.

    ### 🚀 What this app provides
    - A unified interface for **single-nucleus RNA-seq** and **spatial transcriptomics** data.
    - Multiple visualization modes: violin plots, dot plots, heatmaps, embeddings (PCA/UMAP), spatial mapping.
    - Built-in support for **two datasets** in shared cache, enabling fast access and minimal memory overhead.
    - Extensible architecture: you can add new dataset keys, new plot types, new UI pages — easily and modularly.

    ### 🔧 Why we built this
    Bioinformatics tools are often scattered across disparate scripts and formats.  
    With this app, our aim is to bring together standardized loading, caching, and visualization in one place — so researchers can focus on biological insight rather than boilerplate code.

    ### 🙏 Acknowledgements
    This design was **inspired by** the open-source application TFFinder (by Julien Minet on GitHub), whose elegant UI and modular architecture served as a blueprint for this work.

    ### 📝 Cite Our Work
    If you use this tool in your research, please cite our accompanying paper (see sidebar for details).  
    We hope you find it useful, and welcome suggestions & issues via GitHub.

    ---

    *Start by selecting a dataset from the sidebar or navigating through the pages above.*  
    """)
