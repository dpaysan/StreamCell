import streamlit as st

# Inspired by TFinder App from Juiliet Minnet

import hydralit_components as hc
import platform
import pandas as pd
import requests
import streamlit as st
from streamlit_modal import Modal
import streamlit_lottie
import time
import json

from pages.other.contact import contact_page
from pages.other.home import home_page
from pages.scrnaseq.home import scdata_home_page
from pages.spatial.home import spatial_home_page
from utils.layout import footer_style, footer

try:
    from streamlit import rerun as rerun
except ImportError:
    # conditional import for streamlit version <1.27
    from streamlit import experimental_rerun as rerun

import os


st.set_page_config(
    page_title='StreamCell',
    page_icon="",
    initial_sidebar_state="expanded",
    layout="wide"
)

max_width_str = f"max-width: {75}%;"

st.markdown(f"""
        <style>
        .appview-container .main .block-container{{{max_width_str}}}
        </style>
        """,
            unsafe_allow_html=True,
            )

st.markdown("""
        <style>
               .block-container {
                    padding-top: 0rem;
                    padding-bottom: 0rem;

                }
        </style>
        """, unsafe_allow_html=True)

# Footer

st.markdown(footer_style, unsafe_allow_html=True)

st.title("StreamCell")

from streamlit_option_menu import option_menu

with st.sidebar:
    selected = option_menu(
        menu_title="Main Menu",
        options=["Home", "scRNA", "Spatial", "Downloads","Contact"],
        icons=["house", "file-bar-graph", "file-image", "cloud-download", "envelope"],
        default_index=0,
        menu_icon="menu-app-fill"
    )


if selected == "Home":
    home_page()

elif selected == "scRNA":
    scdata_home_page()

elif selected == "Spatial":
    spatial_home_page()

elif selected == "Contact":
    contact_page()

for i in range(4):
    st.markdown('#')
st.markdown(footer, unsafe_allow_html=True)



# Credit

st.sidebar.title("Credit")
# Help
st.sidebar.title("Help")
st.sidebar.markdown("Documentation", unsafe_allow_html=True)


st.sidebar.title("More")
