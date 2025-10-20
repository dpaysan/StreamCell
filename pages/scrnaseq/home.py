import streamlit as st
import hydralit_components as hc

from pages.scrnaseq.data import scdata_data_page
from pages.scrnaseq.embeddings import scdata_embedding_page
from pages.scrnaseq.features import scdata_feature_page

# NavBar

DAT = 'Dataset'
FEATS = 'Feature Plots'
EMBS = 'Embedding Plots'

tabs = [
    DAT,
    FEATS,
    EMBS,
]

option_data = [
    {'icon': "🧬", 'label': DAT},
    {'icon': "🧬", 'label': FEATS},
    {'icon': "🔎", 'label': EMBS},

]


over_theme = {'txc_inactive': 'black', 'menu_background': '#D6E5FA', 'txc_active': 'white', 'option_active': '#749BC2'}
font_fmt = {'font-class': 'h3', 'font-size': '50%'}

def scdata_home_page():
    st.header('Single-cell data')
    page = hc.option_bar(
        option_definition=option_data,
        title='',
        override_theme=over_theme,
        horizontal_orientation=True)

    if page == DAT:
        scdata_data_page()
    elif page == FEATS:
        scdata_feature_page()
    elif page == EMBS:
        scdata_embedding_page()
    else:
        pass


