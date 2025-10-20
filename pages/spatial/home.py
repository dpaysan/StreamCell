import hydralit_components as hc
import streamlit as st

from pages.spatial.data import spatial_data_page
from pages.spatial.embeddings import spatdata_embedding_page
from pages.spatial.features import spatdata_feature_page
from pages.spatial.spatial import spatdata_spatial_page

# NavBar

DAT = 'Data set'
FEATS = 'Feature Plots'
EMBS = 'Embedding Plots'
SPAT = 'Spatial Plots'

tabs = [
    DAT,
    FEATS,
    EMBS,
]

option_data = [
    {'icon': "🧬", 'label': DAT},
    {'icon': "🧬", 'label': FEATS},
    {'icon': "🔎", 'label': EMBS},
    {'icon': "🔎", 'label': SPAT},

]


over_theme = {'txc_inactive': 'black', 'menu_background': '#D6E5FA', 'txc_active': 'white', 'option_active': '#749BC2'}
font_fmt = {'font-class': 'h3', 'font-size': '50%'}


def spatial_home_page():
    st.header('Spatial data')
    page = hc.option_bar(
        option_definition=option_data,
        title='',
        override_theme=over_theme,
        horizontal_orientation=True)

    if page == DAT:
        spatial_data_page()
    elif page == FEATS:
        spatdata_feature_page()
    elif page == EMBS:
        spatdata_embedding_page()
    elif page == SPAT:
        spatdata_spatial_page()
    else:
        pass