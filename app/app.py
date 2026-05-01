import streamlit as st
import sys
import os
import time

# 1. Trouve le dossier où se trouve app.py
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Remonte d'un niveau et va dans 'src'
src_path = os.path.join(current_dir, "..", "src")

# 3. Ajoute-le au système
if src_path not in sys.path:
    sys.path.append(src_path)

import coordchem as cc


# Title of the app
st.title("CoordChemPy", text_alignment="left")
st.subheader("The best python based tool for coordination chemist!", divider="gray",  text_alignment="left")
# 
st.subheader("Coordination compound information finder")

coord_compound  = st.text_input("_Enter the coordination compound formula following the correct format_", placeholder="Type here ...", help="Input format rules are explained in the README.md file")

if st.button("Analysis"):
    if coord_compound:
        try:
            result = cc.show_analysis(coord_compound)
            st.success("Successful analysis")
            st.text(result)
        except ValueError as e:
            st.error(str(e))
    else:
        st.warning("Please enter a formula.")

st.divider()