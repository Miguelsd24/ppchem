import streamlit as st
import sys
import os

# 1. Trouve le dossier où se trouve app.py
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Remonte d'un niveau et va dans 'src'
src_path = os.path.join(current_dir, "..", "src")

# 3. Ajoute-le au système
if src_path not in sys.path:
    sys.path.append(src_path)

# 4. Maintenant, l'import fonctionnera !
import coordchem as cc



# Title of the app
st.title("CoordChemPy")
st.markdown("The best python based tool for coordination chemist!")
# 
st.subheader("Coordination compound Information finder")
st.text("...")

coord_compound  = st.text_input("Enter the coordination compound formula following the correct format", "Type here...")

#
if st.button("Submit"):
    result = cc.analyse_complexe(coord_compound)
    st.success(result)



