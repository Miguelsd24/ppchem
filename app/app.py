import streamlit as st
import streamlit.components.v1 as components
import sys
import os
import py3Dmol
from ase.io import write
import io

# 1. Trouve le dossier où se trouve app.py
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Remonte d'un niveau et va dans 'src'
src_path = os.path.join(current_dir, "..", "src")

# 3. Ajoute-le au système
if src_path not in sys.path:
    sys.path.append(src_path)

import coordchem as cc


def render_molecule(compound):
    xyz_str = io.StringIO()
    write(xyz_str, compound, format='xyz')
    xyz_content = xyz_str.getvalue()

    view = py3Dmol.view(width=400, height=400)
    view.addModel(xyz_content, 'xyz')

    if render_type == "Ball and Stick":
        view.setStyle({'stick': {}, 'sphere': {'scale': atoms_size}})
    elif render_type == "Stick":
        view.setStyle({'stick': {}})
    elif render_type == "Sphere":
        view.setStyle({'sphere': {'scale': atoms_size}})
    elif render_type == "Lines":
        view.setStyle({'line': {}})
    elif render_type == "VDW":  
        view.addSurface(py3Dmol.VDW)
    view.zoomTo()
    
    view_html = view._make_html()
    components.html(view_html, height=400, width=400)



# Title of the app
st.title("CoordChemPy", text_alignment="left")
st.subheader("The best python based tool for coordination chemist! :atom_symbol:", divider="gray",  text_alignment="left")
# 
st.subheader("Coordination compound information finder")

coord_compound  = st.text_input("_Enter the coordination compound formula following the correct format. Try [Pt(Cl)2(NH3)2] !_", placeholder="Type here ...", help="Input format rules are explained in the README.md file")


# Session state
# ----------------------------------------------------------------------------------------
if "analysis_result" not in st.session_state:
    st.session_state.analysis_result = None

if "compound_ase" not in st.session_state:
    st.session_state.compound_ase = None


# ----------------------------------------------------------------------------------------
if st.button("Analysis"):
    if coord_compound:
        try:
            st.session_state.analysis_result = cc.analyze_complexe(coord_compound)[0]
        except ValueError as e:
            st.error(str(e))
            st.session_state.analysis_result = None
            st.session_state.compound_ase = None
    else:
        st.warning("Please enter a formula.")
        st.session_state.analysis_result = None
        st.session_state.compound_ase = None

if st.session_state.analysis_result:
    for line in st.session_state.analysis_result:
        st.write(line)
    st.success("Successful analysis")


st.divider()
# ----------------------------------------------------------------------------------------


st.subheader("Coordination compound 3D rendering")

size_range = [i/10 for i in range(1, 11)]
render_options = ["Ball and Stick", "Stick", "Sphere", "Lines", "VDW"]
with st.form("render_form"):
    atoms_size = st.select_slider("Atom size", options=size_range, value=0.3)
    render_type = st.segmented_control("Rendring type", options=render_options, selection_mode="single", default = "Ball and Stick", help="VDW: Van der Waals" )
    submit = st.form_submit_button("3D render")

if submit:
    if coord_compound:
        try:
            st.session_state.compound_ase = cc.create_compound_render(coord_compound)
        except ValueError as e:
            st.error(str(e))
            st.session_state.analysis_result = None
            st.session_state.compound_ase = None

    else:
        st.warning("Please enter a formula.")
        st.session_state.analysis_result = None
        st.session_state.compound_ase = None

if st.session_state.compound_ase:
    render_molecule(st.session_state.compound_ase)
    st.success("Successful render")


st.divider()



