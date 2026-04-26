# CoordChemPy

![PyPI](https://img.shields.io/pypi/v/coordchempy)
![License](https://img.shields.io/github/license/Miguelsd24/ppchem)

## About CoordChemPy

CoordChemPy is a Python package designed to assist inorganic chemists and chemistry students by providing tools for analysis and modeling of coordination compounds. This package includes :

- Calculation about coordination compounds like electron counting, metal electronic structure
- ...

Some approximations and assumptions in order to yield a fully fonctional chemistry package :
- The ligand database is not exhaustive
- Only classical trasition metals were considered. Lanthanides, actinides and heavy synthetic metals (Rf -> Cn) were excluded
- Coordination complexes with more than two metal centers are not incorporated
- Heterobinuclear complexes are not incorporated and homobinuclear complexes must be symmetric with respect to the two metals

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install CoordChemPy.

```bash
pip install coordchempy
```

## Utilisation
### Importation
```python
import coordchempy
...
```
### Formulas input rules
One of the main inputs of this package is the formula of the coordination complexes. To ensure compatibility with the code, the formula must follow specific rules:
- The coordination sphere must be indicated by brackets
- The coordination sphere charge must be indicated outside the brackets at the end of the formula. Nothing corresponds to a neutral charge and both notations; ±x / x± where x is the charge are equivalent
- The metal come first, inside the brackets, with its stoechiometric coefficient right after. As mentionned before this package contrains this coefficient to 1 or 2, any other values will lead to errors.
- The ligands come after the metals and are indicated inside parenthesis, the stoechiometric coefficient comes outside the parenthesis after the ligand
- Bridging ligand are announced by a prefix "m-" which lies inside the parethesis

For instance, some valid formulas are :
```python
"[V(Cl)4]"
"[Ti(H2O)5(CO)]2+"
"[Co2(m-NH2)(m-OH)(NH3)8]4+"
```
### Coordination compound compational data functions

## Contributing
Contributors are welcome to suggest improvements at https://github.com/Miguelsd24/coordchempy
