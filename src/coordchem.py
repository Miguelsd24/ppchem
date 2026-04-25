import sys
import re
from rdkit import Chem
import json
from pathlib import Path

# Loading the data from the metals/ligands.json file

BASE_DIR = Path(__file__).resolve().parent.parent

with open(BASE_DIR / "data" / "metals.json") as f:
    data_metals = json.load(f)

with open(BASE_DIR / "data" / "ligands.json") as f:
    data_ligands = json.load(f)

# The first part of this file focuses on the formula processing, along with the calculations (e.g. electron counting, oxidation state of the metal, ...)

# Function which extract a list of the metal/s from a given formula. It also verifies if the metal/s is/are in the json database
def extract_metal(formula):

    # 1. The metal string is isolated possibly along with the stoechiometric coefficient. Also, we verify that the input str has the appropriate format
    start = formula.find("[")
    end = formula.find("(")

    if start == -1 or end == -1:
        raise ValueError("Wrong input format.")

    clean_formula = formula[start + 1:end]
    match = re.match(r"^([A-Z][A-Za-z]*)([1-9]\d*)?$", clean_formula) # We use re.match() to find the metal/coefficient pattern

    if clean_formula == "" or match == None:
        raise ValueError("Wrong input format.")

    # 2. We test if the metal is present in the database and treat the output
    metal = match.group(1)
    if metal not in data_metals:
        raise ValueError(f"Metal {metal} is not in the database.")

    # 3. We test if the coefficient has the appropriate value or format and treat the output
    num = match.group(2)
    if num == None:
        coeff = 1
    elif int(num) > 2:
        raise ValueError("Invalid coefficient, the coordination compound has more than two metal centers.")
    else:
        coeff = int(num)

    return metal, coeff


# Function which extract a list of the ligand/s from a given raw formula. It also verifies if the ligand/s is/are in the json database
def extract_ligands(formula):

    # 1. The ligands string is isolated
    clean_formula = formula.replace("[", "").replace("]", "")
    parenthesis_index = clean_formula.find("(")   
    ligands = clean_formula[parenthesis_index:]
    ligands = re.findall(r"\((.*?)\)(\d*)", clean_formula)

    # 2. For each ligand in the ligands list, we isolate the stoechiometric coefficient and test if the ligands are in the database.
    result = []
    for ligand, n in ligands:
        if ligand not in data_ligands:
            sys.exit(f"Error: Ligand {ligand} not in the database")
        n = int(n) if n else 1
        result.extend([ligand] * n)
    return result

# Fonction which transforms a charge string (i.e. "2+, +2, -2, 2-") to a negative or positive int.
def parse_charge(charge):

    if charge is None:
        return 0

    charge = charge.strip()

    # Cas déjà simple
    try:
        return int(charge)
    except ValueError:
        pass

    # Cas avec signe à la fin (ex: 2+, 1-)
    match = re.match(r"(\d+)([+-])", charge)
    if match:
        value = int(match.group(1))
        sign = match.group(2)
        return value if sign == "+" else -value

    # Cas avec signe au début (ex: +3, -2)
    match = re.match(r"([+-])(\d+)", charge)
    if match:
        sign = match.group(1)
        value = int(match.group(2))
        return value if sign == "+" else -value

    raise ValueError(f"Format de charge invalide: {charge}")

# 
def complexe_charge(formula):
    charge = re.search(r"\]([0-9+-]+)", formula)
    if not charge:
        return 0
    return parse_charge(charge.group(1))

# 
def extract_elements(formula):
    elements = []
    elements.extend(extract_ligands(formula))
    elements.extend(extract_metal(formula))
    return elements

# 
def ligands_charge(formula):
    ligands = extract_ligands(formula)    
    charge = 0
    for ligand in ligands:
        charge += data_ligands[ligand]["charge"]
    return charge

# 
def oxidation_state(formula):
    charge = ligands_charge(formula)
    ox_state = complexe_charge(formula)
    metals = extract_metal(formula)
    for metal in metals:
        ox_state += data_metals[metal]["group"]
    ox_state += charge
    return ox_state

#
def electron_count(formula):
    electrons = oxidation_state(formula)
    ligands = extract_ligands(formula)
    for ligand in ligands:
        electrons += data_ligands[ligand]["denticity"]* 2
    return electrons

#
def electronic_structure(formula):
    metals = extract_metal(formula)
    per = 0
    list = []

    for metal in metals:
        per += data_metals[metal]["period"]
    inert_gas = {4: "Ar", 5: "Kr", 6: "Xe"}

    if ligands_charge(formula) + complexe_charge(formula) < 3:
        s = 2 - complexe_charge(formula)
        d = oxidation_state(formula) - s
    else:
        s = 0
        d = oxidation_state(formula) - 2

    list.extend([inert_gas.get(per), per, s, d])
    return list


def analyse_complexe(formula):
    list = electronic_structure(formula)
    print(f"Metal electronic structure : [{list[0]}] {list[1]}s{list[2]} {list[1]-1}d{list[3]}")
    print(f"Electron count : {electron_count(formula)}")

