import sys
import os
import numpy as np
import re
import json
from pathlib import Path
import roman
import string
from IPython.display import display, Markdown

# Allows to display only one line when running into an error in the jupyter notebook
def hide_traceback(exc_tuple=None, filename=None, tb_offset=None, exception_only=False, running_compiled_code=False):
    etype, value, tb = sys.exc_info()
    return print(f"{value}")

try:
    ipython = get_ipython()
except NameError:
    pass # This replaces the default Jupyter error handler with our clean one-liner

# Loading the data from the metals/ligands.json file

BASE_DIR = Path(__file__).resolve().parent.parent

with open(BASE_DIR / "data" / "metals.json") as f:
    data_metals = json.load(f)

with open(BASE_DIR / "data" / "ligands.json") as f:
    data_ligands = json.load(f)

# The first part of this file focuses on the formula format processing, along with the calculations (e.g. electron counting, oxidation state of the metal, ...)



##################################################################### TO DO <
def fomula_input_verification(formula):
    type_input = type(formula)
    if type_input != str:
        raise ValueError("Error: Formula must be a string")
    return True
##################################################################### TO DO <



# Function which process the formula to have a better looking LaTeX formula
def clean_formula(formula):
    clean_formula = re.sub(r"m-", r"μ-", formula)
    start = clean_formula.find("[")
    end = clean_formula.find("]")
    clean_formula_1 = (clean_formula[:start] + re.sub(r"(\d+)", r"_{\1}", clean_formula[start:end + 1]))
    signe = ""
    if complexe_charge(formula) > 0 :
        signe = "+"
    if complexe_charge(formula) == 0 :
        return "$" + clean_formula_1 + "$"
    clean_formula_2 = "^{" +  signe + str(complexe_charge(formula)) + "}"
    return "$" + clean_formula_1 + clean_formula_2 + "$"

# Function which extract the metal and its stoechiometric coefficient from a given formula. It also verifies if the metal is in the json database
def parse_metal(formula):

    # 1. The metal string is isolated possibly along with the stoechiometric coefficient. Also, we verify that the input str has the appropriate format
    start = formula.find("[")
    end = formula.find("(")
    if start == -1 or end == -1  or start >= end:
        raise ValueError("Error: Wrong input format.")
    clean_formula = formula[start + 1:end]

    # 2. We use re.match() to find the metal/coefficient pattern
    match = re.match(r"^(?:(s|d|t))?([A-Z][A-Za-z]*)([1-9]\d*)?$", clean_formula)
    if match == None:
        raise ValueError("Error: Wrong input format.")

    # 3. We test if the metal is present in the database and treat the output
    metal = match.group(2)
    if metal not in data_metals:
        raise ValueError(f"Error: Metal {metal} is not in the database.")

    # 4. We test if the coefficient has the appropriate value or format and treat the output
    num = match.group(3)
    coeff = test_string_for_int(num)
    if coeff > 2:
        raise ValueError("Error: Invalid coefficient, the coordination compound has more than two metal centers.")
    
    # 5. We put the metal and its coefficient in a list and return it
    metals = []
    metals.extend([metal]*coeff)
    return metals

# Function which extracts form the formula if the metals are bonded (single, double ,triple) or not
def bond_order(formula):
    start = formula.find("[")
    end = formula.find("(")
    clean_formula = formula[start + 1:end]
    match = re.match(r"^(?:(s|d|t))?([A-Z][A-Za-z]*)([1-9]\d*)?$", clean_formula) 
    order_dico = {"s": 1, "d":2, "t":3}
    order = order_dico.get(match.group(1), 0)
    if len(parse_metal(formula)) == 1 and order != 0:
        raise ValueError("Error: s,d,t are only to specifiy the bond between two metals center not one")
    return order

# Function which find ligands by name or abbr in the database
def find_ligand(ligand_input):
    # 1. Direct verification
    if ligand_input in data_ligands:
        return ligand_input
    # 2. Verification by abbreviation
    for ligand_key, properties in data_ligands.items():
        if "abbr" in properties and properties["abbr"] == ligand_input:
            return ligand_key
    # 3. Verification for a bridging ligand
    if ligand_input.startswith("m-"):
        for ligand_key in data_ligands:
            if ligand_key == ligand_input[2:]:
                return ligand_input
    return False

# Function which extract a list of the ligand/s from a given raw formula. It also verifies if the ligand/s is/are in the json database
def parse_ligands(formula):

    # 1. The ligand string is isolated possibly along with the stoechiometric coefficient. Also, we verify that the input str has the appropriate format
    start = formula.find("(")
    end = formula.find("]")

    if start == -1 or end == -1  or start >= end:
        raise ValueError("Error: Wrong input format.")

    # 2. Using re.findall we isolate each ligand according to the correct format
    clean_formula = formula[start:end]
    match = re.findall(r"\((.*?)\)(\d*)", clean_formula)

    if match == None:
        raise ValueError("Error: Wrong input format.")

    # 3. For each ligand in the ligands list, we isolate the stoechiometric coefficient and test if the ligands are in the database
    result = []
    coeff = []
    for ligand, n in match:
        if test_string_for_int(n) > 10:
            raise ValueError(f"Error: No ligand can have a coefficient superior to 10")
        elif ligand == "":
            raise ValueError(f"Error: At least one ligand is missing in the formula")
        elif find_ligand(ligand) == False:
            raise ValueError(f"Error: Ligand {ligand} not in the database")
        else:
            coeff.append(test_string_for_int(n))
            result.extend([find_ligand(ligand)])

    return result, coeff 
    
# Function which creates a list with all the ligands times their coeficient. This list is useful to use loops in the calculation function like that There is non need to deal with number later.
def ligands_list(formula):
    ligands_list = []
    ligands = parse_ligands(formula)[0]
    coeff = parse_ligands(formula)[1]
    for n in range(len(ligands)):
        ligands_list.extend([ligands[n]]*coeff[n])     
    return ligands_list

# Function which counts the number of bridging ligands
def count_bridging_ligands(formula):
    num = 0
    for ligand in ligands_list(formula):
        if ligand.startswith("m-"):
            num += 1
    return num

# Function which transforms a charge string (i.e. "+, -, X+, +X, -X, X-") to a negative or positive int
def parse_charge(charge):

    # 1. Case with nothing
    if charge is None:
        return 0
    charge = charge.strip()

    # 2. Case without sign
    try:
        return int(charge)
    except ValueError:
        pass

    # 3. Case only the sign
    if charge == "+":
        return 1
    if charge == "-":
        return -1

    # 4. Case with end sign (ex: 2+, 1-)
    match = re.match(r"(\d+)([+-])", charge)
    if match:
        value = int(match.group(1))
        sign = match.group(2)
        return value if sign == "+" else -value

    # 5. Case with start sign (ex: +3, -2)
    match = re.match(r"([+-])(\d+)", charge)
    if match:
        sign = match.group(1)
        value = int(match.group(2))
        return value if sign == "+" else -value

    raise ValueError(f"Error: Invalid format charge for the coordination sphere: {charge}")

# Function which test if a string is a int or float or str in order to return the good stoechiometric coefficient
def test_string_for_int(num):

    if num is None or num == "":
        return 1
    
    try:
        coeff = float(num)
    except (ValueError, TypeError):
        raise ValueError("Error: Coefficient must be numeric")

    if not coeff.is_integer():
       raise ValueError("Error: Invalid coefficient, the coefficient cannot be a float") 

    coeff = int(num)

    if coeff == 0:
        raise ValueError("Error: Invalid coefficient, coefficients cannot be zero")    
    elif coeff < 0:
        raise ValueError("Error: Invalid coefficient, coefficients cannot be negative")    
    
    return coeff

# Function which return the charge of the coordination sphere as an int
def complexe_charge(formula):
    charge = re.search(r"\]([0-9+-]+)", formula)
    if not charge:
        return 0
    return parse_charge(charge.group(1))

# Function which put the metal and ligands in a same list with their respective stoechiometric coefficient
def parse_elements(formula):
    elements = []
    elements.extend(ligands_list(formula))
    elements.extend(parse_metal(formula))
    # We also verify that there is no bridging ligand if there is only one metal center
    if len(parse_metal(formula)) == 1 and count_bridging_ligands(formula) > 0:
        raise ValueError("Error: A mononuclear complex cannot have bridging ligands")
    return elements

# Function which calulate the sum of all ligands' charges
def ligands_charge(formula):
    ligands = ligands_list(formula)
    charge = 0
    for ligand in ligands:
        if not ligand.startswith("m-"):
           charge += data_ligands[ligand]["charge"] 
        else:
            charge += data_ligands[ligand[2:]]["charge"]
    return charge

# Function which calulate the charge of the metal center (MX+ or MX-)
def metal_charge(formula):
    charge = (complexe_charge(formula) - ligands_charge(formula))//len(parse_metal(formula))
    return charge 

# Fonction which calulate the ox. state of the metal center (dX)
def oxidation_state(formula):
    metals = parse_metal(formula)

    ox_state = (data_metals[metals[0]]["group"] - metal_charge(formula))
    if ox_state < 0 or ox_state > 10:
        raise ValueError(f"Error: There is too extreme charge")
    return ox_state
    
# Function which does the electron counting
def electron_count(formula):
    metal = parse_metal(formula)[0]
    d = oxidation_state(formula)

    electrons = d

    for ligand in ligands_list(formula):
        if ligand.startswith("m-"):
            electrons += data_ligands[ligand[2:]]["bridging_e"]
        else:
            electrons += data_ligands[ligand]["donor_e"]

    electrons += 2 * bond_order(formula)

    return int(electrons)

# Function which calulate the electroni structure of the metal
def electronic_structure(formula):

    metals = parse_metal(formula)
    per = 0
    list = []

    # We first deal with the period and which noble gaz is the base of the electronic structure (both metals are the same so their period is also the same)
    per += data_metals[metals[0]]["period"]
    inert_gas = {4: "Ar", 5: "Kr", 6: "Xe"}

    # Then, we deal with the As^b Cd^e part, also, we seperate if the charge is negative/ between 0-2/ > 2, because it is the s 2 electrons that are first removed
    if metal_charge(formula) == 0 or metal_charge(formula) == 2:
        s = 2 - metal_charge(formula)
        d = oxidation_state(formula) - s
    elif metal_charge(formula) >= 3 or metal_charge(formula) == 1: # Because the only electron in the s orbital fall into the d orbital, because lower in energy
        s = 0
        d = oxidation_state(formula)
    else:
        s = 2 
        d = oxidation_state(formula) - 2 

    list.extend([inert_gas.get(per), per, s, d])
    return list


# =============================================================================================================================================================== #
        
coeff_name1 = {
    1: "",
    2: "di",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
    7: "hepta",
    8: "octa",
    9: "nona",
    10: "deca"
}
        
coeff_name2 = {
    1: "",
    2: "bis",
    3: "tris",
    4: "tetrakis",
    5: "pentakis",
    6: "hexakis",
    7: "heptakis",
    8: "octakis",
    9: "nonakis",
    10: "decakis",
}


def name_ligand(ligand):
    n = 0
    if ligand.startswith("m-"):
        n = 2
    if data_ligands[ligand[n:]].get("nomenclature") is not None: # Generally the name is used in the nomenclature but no always.
        return data_ligands[ligand[n:]]["nomenclature"] 
    else:
        return data_ligands[ligand[n:]]["name"]
    

def should_use_the_coeff_name2(ligand_name):
    for ligand in data_ligands.values():
        if ligand["name"] == ligand_name and ligand.get("coeff") == "yes":
            return True
    return False


def naming_compound(formula):
    # 1. We set the data nedded
    parsed_data = parse_ligands(formula)
    ligands = parsed_data[0]
    coeffs = parsed_data[1]
    metals = parse_metal(formula)
    ligands_with_coeffs = []
    name = ""
    mu = "-" + "\u03bc" + "-"
    last_parenthesis =  False

    # 2. We indentify the the bridging ligands by transforming their coefficient to a negative value (no chimical sense just easier to inditify later)
    for n in range(len(ligands)):
        if ligands[n].startswith("m-"):
            coeffs[n] *= -1
    # 3. Ligands as name ou abbr and are sortes by alphabetic order
        ligand_name = name_ligand(ligands[n])
        ligands_with_coeffs.append((ligand_name,coeffs[n]))
    ligands_with_coeffs.sort(key=lambda x: x[0].lower())


    # 4. Metal as name (primary or secondary)
    if complexe_charge(formula) < 0:
        metal_name = data_metals[metals[0]]["secondary_name"]
    else:
        metal_name = data_metals[metals[0]]["name"]

    # 5. We seperate the cases

    #-------------BRIDGING LIGANDS---------------#
    for i, (ligand_name, coeff) in enumerate(ligands_with_coeffs):
        if coeff < 0:
            bridging = mu
            if should_use_the_coeff_name2(ligand_name) == True:     # We seperate the case where we have to use the second type of prefixes according to the ligand (IUPAC rules)
                prefixe_ligand = coeff_name2[coeff*-1]
                name += (f"{bridging}{prefixe_ligand}({ligand_name})")
            else:
                prefixe_ligand = coeff_name1[coeff*-1]
                name += bridging + prefixe_ligand + ligand_name

    n = 1
    #------------- 2 METALS ---------------#
    
    if len(metals) == 2 and len(ligands_list(formula))-count_bridging_ligands(formula) > 0:
        name += coeff_name2[2] + "("
        n = 2
        last_parenthesis =  True
    elif len(metals) == 2 and len(ligands_list(formula))-count_bridging_ligands(formula) == 0: # If there is only bridging ligands we use the first set of prefixes
        name += coeff_name1[2]

    #----------- TERMINAL LIGANDS---------------#
    for i, (ligand_name, coeff) in enumerate(ligands_with_coeffs):
        if coeff > 0:
            if should_use_the_coeff_name2(ligand_name) == True:
                prefixe_ligand = coeff_name2[coeff/n] 
                name += (f"{prefixe_ligand}({ligand_name})")
            else:
                prefixe_ligand = coeff_name1[coeff/n] 
                name += prefixe_ligand + ligand_name

    # 6. We add the metal name and already put the capital at the begining (avoid interaction between .capitalize and roman numbers)
    name += metal_name
    name = name.capitalize()

    # 7. We add the charge according to preference selected in the site (roman/int)

    #charge_int = complexe_charge(formula)
    #name += (f"({charge_int})")

    charge = metal_charge(formula)
    charge_roman = roman.toRoman(abs(charge))
    if charge == 0:
        charge_roman = 0
    elif charge < 0:
        charge_roman = "-" + roman.toRoman(abs(charge))
    name += (f"({charge_roman})")

    # 8. Last modifs of the name 
    if last_parenthesis == True:
        name += ")"
    if name.startswith("-"):
        name = name[1:]

    return name

# =========================
# STABILITY MODULE (NEW)
# =========================

def ligand_field_strength(formula):
    """
    Estimation simple du champ de ligand (qualitative -> numérique)
    """
    ligands = ligands_list(formula)

    field_score = 0

    for lig in ligands:
        if lig.startswith("m-"):
            lig = lig[2:]

        info = data_ligands.get(lig, {})

        # classification simple (tu peux enrichir plus tard)
        if info.get("field") == "strong":
            field_score += 2
        elif info.get("field") == "medium":
            field_score += 1
        else:
            field_score += 0

    return field_score


def crystal_field_stabilization(formula):
    """
    Approximation CFSE (très simplifiée)
    """
    d = oxidation_state(formula)
    electrons = electron_count(formula)

    # approximation: d electron count influence
    return (electrons - 6) * ligand_field_strength(formula)


def stability_index(formula):
    """
    Score global de stabilité (0-100)
    """
    cfse = crystal_field_stabilization(formula)
    charge = abs(complexe_charge(formula))
    lig_score = ligand_field_strength(formula)

    score = 50

    # CFSE stabilise
    score += cfse * 5

    # charge élevée = moins stable
    score -= charge * 5

    # ligands forts stabilisent
    score += lig_score * 3

    # clamp 0-100
    return max(0, min(100, score))



# =========================
# ISOMERS
# =========================

stereoisomers_dico = {
    "Ma6": 1,
    "Ma5b1": 1,
    "Ma4b2": 2,
    "Ma3b3": 2,
    "Ma4b1c1": 2,
    "Ma3b1c1d1": 5,
    "Ma2b1c1d1e1": 15,
    "Ma1b1c1d1e1f1": 30,
    "Ma2b2c2": 6,
    "Ma2b2c1d1": 8,
    "Ma3b2c1": 3
}

enantiomers_dico = {
    "Ma6": 0,
    "Ma5b1": 0,
    "Ma4b2": 0,
    "Ma3b3": 0,
    "Ma4b1c1": 0,
    "Ma3b1c1d1": 1,
    "Ma2b1c1d1e1": 6,
    "Ma1b1c1d1e1f1": 15,
    "Ma2b2c2": 1,
    "Ma2b2c1d1": 2,
    "Ma3b2c1": 0
}


def isomers(formula):
    key = ""
    number = []
    alphabet = string.ascii_lowercase
    data = (parse_ligands(formula))
    if len(parse_metal(formula)) == 1:
        key += "M"
    else:
        key += "M2"
    for n in range(len(data[1])):
        number.append(int(data[1][n]))
    number.sort(reverse=True)
    for n in range(len(data[1])):
        letter = alphabet[n]
        key += letter + str(number[n])

    stereo = stereoisomers_dico.get(key)
    enantio = enantiomers_dico.get(key)
    return stereo, enantio


# =============================================================================================================================================================== #

# Final function which prints all the relevant information about the coordination compound
def analyze_complexe(formula):
    parse_elements(formula)
    lines = []

    #Formula
    lines.append(f"* **Formula** : {clean_formula(formula)}")

    # Nomenclature
    name = naming_compound(formula)
    lines.append(f"* **IUPAC Name** : {name}")

    # Metal charge
    metals = parse_metal(formula)
    charge = metal_charge(formula) 
    charge_str = f"{charge}+" if charge > 0 else f"{charge}"
    lines.append(f"* **Metal** : {metals[0]} ({charge_str})")
    
    # Electronic structure
    e_list = electronic_structure(formula)
    lines.append(f"* **Electronic structure** : [{e_list[0]}] {e_list[1]}s{e_list[2]} {e_list[1]-1}d{e_list[3]}")

    # Electrons counting
    count = electron_count(formula)
    lines.append(f"* **Electron count** : {count}")

     # Isomers
    if isomers(formula)[0] == None or isomers(formula)[1] == None:
        lines.append("* **Isomers:** The number of isomers of this compound is not specified")
    else:
        lines.append(f"* **Isomers:** This compound has {isomers(formula)[0]} stereoisomers and {isomers(formula)[1]} enantiomeres pairs")

    # Remarks
    lines.append("* **Remarks:** ... ")

        # Stability (NEW)
    stability = stability_index(formula)
    lines.append(f"* **Stability index** : {stability}/100")

    return lines, "\n".join(lines)


def show_analysis(formula):
    return display(Markdown(analyze_complexe(formula)[1]))

resultat = analyze_complexe("[Fe(CN)6]3-")
for ligne in resultat[0]:
    print(ligne)
#formula = "[Co2(m-OH)(m-NH2)(NH3)8]4+"


#print(clean_formula(formula))
#print(parse_metal(formula))
#print(parse_ligands(formula))
#print(ligands_list(formula))
#print(complexe_charge(formula))
#print(parse_elements(formula))
#print(ligands_charge(formula))
#print(oxidation_state(formula))
#print(metal_charge(formula))
#print(count_bridging_ligands(formula))
#print(bond_order(formula))
#print(naming_compound(formula))

#print(analyze_complexe(formula))



# =============================================================================================================================================================== #
# =============================================================================================================================================================== #
# =============================================================================================================================================================== #
# 3D part of the code
# =============================================================================================================================================================== #
# =============================================================================================================================================================== #
# =============================================================================================================================================================== #


from ase import Atoms
from ase.visualize import view
from ase.data import covalent_radii, atomic_numbers
from pyrolite.geochem.ind import get_ionic_radii


#------------------------------------------------------------
#------------------------------------------------------------
#GEOMETRY ---- METAL - LIGAND -------------------------------
#------------------------------------------------------------
#------------------------------------------------------------

def linear(r):
    return [
        (r, 0, 0),
        (-r, 0, 0)
    ]


def tetrahedral(r):
    base = np.array([
        [ 1,  1,  1],
        [-1, -1,  1],
        [-1,  1, -1],
        [ 1, -1, -1]
    ])
    
    base = base / np.linalg.norm(base[0])  # normalisation
    base = r * base
    
    return [tuple(v) for v in base]


def octahedral(r):
    return [
        ( r, 0, 0),
        (-r, 0, 0),
        (0,  r, 0),
        (0, -r, 0),
        (0, 0,  r),
        (0, 0, -r)
    ]


def trigonal_planar(r):
    return [
        ( r, 0, 0),
        (-r/2,  r*np.sqrt(3)/2, 0),
        (-r/2, -r*np.sqrt(3)/2, 0)
    ]


def trigonal_bipyramidal(r):
    return [
        (0, 0, r),
        (0, 0, -r),
        (r, 0, 0),
        (r * np.cos(np.radians(120)), r * np.sin(np.radians(120)), 0),
        (r * np.cos(np.radians(240)), r * np.sin(np.radians(240)), 0)
    ]


def square_planar(r):
    array = [( r, 0, 0),
            (-r, 0, 0),
            (0,  r, 0),
            (0, -r, 0)]
    return array


def find_geometry(formula, r):
    cn = len(ligands_list(formula))
    if cn == 1: 
        return [(r,0,0)]
    elif cn == 2: 
        return linear(r)
    elif cn == 3: 
        return trigonal_planar(r)
    elif cn == 4 and oxidation_state(formula) == 8: 
        return square_planar(r)
    elif cn == 4 and not oxidation_state(formula) == 8:
        return tetrahedral(r)
    elif cn == 5:
        return trigonal_bipyramidal(r)
    elif cn == 6:
        return octahedral(r)
    else:
        raise ValueError("Error: The visualisation 3D does not work for CN over 6")

#------------------------------------------------------------
#------------------------------------------------------------
# GEOMETRY INTERNE ---- LIGANDS -----------------------------
#------------------------------------------------------------
#------------------------------------------------------------

def ligand_linear(ligand, ligand_coord, r):
    ligand_position = np.array(ligand_coord)
    v = ligand_position / np.linalg.norm(ligand_position) 

    inter_distance = data_ligands[ligand]["inter_distance"]
    position = (v*(inter_distance + r))
    return [tuple(float(x) for x in position)]


def ligand_dlinear(ligand, ligand_coord, r):
    ligand_position = np.array(ligand_coord)
    v = ligand_position / np.linalg.norm(ligand_position) 

    inter_distance = data_ligands[ligand]["inter_distance"]
    inter_distance2 = data_ligands[ligand]["inter_distance2"] 
    position1 = (v*(inter_distance + r))
    position2 = (v*(inter_distance2 + inter_distance + r))
    return [tuple(float(x) for x in position1)] + [tuple(float(x) for x in position2)]


def ligand_bent(ligand, ligand_coord, r):
    ligand_position = np.array(ligand_coord)
    v = ligand_position / np.linalg.norm(ligand_position) 

    inter_distance = data_ligands[ligand]["inter_distance"]
    inter_distance2 = data_ligands[ligand]["inter_distance2"]

    if abs(v[0]) > 0.1:
        temp_vec = np.array([0, 1, 0])
    else:
        temp_vec = np.array([1, 0, 0])

    perp = np.cross(v, temp_vec)
    perp /= np.linalg.norm(perp)

    theta = np.deg2rad(60)
    position1 = (np.cos(theta)*inter_distance + r)*v + (np.sin(theta)*inter_distance)*perp
    position2 = (np.cos(theta)*inter_distance2 + r)*v + (-np.sin(theta)*inter_distance2)*perp
    print(position1, position2)
    return [tuple(float(x) for x in position1)] + [tuple(float(x) for x in position2)]


def ligand_tetrahedral(ligand, ligand_coord, r):
    ligand_position = np.array(ligand_coord)
    v = ligand_position / np.linalg.norm(ligand_position) 

    inter_distance = data_ligands[ligand]["inter_distance"]

    if abs(v[0]) > 0.1:
        temp_vec = np.array([0, 1, 0])
    else:
        temp_vec = np.array([1, 0, 0])

    u = np.cross(v, temp_vec)
    u /= np.linalg.norm(u)
    w = np.cross(v, u)

    positions = []
    theta = np.deg2rad(-54.75)
    for i in range(3):
        phi = np.deg2rad(i * 120)
        pos = ((np.cos(theta)*inter_distance + r)*v + np.sin(theta) * np.cos(phi) * inter_distance* u + np.cos(theta) * np.sin(phi) * inter_distance * w)
        positions.append(tuple(float(x) for x in pos))  
    
    return positions


def get_geometry_ligand(ligand_input):
    geometry = data_ligands[ligand_input].get("geometry")
    if geometry != None:
        return geometry
    return False

#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
#-------- 3D visualisation and coumpound creation -----------
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------


def atoms_position_and_bond(formula):
    r = 1.7
    bonding = []
    nb_of_atoms = 0
    position = [(0,0,0)]
    big_array = find_geometry(formula, r)
    ligand_list = ligands_list(formula)
    for i, ligand in enumerate(ligand_list):
        if get_geometry_ligand(ligand) == "sphere":
            nb_of_atoms += 1
            position += [big_array[i]]
            bonding += (0,nb_of_atoms)
        elif get_geometry_ligand(ligand) == "linear":
           nb_of_atoms += 2
           position += [big_array[i]]
           position += ligand_linear(ligand, big_array[i], r)
           bonding += (0,nb_of_atoms-1)
        elif get_geometry_ligand(ligand) == "dlinear":
            nb_of_atoms += 3
            position += [big_array[i]]
            position += ligand_dlinear(ligand, big_array[i], r)
            bonding += (0,nb_of_atoms-2)
        elif get_geometry_ligand(ligand) == "bent":
            nb_of_atoms += 3
            position += [big_array[i]]
            position += ligand_bent(ligand, big_array[i], r)
            bonding += (0,nb_of_atoms-2)
        elif get_geometry_ligand(ligand) == "tetrahedral":
            nb_of_atoms += 4
            position += [big_array[i]]
            position += ligand_tetrahedral(ligand, big_array[i], r)
            bonding += (0,nb_of_atoms-3)
        else:
            raise ValueError("Error: Geometry of the ligand not available in 3D")
    return position, bonding


def get_atoms(ligand_input):
    ligand_info = data_ligands.get(ligand_input)
    donor_atom = ligand_info.get("donor_atoms")
    result = [donor_atom[0]] if donor_atom[0] else []
    atoms = re.findall(r'[A-Z][a-z]?\d*', ligand_input)

    for atom in atoms:
        match = re.match(r'([A-Z][a-z]?)(\d*)', atom)
        symbol = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1

        if symbol == donor_atom[0]:
            continue
        else:
            result.extend([symbol] * count)
    return result


def atom_symbols(formula):
    metal = parse_metal(formula)
    atoms_list = metal

    ligand_list = ligands_list(formula)
    for ligand in ligand_list:
        atoms_list += get_atoms(ligand)
    return atoms_list


def create_compound_render(formula):
    compound = Atoms(atom_symbols(formula), positions= atoms_position_and_bond(formula)[0])
    return compound 



"""
def metal_radii(formula):
    metal = parse_metal(formula)[0]
    coordination = len(ligands_list(formula))
    charge = metal_charge(formula)
    radii = get_ionic_radii(metal, charge, coordination)
    return radii

def ligand_radii(formula):
    ligands = ligands_list(formula)
    r = 0
    for ligand in ligands:
        if data_ligands[ligand].get("radii") != False:
            r += data_ligands[ligand]["radii"]
        else:
            donor_atom = data_ligands[ligand]["donor_atoms"][0]
            r += covalent_radii[atomic_numbers[donor_atom]]
    return r/len(ligands)
"""
