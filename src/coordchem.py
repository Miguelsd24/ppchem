import sys
import re
import json
from pathlib import Path

# Allows to display only one line when running into an error in the jupyter notebook
def hide_traceback(exc_tuple=None, filename=None, tb_offset=None, exception_only=False, running_compiled_code=False):
    etype, value, tb = sys.exc_info()
    return print(f"{value}")
get_ipython().showtraceback = hide_traceback # This replaces the default Jupyter error handler with our clean one-liner

# Loading the data from the metals/ligands.json file

BASE_DIR = Path(__file__).resolve().parent.parent

with open(BASE_DIR / "data" / "metals.json") as f:
    data_metals = json.load(f)

with open(BASE_DIR / "data" / "ligands.json") as f:
    data_ligands = json.load(f)

# The first part of this file focuses on the formula format processing, along with the calculations (e.g. electron counting, oxidation state of the metal, ...)

# Function which extract the metal and its stoechiometric coefficient from a given formula. It also verifies if the metal is in the json database
def parse_metal(formula):

    # 1. The metal string is isolated possibly along with the stoechiometric coefficient. Also, we verify that the input str has the appropriate format
    start = formula.find("[")
    end = formula.find("(")

    if start == -1 or end == -1  or start >= end:
        raise ValueError("Error: Wrong input format.")

    clean_formula = formula[start + 1:end]

    # 2. We use re.match() to find the metal/coefficient pattern

    match = re.match(r"^(?:(b|d|t))?([A-Z][A-Za-z]*)([1-9]\d*)?$", clean_formula)

    if match == None:
        raise ValueError("Error: Wrong input format.")

    # 2. We test if the metal is present in the database and treat the output
    metal = match.group(2)
    if metal not in data_metals:
        raise ValueError(f"Error: Metal {metal} is not in the database.")

    # 3. We test if the coefficient has the appropriate value or format and treat the output
    num = match.group(3)
    coeff = test_string_for_int(num)
    if coeff > 2:
        raise ValueError("Error: Invalid coefficient, the coordination compound has more than two metal centers.")
    
    # 4. We put the metal and its coefficient in a list and return it
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
    for ligand_key in data_ligands:
        if "m-" + ligand_key == ligand_input:
            return "m-" + ligand_key
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
    for ligand, n in match:
        if ligand == "":
            raise ValueError(f"Error: At least one ligand is missing in the formula")
        if find_ligand(ligand) == False:
            raise ValueError(f"Error: Ligand {ligand} not in the database")
        else:
            coeff = test_string_for_int(n)
            result.extend([find_ligand(ligand)])

    return result, coeff 
    
# Function which creates a list with all the ligands times their coeficient. This list is useful to use loops in the calculation function like that There is non need to deal with number later.
def ligands_list(formula):
    ligands_list = []
    ligands = parse_ligands(formula)[0]
	coeff = parse_ligands(formula)[1]
    for n in len(ligands):
        ligands_list.append(ligands[n]*coeff[n])        
    return ligands_list

# Function which counts the number of bridging ligands
def count_bridging_ligands(formula):
    num = 0
    for ligand in ligands_list(formula)]:
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
    return ox_state

# Function which does the electron counting
def electron_count(formula):
    electrons = oxidation_state(formula)*len(parse_metal(formula))
    ligands = ligands_list(formula)

    for ligand in ligands:
        if ligand.startswith("m-"):
            electrons += data_ligands[ligand[2:]]["bridging_e"]
        else:
            electrons += data_ligands[ligand]["donor_e"]
    if bond_order(formula) > 0:
        electrons += 2*bond_order(formula)

    return electrons // len(parse_metal(formula))

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

# Final function which prints all the relevant information about the coordination compound
def analyze_complexe(formula):

    # Metal charge part
    metals = parse_metal(formula)
    charge = metal_charge(formula) 
    if charge > 0:
        print(f"Metal : {metals[0]} {charge}+")
    elif charge == 0:  
        print(f"Metal : {metals[0]} {charge}")
    
    # Electronic structure part
    list = electronic_structure(formula)
    print(f"Metal electronic structure : [{list[0]}] {list[1]}s{list[2]} {list[1]-1}d{list[3]}")

    # Electrons counting part
    print(f"Electron count : {electron_count(formula)}")


def name_or_abbr_ligands(ligand):
    abbr_ligand = data_ligands.get(ligand["abbr"])
    if ligand["abbr"] is None:
        return data_ligands[ligand]["name"]
    else:
        return abbr_ligand
        
coeff_name1 = {
2: "di"
3: "tri"
4: "tetra"
5: "penta"
6: "hexa"
7: "hepta"
8: "octa"
9: "nona"
10: "deca"
11: ""
12: ""

        
coeff_name2 = {
2: "bis"
3: "tris"
4: "tetrakis"
5: "pentakis"
6: "hexakis"
7: "heptakis"
8: "octakis"
9: "


    
    
def naming_compound(formula):
    ligands = parse_ligands(formula)[0]
    ligands_coeff = parse_ligands(formula)[1]
    metals = parse_metal(formula)
    name = ""
    for ligand in ligands:
        ligands_name_or_abbr = [name_or_abbr_ligands(ligand) for ligand in ligands]
        ligands_name_or_abbr.sort(key=str.lower)
    for n in range(len(ligands_name_or_abbr)):
        name += 


    return name
