from rdkit import Chem
import re
import json






# Loading the data from the metals.json file
with open("data/metals.json", "r") as jsonfile:
    data_metals = json.load(jsonfile)

# Loading the data from the metals.json file
with open("data/ligands.json", "r") as jsonfile:
    data_ligands = json.load(jsonfile)





def parse_formula(formula):

    clean_formula = formula.replace("[", "").replace("]", "")
    

    metal_match = re.search(r'^([A-Z][a-z]?)', clean_formula)
    for key in data_metals:
        if key == metal_match
        
        
        print(key)



    metal = metal_match.group(0)
    

    ligand_pattern = r'\(([A-Za-z0-9]+)\)(\d*)|([A-Z][a-z]?)(\d*)'
    matches = re.findall(ligand_pattern, clean_formula[len(metal):])
    
    ligand_list = []
    for m in matches:

        name = m[0] if m[0] else m[2]
        count = m[1] if m[1] else (m[3] if m[3] else 1)
        
        for _ in range(int(count)):
            ligand_list.append(name)
            
    return metal, ligand_list



metal, ligands = parse_formula("[CO(CO)5]")
print(f"Metal: {metal}")    # Output: Fe
print(f"Ligands: {ligands}") # Output: ['CO', 'CO', 'CO', 'CO', 'CO']