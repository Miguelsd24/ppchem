import os
import logging
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# =========================
# LIGAND GENERATION
# =========================
def create_ligand(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    return mol


# =========================
# DONOR DETECTION (simple but robust)
# =========================
def get_donor_atoms(mol):
    donors = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O", "S", "P"]:
            donors.append(atom.GetIdx())

    return donors if donors else [0]


# =========================
# GEOMETRIES (clean chemical sets)
# =========================
def linear(r):
    return [np.array([r, 0, 0]), np.array([-r, 0, 0])]


def tetrahedral(r):
    return [
        np.array([ r,  r,  r]),
        np.array([-r, -r,  r]),
        np.array([-r,  r, -r]),
        np.array([ r, -r, -r]),
    ]


def octahedral(r):
    return [
        np.array([ r, 0, 0]),
        np.array([-r, 0, 0]),
        np.array([0,  r, 0]),
        np.array([0, -r, 0]),
        np.array([0, 0,  r]),
        np.array([0, 0, -r]),
    ]


def square_planar(r):
    return [
        np.array([ r, 0, 0]),
        np.array([-r, 0, 0]),
        np.array([0,  r, 0]),
        np.array([0, -r, 0]),
    ]


# =========================
# GEOMETRY SELECTION
# =========================
def choose_geometry(metal, n):
    if n == 2:
        return "linear"
    if n == 4:
        return "square" if metal == "Cu" else "tetra"
    if n == 6:
        return "octa"
    return "octa"


def get_positions(geom, r):
    if geom == "linear":
        return linear(r)
    if geom == "tetra":
        return tetrahedral(r)
    if geom == "square":
        return square_planar(r)
    return octahedral(r)


# =========================
# METAL–LIGAND DISTANCES (Å)
# =========================
BOND_LENGTH = {
    "Fe": 2.2,
    "Co": 2.1,
    "Ni": 2.0,
    "Cu": 2.0,
    "Zn": 2.1,
    "Cs": 3.5
}


# =========================
# MAIN BUILDER
# =========================
def build_complex(metal="Fe", ligands=[("N", 6)]):

    logger.info(f"Building complex {metal}")

    mol = Chem.RWMol()
    conf = Chem.Conformer()

    # metal center
    metal_idx = mol.AddAtom(Chem.Atom(metal))
    conf.SetAtomPosition(metal_idx, (0.0, 0.0, 0.0))

    # flatten ligand list
    ligand_list = []
    for smi, n in ligands:
        ligand_list.extend([smi] * n)

    n = len(ligand_list)

    geom = choose_geometry(metal, n)
    r = BOND_LENGTH.get(metal, 2.2)

    positions = get_positions(geom, r)

    # =========================
    # PLACE LIGANDS
    # =========================
    for i, smiles in enumerate(ligand_list):

        ligand = create_ligand(smiles)
        donor_atoms = get_donor_atoms(ligand)
        donor = donor_atoms[0]

        lconf = ligand.GetConformer()

        offset = positions[i % len(positions)]

        atom_map = {}

        # atoms
        for atom in ligand.GetAtoms():

            new_idx = mol.AddAtom(Chem.Atom(atom.GetSymbol()))
            atom_map[atom.GetIdx()] = new_idx

            p = lconf.GetAtomPosition(atom.GetIdx())
            v = np.array([p.x, p.y, p.z])

            conf.SetAtomPosition(new_idx, v + offset)

        # bonds inside ligand
        for bond in ligand.GetBonds():

            a1 = atom_map[bond.GetBeginAtomIdx()]
            a2 = atom_map[bond.GetEndAtomIdx()]

            mol.AddBond(a1, a2, bond.GetBondType())

        # metal–donor bond
        mol.AddBond(metal_idx, atom_map[donor], Chem.BondType.SINGLE)

    mol.AddConformer(conf)

    return mol.GetMol()


# =========================
# SAVE
# =========================
def save_sdf(mol, path):
    writer = Chem.SDWriter(path)
    writer.write(mol)
    writer.close()
    logger.info(f"Saved: {path}")


# =========================
# MAIN
# =========================
if __name__ == "__main__":

    mol = build_complex(
        "Fe",
        [("N", 6)]   # ex: Fe(NH3)6
    )

    out = os.path.join(BASE_DIR, "complex.sdf")
    save_sdf(mol, out)