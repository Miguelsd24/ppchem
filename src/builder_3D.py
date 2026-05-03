import os
import logging
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# =========================
# SAFE INPUT NORMALIZATION
# =========================
def normalize_ligands(ligands):
    """
    Convertit tous les formats possibles en liste plate :
    [("N", 6)] → ["N","N","N","N","N","N"]
    """
    if ligands is None:
        raise ValueError("Ligands is None")

    flat = []

    for item in ligands:

        # case tuple ("N", 6)
        if isinstance(item, tuple) or isinstance(item, list):
            if len(item) != 2:
                raise ValueError(f"Invalid ligand tuple: {item}")

            smi, n = item

            if smi is None:
                raise ValueError("Ligand SMILES is None")

            n = int(n)
            flat.extend([smi] * n)

        # case already string
        elif isinstance(item, str):
            flat.append(item)

        else:
            raise ValueError(f"Unknown ligand format: {item}")

    return flat


# =========================
# LIGAND GENERATION
# =========================
def create_ligand(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # SAFE embedding (avoid crash)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except:
        # fallback minimal embedding
        AllChem.EmbedMolecule(mol)

    return mol


# =========================
# DONOR DETECTION
# =========================
def get_donor_atoms(mol):
    donors = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O", "S", "P"]:
            donors.append(atom.GetIdx())

    return donors if donors else [0]


# =========================
# GEOMETRIES
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
# GEOMETRY SELECTOR
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
# DISTANCES
# =========================
BOND_LENGTH = {
    "Fe": 2.2,
    "Co": 2.1,
    "Ni": 2.0,
    "Cu": 2.0,
    "Zn": 2.1
}


# =========================
# CORE BUILDER
# =========================
def _build_complex_3d(metal="Fe", ligands=[("N", 6)]):

    logger.info(f"Building {metal}")

    ligands = normalize_ligands(ligands)  # 🔥 FIX CRUCIAL

    mol = Chem.RWMol()
    conf = Chem.Conformer()

    metal_idx = mol.AddAtom(Chem.Atom(metal))
    conf.SetAtomPosition(metal_idx, (0.0, 0.0, 0.0))

    n = len(ligands)

    geom = choose_geometry(metal, n)
    r = BOND_LENGTH.get(metal, 2.2)

    positions = get_positions(geom, r)

    for i, smi in enumerate(ligands):

        ligand = create_ligand(smi)
        donor_atoms = get_donor_atoms(ligand)
        donor = donor_atoms[0]

        lconf = ligand.GetConformer()
        offset = positions[i % len(positions)]

        atom_map = {}

        for atom in ligand.GetAtoms():
            new_idx = mol.AddAtom(Chem.Atom(atom.GetSymbol()))
            atom_map[atom.GetIdx()] = new_idx

            p = lconf.GetAtomPosition(atom.GetIdx())
            v = np.array([p.x, p.y, p.z])

            conf.SetAtomPosition(new_idx, v + offset)

        for bond in ligand.GetBonds():
            a1 = atom_map[bond.GetBeginAtomIdx()]
            a2 = atom_map[bond.GetEndAtomIdx()]
            mol.AddBond(a1, a2, bond.GetBondType())

        mol.AddBond(metal_idx, atom_map[donor], Chem.BondType.SINGLE)

    mol.AddConformer(conf)

    return mol.GetMol()


# =========================
# SAFE WRAPPER
# =========================
def build_3d_safe(metal, ligands):
    try:
        mol = _build_complex_3d(metal, ligands)

        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError("Invalid molecule generated")

        return mol

    except Exception as e:
        print("[3D ERROR]", e)
        return None