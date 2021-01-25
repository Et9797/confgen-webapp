from rdkit import Chem
from rdkit.Chem import AllChem


def gen_conformers(smiles, no_conformers):

    """Generates conformers for (crystal) ligand"""
    
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol, addCoords=True)
    conf_ids = AllChem.EmbedMultipleConfs(mol, no_conformers, numThreads=0, clearConfs=True) #numThreads=0
    num_of_confs=mol.GetNumConformers()
    return mol, conf_ids


def write_confs_to_pdb(confs_tuple):

    """Writes conformers to PDB files"""

    #conformers[0] = mol object with all conformers, conformers[1] are the conf_ids
    conf_ids=list(confs_tuple[1])
    for cid in conf_ids:
        Chem.MolToPDBFile(confs_tuple[0], "conformer_{}.pdb".format(str(cid)), confId=cid)

