from rdkit import Chem
from rdkit.Chem import AllChem
import os

def generate_conformers(path_to_molecule, no_conformers):

    """Generates conformers for (crystal) ligand"""

    if path_to_molecule.endswith("sdf"):
        sdf_supplier = Chem.SDMolSupplier(path_to_molecule)
        mol = next(sdf_supplier)
    elif path_to_molecule.endswith("mol"):
        mol = Chem.MolFromMolFile(path_to_molecule)
    else:
        mol = Chem.MolFromSmiles(path_to_molecule)

    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMultipleConfs(mol, no_conformers, clearConfs=True) #numThreads=0

    return mol

def write_confs_to_file(molecule, mol_path, output_ext, seperate_files):

    """Writes conformers to PDB/SDF/Mol files"""

    conf_ids = list(range(molecule.GetNumConformers()))

    if seperate_files == "on":
        if output_ext == "pdb":
            for cid in conf_ids:
                Chem.MolToPDBFile(molecule, os.path.join(mol_path, f"conformer_{cid}.pdb"), confId=cid)
        if output_ext == "sdf":
            for cid in conf_ids:
                sd_writer = Chem.SDWriter(os.path.join(mol_path, f"conformer_{cid}.sdf"))
                sd_writer.write(molecule, confId=cid)
                sd_writer.close()
        if output_ext == "mol":
            for cid in conf_ids:
                Chem.MolToMolFile(molecule, os.path.join(mol_path, f"conformer_{cid}.mol"), confId=cid)
    if seperate_files == "off":
        if output_ext == "pdb":
            pdb_writer = Chem.PDBWriter(os.path.join(mol_path, "ConformersMerged.pdb"))
            for cid in conf_ids:
                pdb_writer.write(molecule, confId=cid)
            pdb_writer.close()
        if output_ext == "sdf":
            sd_writer = Chem.SDWriter(os.path.join(mol_path, "ConformersMerged.sdf"))
            for cid in conf_ids:   
                sd_writer.write(molecule, confId=cid)
            sd_writer.close()
                                                                                            
