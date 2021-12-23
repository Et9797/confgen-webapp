from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(molecule, no_conformers):

    """Generates conformers for (crystal) ligand"""

    if molecule.endswith("sdf"):
        sdf_supplier = Chem.SDMolSupplier(molecule)
        mol = next(sdf_supplier)
    elif molecule.endswith("mol"):
        mol = Chem.MolFromMolFile(molecule)
    else:
        mol = Chem.MolFromSmiles(molecule)

    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMultipleConfs(mol, no_conformers, clearConfs=True) #numThreads=0

    return mol

def write_confs_to_file(molecule, output_ext, seperate_files):

    """Writes conformers to PDB/SDF/Mol files"""

    conf_ids = list(range(molecule.GetNumConformers()))

    if seperate_files == "on":
        if output_ext == "pdb":
            for cid in conf_ids:
                Chem.MolToPDBFile(molecule, f"conformer_{cid}.pdb", confId=cid)
        if output_ext == "sdf":
            for cid in conf_ids:
                sd_writer = Chem.SDWriter(f"conformer_{cid}.sdf")
                sd_writer.write(molecule, confId=cid)
                sd_writer.close()
        if output_ext == "mol":
            for cid in conf_ids:
                Chem.MolToMolFile(molecule, f"conformer_{cid}.mol", confId=cid)
    if seperate_files == "off":
        if output_ext == "pdb":
            pdb_writer = Chem.PDBWriter("ConformersMerged.pdb")
            for cid in conf_ids:
                pdb_writer.write(molecule, confId=cid)
            pdb_writer.close()
        if output_ext == "sdf":
            sd_writer = Chem.SDWriter("ConformersMerged.sdf")
            for cid in conf_ids:   
                sd_writer.write(molecule, confId=cid)
            sd_writer.close()
                                                                                            