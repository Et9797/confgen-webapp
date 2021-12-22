from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(molecule, no_conformers):  

    """Generates conformers for (crystal) ligand"""
    
    match molecule.split(".")[-1]:
        case "sdf":
            sdf_supplier = Chem.SDMolSupplier(molecule)
            mol = next(sdf_supplier) 
        case "mol":
            mol = Chem.MolFromMolFile(molecule)
        case _:
            mol = Chem.MolFromSmiles(molecule)

    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMultipleConfs(mol, no_conformers, clearConfs=True) #numThreads=0

    return mol


def write_confs_to_file(molecule, output_ext, seperate_files):

    """Writes conformers to PDB/SDF/Mol files"""

    conf_ids = list(range(molecule.GetNumConformers()))

    match seperate_files:
        case "on":
            match output_ext:
                case "pdb":
                    for cid in conf_ids:
                        Chem.MolToPDBFile(molecule, f"conformer_{cid}.pdb", confId=cid)
                case "sdf":
                    for cid in conf_ids:
                        sd_writer = Chem.SDWriter(f"conformer_{cid}.sdf")
                        sd_writer.write(molecule, confId=cid)
                        sd_writer.close()
                case "mol":
                    for cid in conf_ids:
                        Chem.MolToMolFile(molecule, f"conformer_{cid}.mol", confId=cid)
        case "off":
            match output_ext:
                case "pdb":
                    pdb_writer = Chem.PDBWriter("ConformersMerged.pdb")
                    for cid in conf_ids:
                        pdb_writer.write(molecule, confId=cid)
                    pdb_writer.close()
                case "sdf":
                    sd_writer = Chem.SDWriter("ConformersMerged.sdf")
                    for cid in conf_ids:
                        sd_writer.write(molecule, confId=cid)
                    sd_writer.close()
