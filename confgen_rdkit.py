from rdkit import Chem
from rdkit.Chem import AllChem
#import Mol2Writer

def gen_conformers(*args, no_conformers):  

    """Generates conformers for (crystal) ligand"""

    molecule = args[0]
    if molecule.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(molecule)
        mol = next(suppl) 
    else:
        mol = Chem.MolFromSmiles(molecule)
    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMultipleConfs(mol, no_conformers, clearConfs=True) #numThreads=0
    return mol


def write_confs_to_file(molecule, output_ext, seperate_files):

    """Writes conformers to files"""

    conf_ids = list(range(molecule.GetNumConformers()))
    
    if output_ext == "PDB":
        if seperate_files == "True":
            for cid in conf_ids:
                Chem.MolToPDBFile(molecule, f"conformer_{cid}.pdb", confId=cid)
        else:
            pdb_writer = Chem.PDBWriter("ConformersMerged.pdb")
            for cid in conf_ids:
                pdb_writer.write(molecule, confId=cid)
            pdb_writer.close()
    if output_ext == "SDF":
        if seperate_files == "True":
            for cid in conf_ids:
                sd_writer = Chem.SDWriter(f"conformer_{cid}.sdf")
                sd_writer.write(molecule, confId=cid)
                sd_writer.close()
        else:
            sd_writer = Chem.SDWriter("ConformersMerged.sdf")
            for cid in conf_ids:
                sd_writer.write(molecule, confId=cid)
            sd_writer.close()
