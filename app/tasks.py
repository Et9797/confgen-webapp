from pathlib import Path
from app import make_celery, app
from ._rdkit import pdb_to_smiles, confgen

celery = make_celery(app)

@celery.task()
def generate_confs(smiles, mol_filename, mol_path, no_conformers, output_ext, output_separate):
    # Generate conformers
    if smiles:
        conformers = confgen.generate_conformers(smiles, no_conformers=no_conformers)
    elif mol_filename.split(".")[-1] == 'pdb':
        smiles = pdb_to_smiles.convert(Path(mol_path, mol_filename))
        conformers = confgen.generate_conformers(smiles, no_conformers=no_conformers)
    else:
        conformers = confgen.generate_conformers(Path(mol_path, mol_filename),
                                                 no_conformers=no_conformers
                                                 )
    # Write conformers to disk
    confgen.write_confs_to_file(conformers, mol_path, output_ext, output_separate)