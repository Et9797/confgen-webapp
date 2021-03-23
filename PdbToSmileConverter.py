import re
import requests

#temp fix for converting pdb to smiles, need local solution

def pdb_to_smiles(ligand_pdb):
    
    """Convert PDB to smiles using NIH's PDB -> smiles converter"""

    url = 'https://cactus.nci.nih.gov/cgi-bin/translate.tcl'

    with open(ligand_pdb, "rb") as file:
        try:
            pdb_file = {'file': file}
            data = {
                "smiles": "C12C3C4C1C5C4C3C25",
                "format": "screen",
                "astyle": "kekule",
                "dim": "2D"
            }
            r = requests.post(url, files=pdb_file, data=data)
        except requests.exceptions.HTTPError as e:
            raise e 

    pattern = re.compile(r"<B>.+</B>")
    smiles = re.search(pattern, r.text).group(0)[3:-4]
    return smiles
