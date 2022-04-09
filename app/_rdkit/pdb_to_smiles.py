import re
import requests


def convert(ligand_pdb):
    
    """Convert PDB to Smiles in order to assign bond orders"""

    url = 'https://cactus.nci.nih.gov/cgi-bin/translate.tcl'

    with open(ligand_pdb, "rb") as file:
        pdb_file = {'file': file}
        data = {
            "smiles": "C12C3C4C1C5C4C3C25",
            "format": "screen",
            "astyle": "kekule",
            "dim": "2D"
        }
        r = requests.post(url, files=pdb_file, data=data)

    pattern = re.compile(r"<B>.+</B>")
    smiles = re.search(pattern, r.text).group(0)[3:-4]

    return smiles