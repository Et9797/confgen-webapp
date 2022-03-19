import re
import requests

# For converting pdb to smiles to assign bond orders

def convert(ligand_pdb):
    
    """Convert PDB to smiles using NIH's PDB -> Smiles converter"""

    url = 'https://cactus.nci.nih.gov/cgi-bin/translate.tcl'

    try:
        with open(ligand_pdb, "rb") as file:
            pdb_file = {'file': file}
            data = {
                "smiles": "C12C3C4C1C5C4C3C25",
                "format": "screen",
                "astyle": "kekule",
                "dim": "2D"
            }
            r = requests.post(url, files=pdb_file, data=data)
    except Exception as e:
        raise e 

    pattern = re.compile(r"<B>.+</B>")
    smiles = re.search(pattern, r.text).group(0)[3:-4]

    return smiles
