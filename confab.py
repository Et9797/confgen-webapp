from openbabel import openbabel as ob
from openbabel import pybel


def generate_conformers(smiles, force_field):
    mole = pybel.readstring("can", smiles)
    mole.addh()
    mole.make3D()
    mole = mole.OBMol
    mole.SetChainsPerceived(True) #temp fix for github problem # 
    
    ff = ob.OBForceField_FindType(force_field) #amber
    ff.Setup(mole)
    ff.DiverseConfGen(0.5, 100000, 50.0, True)
    ff.GetConformers(mole)
    return mole


def write_conformers(mole, n_conformers):
    #n_conformers is a list with conf_ids 
    cv = ob.OBConversion()
    cv.SetOutFormat("pdb")
    for index, conf_nr in enumerate(n_conformers):
        print(conf_nr)
        mole.SetConformer(conf_nr)
        cv.WriteFile(mole, f"conformer_{index}.pdb")