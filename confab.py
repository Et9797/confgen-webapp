import os
from openbabel import openbabel as ob
from openbabel import pybel


def generate_conformers(*args, force_field):

    """ Generate conformers"""

    molecule = args[0]
    if molecule.endswith(".sdf"):
        mole = next(pybel.readfile("sdf", molecule))
    else:
        mole = pybel.readstring("can", molecule)
    mole.addh()
    mole.make3D()
    mole = mole.OBMol
    mole.SetChainsPerceived(True) #temp fix for github problem # 
    
    ff = ob.OBForceField_FindType(force_field)
    ff.Setup(mole)
    ff.DiverseConfGen(0.5, 100000, 50.0, True)
    ff.GetConformers(mole)
    return mole


def write_conformers(mole, confid_list, output_ext, seperate_files):

    """Write conformers to file"""

    cv = ob.OBConversion()
    cv.SetOutFormat(output_ext.lower())
    if seperate_files == "True":
        for index, conf_nr in enumerate(confid_list):
            mole.SetConformer(conf_nr)
            cv.WriteFile(mole, f"conformer_{index}.{output_ext.lower()}")
    else:
        for index, conf_nr in enumerate(confid_list):
            mole.SetConformer(conf_nr)
            cv.WriteFile(mole, f"conformer_{index}.{output_ext.lower()}")
        with open(f"ConformersMerged.{output_ext.lower()}", "a") as merged:
            for f in os.listdir():
                if f.startswith("conformer_"):
                    with open(f, "r") as conf:
                        conformer = conf.read()
                        merged.write(conformer)
        [os.remove(f) for f in os.listdir() if f.startswith("conformer_")]


