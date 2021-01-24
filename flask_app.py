import io
import os
import sys
import shutil
import re
import zipfile
import requests
import logging
import traceback
import random
from openbabel import openbabel as ob
from openbabel import pybel
import confab 
from rdkit import Chem
from rdkit.Chem import AllChem
import conf_gen_rdkit
import pdbToSmileConverter
from flask import Flask, Response, render_template, request, redirect, url_for, send_file

BASE_DIR = '/home/et/personal_projects/rdkit-obabel-confgen/'
MOLECULE_UPLOADS = '/home/et/personal_projects/rdkit-obabel-confgen/MOLECULE_UPLOADS/'
#change to '/var/www/html/obabel_confgen/MOLECULE_UPLOADS/'
app = Flask(__name__)
app.config["BASE_DIR"] = BASE_DIR
app.config["MOLECULE_UPLOADS"] = MOLECULE_UPLOADS
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0 #don't cache the files


@app.errorhandler(Exception)
def internal_error(exception):
    with open(os.path.join(app.config["BASE_DIR"], "log.txt"), "a") as f:
        f.write(str(exception) + "\n")
        f.write(traceback.format_exc())


@app.route("/reset/<method>/<mol>")
def reset(method, mol):
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], mol)):
        shutil.rmtree(os.path.join(app.config["MOLECULE_UPLOADS"], mol))
        if method == "confab":
            return redirect(url_for("confab_page"))
        else:
            return redirect(url_for("rdkit"))
    else:
        if method == "confab":
            return redirect(url_for("confab_page"))
        else:
            return redirect(url_for("rdkit"))


@app.route("/")
def index():
    return redirect(url_for("confab_page"))


@app.route("/confab")
def confab_page():
    return render_template("index.html")


@app.route("/rdkit")
def rdkit():
    return render_template("rdkit.html")


@app.route("/<method>/<mol>") 
def serve_pdbs(method, mol): 
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], mol)):
        zipfolder = zipfile.ZipFile("Conformers.zip", "w", zipfile.ZIP_STORED)
        for f in os.listdir():
            if f != "Conformers.zip":
                zipfolder.write(f)
        zipfolder.close()
        zip_mem = io.BytesIO()
        with open(zipfolder.filename, "rb") as fo:
            zip_mem.write(fo.read())
            zip_mem.seek(0)
        fo.close()

        shutil.rmtree(os.path.join(app.config["MOLECULE_UPLOADS"], mol))

        return send_file(zip_mem, mimetype="application/zip", as_attachment=True, 
        attachment_filename="Conformers.zip", cache_timeout=0)
    else:
        if method == "confab":
            return redirect(url_for("confab_page"))
        else:
            return redirect(url_for("rdkit"))


@app.route("/<method>", methods=["POST", "GET"])
def form_handler(method):
    if request.method == "POST":
        if method == "confab":
            force_field = request.form["force_field"]
        smiles = request.form["smiles_molecule"]
        pdb_file = request.files["pdb_molecule"]
        no_conformers = int(request.form["no_conformers"])
        pattern = re.compile('[^A-Za-z0-9]+')
        if smiles:
            smiles_no_special_chars = re.sub(pattern, "", smiles)
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], smiles_no_special_chars[0:10])
        else: #PDB was provided
            assert pdb_file.filename.split(".")[-1] == "pdb"
            pdb_temp_path = os.path.join(app.config["MOLECULE_UPLOADS"], pdb_file.filename)
            pdb_file.save(pdb_temp_path)
            smiles = pdbToSmileConverter.pdb_to_smiles(pdb_temp_path)
            smiles_no_special_chars = re.sub(pattern, "", smiles)
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], smiles_no_special_chars[0:10])
            os.remove(pdb_temp_path)
            
        if os.path.exists(mol_path):
            shutil.rmtree(mol_path)
        os.mkdir(mol_path)
        os.chdir(mol_path)
        print(method)

        if method == "confab":
            mole = confab.generate_conformers(smiles, force_field)
            if mole.NumConformers() > no_conformers:
                conf_sample = random.sample(range(mole.NumConformers()), no_conformers)
                confab.write_conformers(mole, conf_sample)
            else:
                confab.write_conformers(mole, range(mole.NumConformers()))
            return render_template("index.html", method="confab", mol=smiles_no_special_chars[0:10])
        else:
            conformers = conf_gen_rdkit.gen_conformers(smiles, no_conformers)
            conf_gen_rdkit.write_confs_to_pdb(conformers)
            return render_template("rdkit.html", method="rdkit", mol=smiles_no_special_chars[0:10])
        

if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)

