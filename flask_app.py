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
import time
import uuid
from openbabel import openbabel as ob
from openbabel import pybel
import confab 
from rdkit import Chem
from rdkit.Chem import AllChem
import confgen_rdkit
import pdbToSmileConverter
from flask import Flask, Response, render_template, request, redirect, url_for, send_file
from flask_mail import Mail, Message
from config import mail_username, mail_password

#change to '/var/www/html/private_confgen/MOLECULE_UPLOADS/'
BASE_DIR = '/home/et/personal_projects/private-confgen/'
MOLECULE_UPLOADS = '/home/et/personal_projects/private-confgen/MOLECULE_UPLOADS/'

app = Flask(__name__)
app.config["BASE_DIR"] = BASE_DIR
app.config["MOLECULE_UPLOADS"] = MOLECULE_UPLOADS
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0 #no caching of the files
app.config["MAIL_SERVER"] = "smtp-mail.outlook.com"
app.config["MAIL_PORT"] = 587
app.config["MAIL_USE_TLS"] = True
app.config["MAIL_USE_SSL"] = False
app.config["MAIL_USERNAME"] = mail_username
app.config["MAIL_PASSWORD"] = mail_password

mail = Mail(app)

@app.errorhandler(Exception)
def internal_error(exception):
    with open(os.path.join(app.config["BASE_DIR"], "error.log"), "a") as f:
        f.write(time.strftime("%d/%m %H:%M:%S") + "\n")
        exc_type, _ , _ = sys.exc_info()
        if exc_type.__name__ == "ArgumentError":
            exc_type.__name__= "SMILES Parse Error"
        f.write(f"{exc_type.__name__}\n\n")
    current_page = request.path.strip("/")
    return render_template(f"{current_page}.html", err=str(exc_type.__name__))


@app.route("/contact", methods=["POST","GET"])
def contact():
    if request.method == "POST":
        email = request.form["email"]
        message = request.form["message"].strip()
        msg = Message(subject="Conformer Webapp", body=f"Email: {email} \n \n{message}",
        sender=mail_username, recipients=[mail_username])
        mail.send(msg)
        return redirect(url_for("rdkit_page")) 
    

@app.route("/")
def index():
    return redirect(url_for("rdkit_page"))
    

@app.route("/rdkit")
def rdkit_page():
    return render_template("rdkit.html")


@app.route("/confab")
def confab_page():
    return render_template("confab.html")


@app.route("/reset/<method>/<job_id>")
def reset(method, job_id):
    try:
        shutil.rmtree(os.path.join(app.config["MOLECULE_UPLOADS"], job_id))
    except:
        pass #no job_id dir since user pressed download and then reset (download deletes job_id dir) 
    finally:
        if method == "confab":
            return redirect(url_for("confab_page"))
        else:
            return redirect(url_for("rdkit_page"))


@app.route("/<method>/<job_id>") 
def serve_files(method, job_id): 
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], job_id)):
        #check if a merged pdb/sdf file exists  
        match = [f for f in os.listdir() if f.startswith("ConformersMerged")]
        if match:
            name_file = match[0]
            mol_mem = io.BytesIO()
            with open(name_file, "rb") as fo:
                mol_mem.write(fo.read())
                mol_mem.seek(0)
            shutil.rmtree(os.path.join(app.config["MOLECULE_UPLOADS"], job_id))
            return send_file(mol_mem, as_attachment=True, attachment_filename=name_file, cache_timeout=0)
        else:
            zipfolder = zipfile.ZipFile("Conformers.zip", "w", zipfile.ZIP_STORED)
            for f in os.listdir():
                if f.startswith("conformer_"):
                    zipfolder.write(f)
            zipfolder.close()
            zip_mem = io.BytesIO()
            with open(zipfolder.filename, "rb") as fo:
                zip_mem.write(fo.read())
                zip_mem.seek(0)
            shutil.rmtree(os.path.join(app.config["MOLECULE_UPLOADS"], job_id))
            return send_file(zip_mem, mimetype="application/zip", as_attachment=True, 
            attachment_filename="Conformers.zip", cache_timeout=0)
    else:
        if method == "confab":
            return redirect(url_for("confab_page"))
        else:
            return redirect(url_for("rdkit_page"))


@app.route("/<method>", methods=["POST", "GET"])
def form_handler(method):
    if request.method == "POST":
        unique_id = str(uuid.uuid4())
        if method == "confab":
            force_field = request.form["force_field"]
        smiles = request.form["smiles_molecule"]
        mol_file = request.files["molecule_file"]
        no_conformers = int(request.form["no_conformers"])
        output_ext = request.form["output_ext"]
        output_seperate = request.form["output_seperate"]

        with open(os.path.join(app.config["BASE_DIR"], "molecules.txt"), "a") as f:
            f.write(time.strftime("%d/%m %H:%M:%S") + "\n")
            f.write(f"\t \t Smiles: {smiles} \n")
            f.write(f"\t \t PDB: {mol_file.filename} \n")
            f.write(f"\t \t N_conformers: {no_conformers} \n \n")

        if smiles:
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], unique_id)
            os.mkdir(mol_path)
        else:
            #else file was provided
            allowed_ext = ["pdb", "sdf"]
            extension = mol_file.filename.split(".")[-1] 
            assert extension in allowed_ext
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], unique_id)
            os.mkdir(mol_path)
            mol_file.save(os.path.join(mol_path, mol_file.filename))
            
            #use the NIH converter to get SMILES from PDB 
            if extension == "pdb":
                smiles = pdbToSmileConverter.pdb_to_smiles(os.path.join(mol_path, mol_file.filename))

        os.chdir(mol_path)

        if method == "confab":
            if smiles:
                mole = confab.generate_conformers(smiles, force_field=force_field)
            else:
                mole = confab.generate_conformers(mol_file.filename, force_field=force_field)
            if mole.NumConformers() > no_conformers:
                conf_sample = random.sample(range(mole.NumConformers()), no_conformers)
                confab.write_conformers(mole, conf_sample, output_ext, output_seperate)
            else:
                confab.write_conformers(mole, range(mole.NumConformers()), output_ext, output_seperate)

            return render_template("confab.html", method="confab", job_id=unique_id)

        #else method is rdkit
        else:
            if smiles:
                conformers = confgen_rdkit.gen_conformers(smiles, no_conformers=no_conformers)
            else:
                #sdf was provided (or other file ext, but only sdf for now)
                conformers = confgen_rdkit.gen_conformers(mol_file.filename, no_conformers=no_conformers)
            
            confgen_rdkit.write_confs_to_file(conformers, output_ext, output_seperate)

            return render_template("rdkit.html", method="rdkit", job_id=unique_id)
        

if __name__ == "__main__":
    app.run()

