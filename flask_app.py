import io
import os
import zipfile
import traceback
import uuid
import confgen_rdkit
import PdbToSmileConverter
from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from flask_mail import Mail, Message
import logging

BASE_DIR = '/home/et/personal_projects/confgen-webapp/'
MOLECULE_UPLOADS = '/home/et/personal_projects/confgen-webapp/MOLECULE_UPLOADS/'

# Flask
app = Flask(__name__)
app.config["BASE_DIR"] = BASE_DIR
app.config["MOLECULE_UPLOADS"] = MOLECULE_UPLOADS
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0 # No caching of served conformer files

# Logging
app.logger.setLevel(logging.INFO)
fh = logging.FileHandler("appLog.log")
formatter = logging.Formatter('[%(asctime)s] - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%m-%Y %H:%M:%S')
fh.setFormatter(formatter)
app.logger.addHandler(fh)

# Contact
app.config["MAIL_SERVER"] = "smtp-mail.outlook.com"
app.config["MAIL_PORT"] = 587
app.config["MAIL_USE_TLS"] = True
app.config["MAIL_USE_SSL"] = False
app.config["MAIL_USERNAME"] = os.environ.get("EMAIL")
app.config["MAIL_PASSWORD"] = os.environ.get("EMAIL_PASS")
mail = Mail(app)

@app.route("/contact", methods=["POST","GET"])
def contact():
    if request.method == "POST":
        email = request.form["email"]
        message = request.form["message"].strip()
        msg = Message(subject = "Conformer Webapp", body = f"Email: {email} \n \n{message}",
                      sender = app.config["MAIL_USERNAME"], recipients = [app.config["MAIL_USERNAME"]]
                     )
        mail.send(msg)
        return redirect(url_for("index"))
    
@app.route("/")
def index():
    return render_template("index.html")

@app.route("/generate/<job_id>") 
def serve_files(job_id): 
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], job_id)):
        match = [f for f in os.listdir() if f.startswith("ConformersMerged")]
        if match:
            name_file = match[0]
            mol_mem = io.BytesIO()
            with open(name_file, "rb") as fo:
                mol_mem.write(fo.read())
                mol_mem.seek(0)
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
            return send_file(zip_mem, mimetype="application/zip", as_attachment=True, 
            attachment_filename="Conformers.zip", cache_timeout=0)

@app.route("/generate", methods=["POST", "GET"])
def form_handler():
    if request.method == "POST":
        unique_id = str(uuid.uuid4())
        smiles = request.form["SMILES"]
        mol_file = request.files["molFile"]
        no_conformers = int(request.form["noConfs"])
        output_ext = request.form["outputFormat"]
        try:
            output_seperate = request.form["seperateFiles"]
        except:
            output_seperate = "off"

        app.logger.info(f"ID: {unique_id}, SMILES: {smiles}, MolFile: {mol_file.filename}," 
                        f" N_conformers: {no_conformers}, Output: {output_ext}, Merged: {output_seperate}")

        if smiles:
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], unique_id)
            os.mkdir(mol_path)
        else:
            # Else a file was provided
            allowed_extensions = ["pdb", "sdf", "mol"]
            extension = mol_file.filename.split(".")[-1] 
            assert extension in allowed_extensions
            mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], unique_id)
            os.mkdir(mol_path)
            mol_file.save(os.path.join(mol_path, mol_file.filename))
            
            # Use the NIH converter to get SMILES from PDB 
            if extension == "pdb":
                try:
                    smiles = PdbToSmileConverter.pdb_to_smiles(os.path.join(mol_path, mol_file.filename))
                except Exception as e:
                    app.logger.error(traceback.format_exc())
                    return jsonify({'exception': str(e)}), 500

        os.chdir(mol_path)

        if smiles:
            try:
                conformers = confgen_rdkit.generate_conformers(smiles, no_conformers=no_conformers)
            except Exception as e:
                app.logger.error(traceback.format_exc())
                return jsonify({'exception': str(e)}), 500
        else:
            try:
                conformers = confgen_rdkit.generate_conformers(mol_file.filename, no_conformers=no_conformers)
            except Exception as e:
                app.logger.error(traceback.format_exc())
                return jsonify({'exception': str(e)}), 500

        confgen_rdkit.write_confs_to_file(conformers, output_ext, output_seperate)

        return jsonify({"job_id": unique_id})

        
if __name__ == "__main__":
    app.run()
