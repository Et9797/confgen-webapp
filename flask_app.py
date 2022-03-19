from enum import unique
import io
import os
import zipfile
import traceback
import uuid
import confgen_rdkit
import pdb_to_smiles
from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from flask_mail import Mail, Message
import logging
import tasks

BASE_DIR = '/home/et/personal_projects/confgen-webapp/'
MOLECULE_UPLOADS = '/home/et/personal_projects/confgen-webapp/MOLECULE_UPLOADS/'

# Flask & Celery setup
app = Flask(__name__)
app.config["BASE_DIR"] = BASE_DIR
app.config["MOLECULE_UPLOADS"] = MOLECULE_UPLOADS
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0 # No caching of served conformer files
app.config.update(
    CELERY_BACKEND_URL='redis://localhost:6379/0',
    CELERY_BROKER_URL='redis://localhost:6379',
    CELERY_ACCEPT_CONTENT = ['pickle', 'application/json']
)
celery = tasks.make_celery(app)

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
        return ('', 204)
    
@app.route("/")
def index():
    return render_template("index.html")

@app.route("/<job_id>") 
def serve_files(job_id): 
    mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], job_id)
    if os.path.exists(mol_path):
        match = [f for f in os.listdir(mol_path) if f.startswith("ConformersMerged")]
        if match:
            name_file = match[0]
            mol_mem = io.BytesIO()
            with open(os.path.join(mol_path, name_file), "rb") as fo:
                mol_mem.write(fo.read())
                mol_mem.seek(0)
            return send_file(mol_mem, as_attachment=True, attachment_filename=name_file, cache_timeout=0)
        else:
            zipfolder = zipfile.ZipFile(os.path.join(mol_path, "Conformers.zip"), "w", zipfile.ZIP_STORED)
            for f in os.listdir(mol_path):
                if f.startswith("conformer_"):
                    zipfolder.write(os.path.join(mol_path, f), f)
            zipfolder.close()
            zip_mem = io.BytesIO()
            with open(os.path.join(mol_path, zipfolder.filename), "rb") as fo:
                zip_mem.write(fo.read())
                zip_mem.seek(0)
            return send_file(zip_mem, mimetype="application/zip", as_attachment=True, 
                             attachment_filename="Conformers.zip", cache_timeout=0)
    else:
        return ('', 404)

@app.route("/generate", methods=["POST", "GET"])
def form_handler():
    if request.method == "POST":
        # Extract form data
        unique_id = str(uuid.uuid4())
        smiles = request.form["SMILES"]
        mol_file = request.files["molFile"]
        no_conformers = int(request.form["noConfs"])
        output_ext = request.form["outputFormat"]
        try:
            output_separate = request.form["separateFiles"]
        except:
            output_separate = "off"

        # Log form data 
        app.logger.info(f"ID: {unique_id}, SMILES: {smiles}, MolFile: {mol_file.filename}," 
                        f" N_conformers: {no_conformers}, Output: {output_ext}, Merged: {output_separate}")

        # Create folder to store conformers
        mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], unique_id)
        os.mkdir(mol_path)
        if not smiles:
            # A file was provided
            allowed_extensions = ["pdb", "sdf", "mol"]
            extension = mol_file.filename.split(".")[-1] 
            assert extension in allowed_extensions
            mol_file.save(os.path.join(mol_path, mol_file.filename))
            
            # Assign bond orders if ligand is provided as PDB
            try:
                smiles = pdb_to_smiles.convert(os.path.join(mol_path, mol_file.filename))
            except Exception as e:
                app.logger.error(traceback.format_exc())
                return render_template("index.html", error="NIH")
        
        # Generate conformers
        task = gen_confs.delay(smiles, mol_file.filename, mol_path, no_conformers, 
                        output_ext, output_separate)
        
        return jsonify({"uniq_id": unique_id, "task_id": task.id})

@celery.task()
def gen_confs(smiles, mol_filename, mol_path, no_conformers, output_ext, output_seperate):
    if smiles:
        conformers = confgen_rdkit.generate_conformers(smiles, no_conformers=no_conformers)
    else:
        conformers = confgen_rdkit.generate_conformers(os.path.join(mol_path, mol_filename), 
                                                       no_conformers=no_conformers)
    confgen_rdkit.write_confs_to_file(conformers, mol_path, output_ext, output_seperate)

@app.route("/task_status/<task_id>", methods=["POST", "GET"])
def task_status(task_id):
    if request.method == "POST":
        status = celery.AsyncResult(task_id).state
        return jsonify({"state": status})

#handling error ook met polling

if __name__ == "__main__":
    app.run(debug=True)
