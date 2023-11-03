from . import app, mail, Message
from flask import render_template, redirect, url_for, request, send_file, jsonify
from .tasks import celery, generate_confs
import os
import io
import zipfile
import uuid


@app.route("/contact", methods=["POST"])
def contact():
    if request.method == "POST":
        email = request.form["email"]
        message = request.form["message"].strip()
        msg = Message(
            subject = "Conformer Webapp", 
            body = f"Email: {email}\n\n{message}",
            sender = app.config["MAIL_USERNAME"],
            recipients = [app.config["MAIL_USERNAME"]]
        )
        mail.send(msg)

        return ('', 204)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/generate", methods=["POST"])
def form_handler():
    if request.method == "POST":
        # Extract form data
        uniq_id = str(uuid.uuid4())
        smiles = request.form["SMILES"]
        mol_file = request.files["molFile"]
        no_conformers = int(request.form["noConfs"])
        output_ext = request.form["outputFormat"]
        mail_address = request.form["emailAddress"]
        output_separate = request.form.get("separateFiles", "off")

        # Log form data 
        app.logger.info(
            f"ID: {uniq_id}, SMILES: {smiles}, MolFile: {mol_file.filename}, " 
            f"N_conformers: {no_conformers}, Output: {output_ext}, "
            f"OutputSeparate: {output_separate}, E-mail: {mail_address}"
        )

        # Create folder to store conformers
        mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], uniq_id)
        os.mkdir(mol_path)
        if not smiles:
            # A file was provided
            allowed_extensions = ["pdb", "sdf", "mol"]
            extension = mol_file.filename.split(".")[-1]
            assert extension in allowed_extensions
            mol_file.save(os.path.join(mol_path, mol_file.filename))

        # Generate conformers
        task = generate_confs.apply_async(
            args = [
                smiles, mol_file.filename, mol_path, no_conformers, 
                output_ext, output_separate, mail_address, uniq_id
            ],
            task_id = uniq_id
        )
        
        return redirect(url_for("results", task_id=task.id))


@app.route("/results")
def results():
    task_id = request.args.get("task_id")
    if not task_id:
        return redirect(url_for("index"))
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], task_id)):
        return render_template("results.html")
    

@app.route("/task_status/<task_id>")
def task_status(task_id):
    status = celery.AsyncResult(task_id).state
    return jsonify({"state": status})


@app.route("/results/<task_id>")
def results_outcome(task_id):
    if os.path.exists(os.path.join(app.config["MOLECULE_UPLOADS"], task_id)):
        status = celery.AsyncResult(task_id).state
        info = celery.AsyncResult(task_id).info
        if info:
            # Task failed and threw an error. When task is successful info = None
            return render_template("failure.html", error_message=info)
        elif status == "SUCCESS":
            return render_template("success.html", task_id=task_id)
        else:
            return redirect(url_for("index"))
    else:
        return redirect(url_for("index"))


@app.route("/results/job/<task_id>") 
def serve_files(task_id): 
    mol_path = os.path.join(app.config["MOLECULE_UPLOADS"], task_id)
    if os.path.exists(mol_path):
        confs_merged = [f for f in os.listdir(mol_path) if f.startswith("ConformersMerged")]
        if confs_merged:
            file_name = confs_merged.pop()
            mol_buffer = io.BytesIO()
            with open(os.path.join(mol_path, file_name), "rb") as fo:
                mol_buffer.write(fo.read())
                mol_buffer.seek(0)
            return send_file(
                mol_buffer, 
                as_attachment=True,
                attachment_filename=file_name, 
                cache_timeout=0
            )
        else:
            zipfolder = zipfile.ZipFile(os.path.join(mol_path, "Conformers.zip"), "w", zipfile.ZIP_STORED)
            for f in os.listdir(mol_path):
                if f.startswith("conformer_"):
                    zipfolder.write(os.path.join(mol_path, f), f)
            zipfolder.close()
            zip_buffer = io.BytesIO()
            with open(os.path.join(mol_path, zipfolder.filename), "rb") as fo:
                zip_buffer.write(fo.read())
                zip_buffer.seek(0)
            return send_file(
                zip_buffer, 
                mimetype="application/zip", 
                as_attachment=True, 
                attachment_filename="Conformers.zip", 
                cache_timeout=0
            )
            