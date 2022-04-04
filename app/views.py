from app import app, mail, Message
from flask import render_template, redirect, url_for, request, send_file, jsonify
import io
from pathlib import Path
import zipfile
import uuid
from .tasks import celery, generate_confs

@app.route("/contact", methods=["POST"])
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

@app.route("/<uniq_id>") 
def serve_files(uniq_id): 
    mol_path = Path(app.config["MOLECULE_UPLOADS"], uniq_id)
    if Path(mol_path).exists():
        match = [f for f in Path(mol_path).iterdir() if f.startswith("ConformersMerged")]
        if match:
            name_file = match[0]
            mol_mem = io.BytesIO()
            with open(Path(mol_path, name_file), "rb") as fo:
                mol_mem.write(fo.read())
                mol_mem.seek(0)
            return send_file(mol_mem, as_attachment=True, attachment_filename=name_file, 
                             cache_timeout=0
                             )
        else:
            zipfolder = zipfile.ZipFile(Path(mol_path, "Conformers.zip"), 
                                        "w", zipfile.ZIP_STORED)
            for f in Path(mol_path).iterdir():
                if f.startswith("conformer_"):
                    zipfolder.write(Path(mol_path, f), f)
            zipfolder.close()
            zip_mem = io.BytesIO()
            with open(Path(mol_path, zipfolder.filename), "rb") as fo:
                zip_mem.write(fo.read())
                zip_mem.seek(0)
            return send_file(zip_mem, mimetype="application/zip", as_attachment=True, 
                             attachment_filename="Conformers.zip", cache_timeout=0
                             )
    return ('', 404)

@app.route("/generate", methods=["POST"])
def form_handler():
    if request.method == "POST":
        # Extract form data
        uniq_id = str(uuid.uuid4())
        smiles = request.form["SMILES"]
        mol_file = request.files["molFile"]
        no_conformers = int(request.form["noConfs"])
        output_ext = request.form["outputFormat"]
        email_address = request.form["emailAddress"]
        try:
            output_separate = request.form["separateFiles"]
        except:
            output_separate = "off"
        
        # Log form data 
        app.logger.info(f"ID: {uniq_id}, SMILES: {smiles}, MolFile: {mol_file.filename}," 
                        f" N_conformers: {no_conformers}, Output: {output_ext},"
                        f" OutputSeparate: {output_separate}, E-mail: {email_address}" 
                        )

        return redirect(url_for("results", task_id=uniq_id))
        # Create folder to store conformers
        # mol_path = Path(app.config["MOLECULE_UPLOADS"], uniq_id)
        # Path(mol_path).mkdir()
        # if not smiles:
        #     # A file was provided
        #     allowed_extensions = ["pdb", "sdf", "mol"]
        #     extension = mol_file.filename.split(".")[-1]
        #     assert extension in allowed_extensions
        #     mol_file.save(Path(mol_path, mol_file.filename))

        # # Generate conformers
        # task = generate_confs.apply_async(args = [smiles, mol_file.filename, mol_path,
        #                                           no_conformers, output_ext, output_separate],
        #                                   task_id = uniq_id
        #                                   )
        
        # #custom task id... nog kijken hoe file saven
        # # return redirect(url_for('results', task_id=task.id))
        # # return render template results html met vars uniq id en taskid
        # return jsonify({"uniq_id": uniq_id, "task_id": task.id})

@app.route("/results")
def results():
    args = request.args
    if not args:
        return ('', 404)
    
    # if Path(app.config["MOLECULE_UPLOADS"], args.get(task_id)).exists():
    #     pass

    # if pathtotask-id folder exists..
    # else return 404
    print(args.get("task_id"))
    # pass
    return render_template("results.html")
    # render template met vars args[task_id]
    # if uniq_id molpath exists..
    # if request.method == "GET":
    #     return render_template("results.html")
    # elif request.method == "POST":
    #     # serve_files code 
    #     pass
    
    
@app.route("/task_status/<task_id>")
def task_status(task_id):
    status = celery.AsyncResult(task_id).state
    return jsonify({"state": status})