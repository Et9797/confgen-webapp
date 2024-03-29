from . import app, make_celery, mail, Message
from ._rdkit import pdb_to_smiles, confgen
import os
import sys
from io import StringIO
from smtplib import SMTPRecipientsRefused


celery = make_celery(app)


@celery.task()
def generate_confs(
    smiles, 
    mol_filename,
    mol_path,
    no_conformers, 
    output_ext, 
    output_separate, 
    mail_address, 
    task_id
):
    # Capture stderr in case exception is thrown
    sio = sys.stderr = StringIO()
    exception_occurred = False
    try:
        # Generate conformers
        if smiles:
            conformers = confgen.generate_conformers(smiles, no_conformers=no_conformers)
        elif mol_filename.split(".")[-1] == 'pdb':
            smiles = pdb_to_smiles.convert(os.path.join(mol_path, mol_filename))
            conformers = confgen.generate_conformers(smiles, no_conformers=no_conformers)
        else:
            conformers = confgen.generate_conformers(os.path.join(mol_path, mol_filename), no_conformers=no_conformers)
    except:
        exception_occurred = True
        return sio.getvalue()
    else:
        # Write conformers to disk
        confgen.write_confs_to_file(conformers, mol_path, output_ext, output_separate)
    finally:
        if exception_occurred:
            app.logger.warn(sio.getvalue())
        
        # Send mail to user if one was provided
        if mail_address:
            msg = Message(
                subject = "Conformer generation job completed.", 
                body = (
                    "Your job has been completed. URL for the results page:\n"
                    f"http://confgen.net/results/{task_id}"
                ),
                sender = app.config["MAIL_USERNAME"], 
                recipients = [mail_address]
            )   
            try:
                mail.send(msg)
            except SMTPRecipientsRefused: 
                pass # Recipient email address may not be valid -> ignore