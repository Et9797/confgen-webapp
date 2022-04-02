from flask import Flask
from flask_mail import Mail, Message
import logging
from celery import Celery

app = Flask(__name__)
if app.config["ENV"] == "development":
    app.config.from_object("config.Config")
else:
    app.config.from_object("config.ConfigProduction")

# Email setup
mail = Mail(app)

# Logging setup
app.logger.setLevel(logging.INFO)
fh = logging.FileHandler(app.config["FH"])
formatter = logging.Formatter('[%(asctime)s] - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%m-%Y %H:%M:%S')
fh.setFormatter(formatter)
app.logger.addHandler(fh)

# Celery setup
def make_celery(app):
    celery = Celery(
        app.import_name, 
        backend=app.config['CELERY_BACKEND_URL'],
        broker=app.config['CELERY_BROKER_URL'],
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    celery.finalize()
    return celery


from app import views