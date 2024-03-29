import os
import _email


class Config:
    BASE_DIR = '/home/et/personal-projects/confgen-webapp/' 
    MOLECULE_UPLOADS = '/home/et/personal-projects/confgen-webapp/app/MOLECULE_UPLOADS/' 
    FH = '/home/et/personal-projects/confgen-webapp/confgen.log'

    SEND_FILE_MAX_AGE_DEFAULT = 0
    
    MAIL_SERVER = "smtp.gmail.com"
    MAIL_PORT = 587
    MAIL_USE_TLS = True
    MAIL_USE_SSL = False
    MAIL_USERNAME = _email.EMAIL
    MAIL_PASSWORD = _email.EMAIL_PASS
    
    CELERY_BACKEND_URL = "redis://localhost:6379/0"
    CELERY_BROKER_URL = "redis://localhost:6379/0"
    
class ProductionConfig(Config):
    BASE_DIR = '/var/www/html/confgen-webapp/'
    MOLECULE_UPLOADS = '/var/www/html/confgen-webapp/app/MOLECULE_UPLOADS/'
    FH = "/var/www/html/confgen-webapp/confgen.log"
