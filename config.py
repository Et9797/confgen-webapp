import os

# Try to load email credentials (only available in production)
try:
    import _email
    _EMAIL = getattr(_email, 'EMAIL', None)
    _EMAIL_PASS = getattr(_email, 'EMAIL_PASS', None)
except (ImportError, AttributeError):
    _EMAIL = None
    _EMAIL_PASS = None


class Config:
    """Development config - no email required"""
    BASE_DIR = '/home/et/personal-projects/confgen-webapp/'
    MOLECULE_UPLOADS = '/home/et/personal-projects/confgen-webapp/app/MOLECULE_UPLOADS/'
    FH = '/home/et/personal-projects/confgen-webapp/confgen.log'

    SEND_FILE_MAX_AGE_DEFAULT = 0

    # Mail disabled in development
    MAIL_SERVER = "smtp.gmail.com"
    MAIL_PORT = 587
    MAIL_USE_TLS = True
    MAIL_USE_SSL = False
    MAIL_USERNAME = None
    MAIL_PASSWORD = None
    MAIL_SUPPRESS_SEND = True

    CELERY_BACKEND_URL = "redis://localhost:6379/0"
    CELERY_BROKER_URL = "redis://localhost:6379/0"


class ProductionConfig(Config):
    """Production config - uses email credentials if available"""
    BASE_DIR = '/var/www/html/confgen-webapp/'
    MOLECULE_UPLOADS = '/var/www/html/confgen-webapp/app/MOLECULE_UPLOADS/'
    FH = "/var/www/html/confgen-webapp/confgen.log"

    # Email enabled in production (if credentials exist)
    MAIL_USERNAME = _EMAIL
    MAIL_PASSWORD = _EMAIL_PASS
    MAIL_SUPPRESS_SEND = False
