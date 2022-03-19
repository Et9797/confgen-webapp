from celery import Celery

def make_celery(app):
    celery = Celery(
        app.import_name, 
        backend='redis://localhost:6379/0',
        broker='redis://localhost:6379/0',
        accept_content=['pickle', 'application/json']
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    celery.finalize()
    return celery