from gluon.storage import Storage
settings = Storage()

settings.migrate = True
settings.title = 'MyBioDb'
settings.subtitle = 'A Web UI to a BioSQL database'
settings.author = 'Andrea Pierleoni'
settings.author_email = 'apierleoni.dev@gmail.com'
settings.keywords = 'Web BioSQL web2py '
settings.description = 'A web-based user inteface to a BioSQL database'
settings.layout_theme = 'Default'
settings.database_uri = 'sqlite://storage.sqlite'
settings.security_key = 'dd2ad7b5-b721-4ea6-847e-13ef0ba6977a'
settings.email_server = 'localhost'
settings.email_sender = 'webbiosql@noreply.com'
settings.email_login = ''
settings.login_method = 'local'
settings.login_config = ''
settings.plugins = []
