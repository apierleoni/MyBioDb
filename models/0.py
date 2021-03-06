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


#MyBioDb specific settings
settings.biodb_conn_string = 'mysql://root:@127.0.0.1/mybiodb'
settings.timestamps_db = 'mysql://root:@127.0.0.1/mybiodb'
settings.search_engine = 'elasticsearch' # available choices: 'whoosh','solr', elasticsearch, or None
settings.solr_uri = 'http:/localhost:8983'
settings.solr_uri = 'http:/localhost:9200'
'''seconds to store a query in cache '''
settings.search_cache_type = 300
'''maximum number of results to report upon a search '''
settings.max_query_results=500000
'''maximum number of results to report upon a search '''
settings.max_dbname_search=1000
'''maximum number of results to report upon a search '''
settings.max_qualifier_search=10000
'''maximum number of entry loadable at one time'''
settings.max_entry_load = 1000000
