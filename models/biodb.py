
'''Initialize biosql db '''
biodb_conn_string = 'sqlite://biodb.sqlite'
biodb_handler = BioSQLHandler(biodb_conn_string, compatibility_mode = False, time_stamps = biodb_conn_string)
biodb = biodb_handler.adaptor
if biodb_handler._build_error:
    response.flash = 'DB model building error'
    
    
    
'''Initial configuration data '''
biodatabases = [row.name for row in biodb(biodb.biodatabase.id>0).select()]
if 'UniProt' not in biodatabases:
    biodb_handler.make_new_db(dbname = 'UniProt',
                             description = 'Entries loaded from UniProt',)
