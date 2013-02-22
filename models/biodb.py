
'''Initialize biosql db '''
biodb_conn_string = 'sqlite://biodb.sqlite'
#biodb_conn_string = 'mysql://root:@127.0.0.1/mybiodb'
#biodb_conn_string = 'postgres://mybiodb:mypass@127.0.0.1/mybiodb'
timestamps_db = 'sqlite://biodb.sqlite'
biodb_handler = BioSQLHandler(biodb_conn_string, compatibility_mode = False, time_stamps = timestamps_db)
biodb = biodb_handler.adaptor
if biodb_handler._build_error:
    response.flash = 'DB model building error'
    
    
    
'''Initial configuration data '''
biodatabases = [row.name for row in biodb(biodb.biodatabase.id>0).select()]
if 'UniProt' not in biodatabases:
    biodb_handler.make_new_db(dbname = 'UniProt',
                             description = 'Entries loaded from UniProt',)
