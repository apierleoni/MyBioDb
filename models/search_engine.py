from applications.MyBioDb.modules.search_engine import SolrBackend

__author__ = 'pierleonia'
from plugin_haystack import Haystack,WhooshBackend
from search_engine import BioentrySearchEngine
import os

#index_bioentry = Haystack(biodb.bioentry,
#                          backend=WhooshBackend,
#                          indexdir= os.path.join(request.folder, 'databases'))
#index_bioentry.indexes('name','description', 'accession', 'identifier' )
#
#index_bioentry_qualifier_value = Haystack(biodb.bioentry_qualifier_value,
#                          backend=WhooshBackend,
#                          indexdir= os.path.join(request.folder, 'databases'))
#index_bioentry.indexes('name','description', 'accession', 'identifier' )

# biodb_index = BioentrySearchEngine(biodb_handler,
#                                    indexdir= os.path.join(request.folder, 'databases'))

biodb_index = BioentrySearchEngine(biodb_handler,
                                   backend=SolrBackend,
                                   url="http://localhost:8983/solr",
                                   schema= os.path.join(request.folder, 'static', 'solr', 'schema.xml'))

biodb_index.indexes()
