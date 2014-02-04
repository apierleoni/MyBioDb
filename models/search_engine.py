from applications.MyBioDb.modules.search_engine import SolrBackend, ElasticSearchBackend

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
#
if settings.search_engine == 'elasticsearch':
    biodb_index = BioentrySearchEngine(biodb_handler,
                                       backend=ElasticSearchBackend)

elif settings.search_engine == 'solr':
    biodb_index = BioentrySearchEngine(biodb_handler,
                                   backend=SolrBackend,
                                   url="http://localhost:8983/solr",
                                   schema= os.path.join(request.folder, 'static', 'solr', 'schema.xml'))
elif settings.search_engine == 'whoosh':
    biodb_index = BioentrySearchEngine(biodb_handler,
                                        indexdir= os.path.join(request.folder, 'databases'))

if settings.search_engine:
    biodb_index.indexes()




