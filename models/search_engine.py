__author__ = 'pierleonia'
from plugin_haystack import Haystack,WhooshBackend
import os

index_bioentry = Haystack(biodb.bioentry,
                          backend=WhooshBackend,
                          indexdir= os.path.join(request.folder, 'databases'))
index_bioentry.indexes('name','description', 'accession', 'identifier' )
