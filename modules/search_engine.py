
__author__ = 'pierleonia'

DEBUG=True

import os, traceback
from multiprocessing import Pool

class BioentrySearchEngineBackend(object):

    def rebuild(self, bioentry_ids=[], **kwargs):
        raise NotImplementedError()
    def indexes(self, **kwargs):
        raise NotImplementedError()

    def after_insert(self, **kwargs):
        raise NotImplementedError()

    def after_update(self, **kwargs):
        raise NotImplementedError()

    def get_ids(self, **kwargs):
        raise NotImplementedError()

    def after_delete(self, **kwargs):
        raise NotImplementedError()

    def search(self, query, **kwargs):
        raise NotImplementedError()

    def quick_search(self, query):
        raise NotImplementedError()

    def create_loading_Pool(self):
        self.pool = Pool(processes=getCPUs())

    def add_bioentry_id_to_index(self, counter,bioentry_id):
        raise NotImplementedError()


    def map_to_index(self,handler, bioentry_id):

        def add_element(element, container):
            if isinstance(element,str):
                    container.append(unicode(element))
            elif isinstance(element,list):
                for i in element:
                    if isinstance(i, list):
                        container.append(unicode(i[0]))
                    elif isinstance(i, str):
                        container.append(unicode(i))
            return container

        seqrecord = handler._retrieve_seqrecord(bioentry_id)
        annotation_types, annotation_values = [],[]
        feature_types, feature_values = [],[]
        comments = []
        accessions = []
        keywords = []
        pubids, pubauths, pubtitles, pubjournals= [],[],[],[]
        taxonomy = [] #TODO: add taxonomy

        for k,v in seqrecord.annotations.items():
            if k == 'accessions':
                accessions = add_element(v,accessions)
            elif k.startswith('comment'):
                comments = add_element(v,comments)
            elif k == 'keywords':
                keywords = add_element(v,keywords)
            elif k=='references':
                if isinstance(v,list):
                    for ref in v:
                        pubids.append(ref.pubmed_id)
                        pubtitles.append(ref.title)
                        pubauths.append(ref.authors.strip())
                        pubjournals.append(ref.journal)

            else:
                annotation_values = add_element(v,annotation_values)
                annotation_types = add_element(k,annotation_types)
        for feature in seqrecord.features:
            feature_types.append(feature.type)
            for k,v in feature.qualifiers.items():
                feature_values = add_element(v,feature_values)

        kwargs = dict(id = unicode(bioentry_id),
                      db = unicode(handler.adaptor.biodatabase[handler.adaptor.bioentry[bioentry_id].biodatabase_id].name),
                      name=unicode(seqrecord.name),
                      accession=accessions,
                      identifier=unicode(seqrecord.id),
                      description=unicode(seqrecord.description),
                      keyword=keywords,
                      annotation=annotation_values,
                      annotationtype=annotation_types,
                      comment=comments,
                      feature=feature_values,
                      featuretype=feature_types,
                      lenght=unicode(len(seqrecord)),
                      dbxref=seqrecord.dbxrefs,
                      pubid=pubids,
                      pubtitle=pubtitles,
                      pubauth=pubauths,
                      pubjournal=pubjournals)
        return kwargs



class SearchEngineResult(object):
    def __init__(self, ids, handler):

        self.biodb = handler.adaptor
        self.db_query = self.biodb.bioentry._id.belongs(ids) # to be used in DAL queries
        self.selected_ids = ids
        self.count = len(ids)
        self.select_sql = self.biodb(self.biodb.bioentry.id.belongs(ids))._select() #raw sql to retrieve data from the bioentry table


def getCPUs():
    import multiprocessing
    try:
        return multiprocessing.cpu_count()
    except:
        return 1

def picklable_call(instance, name, args=(), kwargs=None):
    "indirect caller for instance methods and multiprocessing"
    if kwargs is None:
        kwargs = {}
    return getattr(instance, name)(*args, **kwargs)



class WhooshBackend(BioentrySearchEngineBackend):


    def __init__(self, handler, indexdir):
        self.handler = handler
        self.biodb = handler.adaptor
        self.indexdir = indexdir
        if not os.path.exists(indexdir):
            os.mkdir(indexdir)
    def indexes(self):
        try:
            from whoosh.index import create_in,open_dir
        except ImportError:
            raise ImportError("Cannot find Whoosh")
        self.indexname =".".join([self.biodb._uri_hash, 'whoosh'])
        try:
            self.ix = open_dir(self.indexdir, indexname=self.indexname)
        except:
            self.ix = create_in(self.indexdir, self._get_schema(), indexname=self.indexname)

    def rebuild(self, **kwargs):
        cpus = getCPUs()
        writer = self.ix.writer(procs=cpus, multisegment=True)
        bioentries =  kwargs.get('bientry_ids',[])
        if DEBUG: print "starting global index rebuilding with %i CPUs"%cpus
        if not bioentries:
            bioentries = [row.id for row in self.biodb(self.biodb.bioentry.id >0
                                                      ).select(self.biodb.bioentry.id)]
        if DEBUG: print "starting indexing of %i bioentries"%len(bioentries)
        #iterate over all bioentries at 100 max a time
        i, m = 0, 1000
        while True:
            start = i*m
            end = (i+1)*m
            if DEBUG:
                print "searching for round ",start,end
                #print "searching query: " + self.biodb(self.biodb.bioentry.id.belongs(bioentries[start:end]))._select()

            rows = self.biodb(self.biodb.bioentry.id.belongs(bioentries[start:end])).select(self.biodb.bioentry.id)
            #if DEBUG: print "round size found: ",len(rows)

            for row in rows:
                try:
                    #if DEBUG: print "Indexing bioentry %s - %i"%(row.name, i+1)
                    writer.update_document(**self.map_to_index(self.handler,row.id))
                    #if DEBUG:
                    #    print "Indexed bioentry %s - %i"%(row.name, start)
                except:
                    if DEBUG:
                        print "error building index for id: ",row.id
                        traceback.print_exc()

            if len(rows)<m: break
            i+=1
        writer.commit()



    def search(self, query, **kwargs):

        from whoosh.qparser import QueryParser,MultifieldParser

        fieldnames =  kwargs.pop('fieldnames', self.ix.schema.names())
        qp = MultifieldParser( fieldnames, schema=self.ix.schema)
        q = qp.parse(query)
        with self.ix.searcher() as s:
            results = s.search(q, **kwargs)
            if DEBUG: print "found %i hits in %.2fms"%(len(results.top_n),results.runtime*1000)
            ids = list(set(long(result['id']) for result in results))

        result = SearchEngineResult(ids, self.handler)
        return result
        #return ids

    def _get_schema(self):
        from whoosh.fields import Schema, TEXT, ID, KEYWORD, NUMERIC
        from whoosh.analysis import StemmingAnalyzer
        return Schema(id=ID(unique=True,stored=True),
                            db=ID(stored=True),
                            name=ID(stored=True),
                            accession=KEYWORD(scorable=True),
                            identifier=ID(stored=True),
                            description=TEXT(stored=True),
                            taxonomy=KEYWORD(lowercase=True,
                                            commas=True,
                                            scorable=True),
                            keyword=KEYWORD(lowercase=True,
                                            commas=True,
                                            scorable=True),
                            annotation=TEXT(analyzer=StemmingAnalyzer()),
                            annotationtype=KEYWORD(lowercase=True,
                                            scorable=True),
                            comment=TEXT(analyzer=StemmingAnalyzer()),
                            feature=TEXT(analyzer=StemmingAnalyzer()),
                            featuretype=KEYWORD(lowercase=True,
                                                commas=True,
                                                scorable=True),
                            lenght=NUMERIC(),
                            dbxref=KEYWORD(scorable=True),
                            pubid=KEYWORD(scorable=True),
                            pubtitle=TEXT(analyzer=StemmingAnalyzer()),
                            pubauth=KEYWORD(lowercase=True,
                                            commas=True,
                                            scorable=True),
                            pubjournal=KEYWORD(lowercase=True,
                                            commas=True,
                                            scorable=True),
                            )

    def map_to_index(self, handler, bioentry_id):
        documents = super(WhooshBackend, self).map_to_index(handler, bioentry_id)

        for k,v in documents.items():
            if isinstance(v, list):
                documents[k]=unicode(" ".join(v))
        return documents

    def quick_search(self, query, limit = 0):
        if limit > 0:
            return self.search(query, limit = limit,
                                 scored=True,
                                 fieldnames = ['accession',
                                               'description',
                                               'name'])
        else:
            return self.search(query,scored=True,
                                 fieldnames = ['accession',
                                               'description',
                                               'name'])


class SolrBackend(BioentrySearchEngineBackend):


    def __init__(self, handler, url="http://localhost:8983",schema=""):

        self.handler = handler
        self.biodb = handler.adaptor
        self.url = url
        if not schema:
            schema = self._get_default_schema()
        self.schemadoc = schema
        # if DEBUG: print schema

    def _get_default_schema(self):#TODO: update schema.xml to make strings not case sensitive
        from gluon import request

        return os.path.join(request.folder, 'databases', 'solr_schema.xml')


    def indexes(self):
        import sunburnt
        # import pysolr
        # if 1:
        try:
            self.interface = sunburnt.SolrInterface(self.url, self.schemadoc)
        except:
            raise RuntimeError("Cannot connect to Solr: %s" % self.url)

        # self.interface = pysolr.Solr(self.url)


    def rebuild(self, bioentry_ids=[], **kwargs):
        bioentries =  kwargs.get('bientry_ids',[])
        if DEBUG: print "starting global index rebuilding"
        if not bioentries:
            bioentries = [row.id for row in self.biodb(self.biodb.bioentry.id >0
                                                      ).select(self.biodb.bioentry.id)]
        if DEBUG: print "starting indexing of %i bioentries"%len(bioentries)
        #iterate over all bioentries at 100 max a time
        i, m = 0, 100
        while True:
            start = i*m
            end = (i+1)*m
            if DEBUG:
                print "searching for round ",start,end
            rows = self.biodb(self.biodb.bioentry.id.belongs(bioentries[start:end])).select(self.biodb.bioentry.id)
            documents = []
            for row in rows:
                try:
                    documents.append(self.map_to_index(self.handler,row.id))
                except:
                    if DEBUG:
                        print "error building index for id: ",row.id
                        traceback.print_exc()
            self.interface.add(documents)
            # self.interface.add_many(documents)
            if len(rows)<m: break
            i+=1

        self.interface.commit()


    def search(self, query, **kwargs):
        # results = self.interface.query(**fieldkeys).paginate(0,limit)
        # ids = [r['id'] for r in results]
        # return ids


        fieldnames =  kwargs.pop('fieldnames', self.interface.schema.fields.keys())
        search_all_fields =  kwargs.pop('search_all_fields', False)
        if search_all_fields:
            fieldnames = self.interface.schema.fields.keys()
        qd = dict()
        fields = self.interface.schema.fields
        for fname in fieldnames:
            field = fields[fname]
            if getattr(field, "class") == 'solr.StrField' :
                qd[fname] = query

            elif getattr(field, "class") == 'solr.TriIntField' :
                try:
                    qd[fname] = int(query)
                except:
                    pass


        results = self.interface.query(**qd).field_limit("id") .execute()#TODO: modify to get the OR by default
        if DEBUG: print "found %i hits in %.2fms"%(len(results),results.QTime)

        ids = list(set(long(result['id']) for result in results))

        result = SearchEngineResult(ids, self.handler)#TODO: return the stored data to avoid querying the db again if possible. use a Storage object and try to get the required fields, otherwise fallback to db query.
        return result

    def quick_search(self, query, limit = 0):
        if limit > 0:
            return self.search(query,rows=limit,
                                 fieldnames = ['accession',
                                               'description',
                                               'name'])
        else:
            return self.search(query,fieldnames = ['accession',
                                               'description',
                                               'name'])


    def map_to_index(self, handler, bioentry_id):
        document = super(SolrBackend, self).map_to_index(handler, bioentry_id)
        try:
            document['lenght'] = int(document['lenght'])
        except:
            pass
        return document


class ElasticSearchBackend(BioentrySearchEngineBackend):
    def __init__(self, handler, nodes = [], index_name = 'mybiodb'):
        self.handler = handler
        self.biodb = handler.adaptor
        self.nodes = nodes
        self.index_name = index_name

    def rebuild(self, bioentry_ids=[], **kwargs):
        # self.create_loading_Pool()


        bioentries =  kwargs.get('bientry_ids',[])
        if DEBUG: print "starting global index rebuilding"
        if not bioentries:
            bioentries = [row.id for row in self.biodb(self.biodb.bioentry.id >0
                                                      ).select(self.biodb.bioentry.id)]
        if DEBUG: print "starting indexing of %i bioentries"%len(bioentries)
        #iterate over all bioentries at 100 max a time
        # self.pool.apply_async(picklable_call, args = (self, 'add_bioentry_id_to_index',  zip(range(len(bioentries)),bioentries)))
        # self.pool.close()
        # self.pool.join()
        for i,bioentry_id in enumerate(bioentries):
            self.add_bioentry_id_to_index(i,bioentry_id)

    def add_bioentry_id_to_index(self, counter, bioentry_id):

        if counter%100 ==0 and DEBUG:
             print "\tadded %i bioentries to index"%counter

        try:
            self.interface.index(index=self.index_name,
                            doc_type="full_bioentry",
                            id=bioentry_id,
                            body=self.map_to_index(self.handler,bioentry_id)
                        )
        except:
            if DEBUG:
                print "error building index for id: ",bioentry_id
                traceback.print_exc()


    def search(self, query, **kwargs):
        if DEBUG:
            from datetime import datetime
            start_time = datetime.now()
        fieldnames =  kwargs.pop('fieldnames', "_all")

        results = self.interface.search(index=self.index_name, body={"query": {
                                                                            "query_string": {
                                                                                "query": query},
                                                                            "term": {
                                                                                "fields": fieldnames}

                                                                            }})


        if DEBUG:
            print "found %i hits in "%(results['_shards']['successful']), (datetime.now()-start_time)
        ids = []
        if results['_shards']['successful']:
            ids = [r['_id'] for r in results['hits']['hits']]

        return SearchEngineResult(ids, self.handler)


    def quick_search(self, query, limit = 0):
        if limit > 0:
            return self.search(query,
                               limit = limit,
                               fieldnames = ['accession',
                                               'description',
                                               'name'],
                               **{'from':0 })
        else:
            return self.search(query,fieldnames = ['accession',
                                               'description',
                                               'name'])

    def map_to_index(self, handler, bioentry_id):
        return super(ElasticSearchBackend, self).map_to_index(handler, bioentry_id)

    def indexes(self, **kwargs):
        import elasticsearch
        if self.nodes:
            try:
                self.interface = self.interface =  elasticsearch.Elasticsearch(self.nodes, **kwargs)
            except:
                raise RuntimeError("Cannot connect to ElasticSearch nodes: %s" % ", ".join(self.nodes))
        else:
            try:
                self.interface =  elasticsearch.Elasticsearch(**kwargs)
            except:
                raise RuntimeError("Cannot connect to ElasticSearch on localhost")




class BioentrySearchEngine(object):
    def __init__(self,handler, backend=WhooshBackend,**kwargs):
        self.handler = handler
        self.backend = backend(handler, **kwargs)
    def indexes(self):
        '''init indexes '''
        self.backend.indexes()
        #self.table._after_insert.append(
        #    lambda fields,id: self.backend.after_insert(fields,id))
        #self.table._after_update.append(
        #    lambda queryset,fields: self.backend.after_update(queryset,fields))
        #self.table._after_delete.append(
        #    lambda queryset: self.backend.after_delete(queryset))
    def rebuild(self, **kwargs):
        self.backend.rebuild( **kwargs)
    def search(self, query, **kwargs):
        return self.backend.search(query, **kwargs)
    def quick_search(self, query, limit = 0):
        return self.backend.quick_search(query, limit)