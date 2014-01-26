__author__ = 'pierleonia'

DEBUG=True

import os, traceback

class BioentrySearchEngineBackend(object):

    def rebuild(self, bioentry_ids, **kwargs):
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
    def map_to_index(self, seqrecord):
        pass


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
            from whoosh.fields import Schema, TEXT, ID, KEYWORD, NUMERIC
            from whoosh.analysis import StemmingAnalyzer
        except ImportError:
            raise ImportError("Cannot find Whoosh")
        self.indexname =".".join([self.biodb._uri_hash, 'whoosh'])
        try:
            self.ix = open_dir(self.indexdir, indexname=self.indexname)
        except:
            schema = Schema(id=ID(unique=True,stored=True),
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
            self.ix = create_in(self.indexdir, schema, indexname=self.indexname)

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

            rows = self.biodb(self.biodb.bioentry.id.belongs(bioentries[start:end])).select()
            #if DEBUG: print "round size found: ",len(rows)

            for row in rows:
                try:
                    #if DEBUG: print "Indexing bioentry %s - %i"%(row.name, i+1)
                    seqrecord = self.handler._retrieve_seqrecord(row.bioentry_id)
                    writer.update_document(id = unicode(row.id),
                                           db = unicode(self.biodb.biodatabase[row.biodatabase_id].name),
                                           **self.map_to_index(seqrecord))
                    #if DEBUG:
                    #    print "Indexed bioentry %s - %i"%(row.name, start)
                except:
                    if DEBUG:
                        print "error building index for id: ",row.id
                        traceback.print_exc()

            if len(rows)<m: break
            i+=1
        writer.commit()


    def map_to_index(self, seqrecord): #TODO: add taxonomy parsing

        def add_element(element, container):
            if isinstance(element,str):
                    container.append(element)
            elif isinstance(element,list):
                for i in element:
                    if isinstance(i, list):
                        container.append(i[0])
                    elif isinstance(i, str):
                        container.append(i)
            return container

        annotation_types, annotation_values = [],[]
        feature_types, feature_values = [],[]
        comments = []
        accessions = []
        keywords = []
        pubids, pubauths, pubtitles, pubjournals= [],[],[],[]

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

        kwargs = dict(name=unicode(seqrecord.name),
                    accession=unicode(" ".join(accessions)),
                    identifier=unicode(seqrecord.id),
                    description=unicode(seqrecord.description),
                    keyword=unicode(",".join(keywords)),
                    annotation=unicode( " ".join(annotation_values)),
                    annotationtype=unicode(",".join(annotation_types)),
                    comment=unicode(" ".join(comments)),
                    feature=unicode(" ".join(feature_values)),
                    featuretype=unicode(",".join(feature_types)),
                    lenght=unicode(len(seqrecord)),
                    dbxref=unicode(" ".join(seqrecord.dbxrefs)),
                    pubid=unicode(",".join(pubids)),
                    pubtitle=unicode(" ".join(pubtitles)),
                    pubauth=unicode(",".join(pubauths)),
                    pubjournal=unicode(",".join(pubjournals)))
        return kwargs

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