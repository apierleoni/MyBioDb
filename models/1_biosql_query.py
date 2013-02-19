

class BioSQLQuery(object): 
    '''
    template class
    
    These are the so-called �rich comparison� methods, and are called for comparison operators in preference to __cmp__() below. 
    The correspondence between operator symbols and method names is as follows: x<y calls x.__lt__(y), x<=y calls x.__le__(y), 
    x==y calls x.__eq__(y), x!=y and x<>y call x.__ne__(y), x>y calls x.__gt__(y), and x>=y calls x.__ge__(y).'''
    
    def __init__(self, handler, biodb = None):
        self.handler = handler
        self.adaptor = handler.adaptor
        self.biodb = biodb
        self.query = None
        if self.biodb:
            self.dbid = self.handler.dbs[biodb]
            self.biodbquery = (self.adaptor.bioentry.biodatabase_id == self.dbid)
        else:
            self.biodbquery = (self.adaptor.bioentry.biodatabase_id > 0)
    
    def __lt__(self, other):
        raise NotImplementedError('This method is not available')
    def __le__(self, other):
        raise NotImplementedError('This method is not available')
    def __eq__(self, other):
        raise NotImplementedError('This method is not available')
    def __ne__(self, other):
        raise NotImplementedError('This method is not available')
    def __gt__(self, other):
        raise NotImplementedError('This method is not available')
    def __ge__(self, other):
        raise NotImplementedError('This method is not available')
    def contains(self, other):
        raise NotImplementedError('This method is not available')
    def startswith(self, other):
        raise NotImplementedError('This method is not available')
    def like(self, other):
        raise NotImplementedError('This method is not available')
    
    def select(self, mode = 'list', fields = ['bioentry_id']):
        if fields == ['bioentry_id']:
            ids = self.adaptor(self.biodbquery & self.query).select(self.adaptor.bioentry.bioentry_id,
                                                                    groupby = self.adaptor.bioentry.bioentry_id, 
                                                                    cache=(cache.ram, CacheTimes.search))
        results = []
        try:
            results = [row.id for row in ids]
        except:
            results = [row.bioentry.id for row in ids]
        
        if mode == 'set':
            return set(results)
        else:
            return results
        
    def count(self,):
        return self.adaptor(self.biodbquery & self.query).count()

        
    '''
    def __and__(self, other):
        return Query('(%s AND %s)' % (self, other))

    def __or__(self, other):
        return Query('(%s OR %s)' % (self, other))

    def __invert__(self):
        return Query('(NOT %s)' % self)

    def __str__(self):
        return self.query
    '''
   
   
class BioSQLQueryQuick(BioSQLQuery):
    '''Handles search by name'''
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry.name == other) | \
                                  (self.adaptor.bioentry.accession == other) | \
                                  (self.adaptor.bioentry.description == other) | \
                                  (self.adaptor.bioentry.identifier == other) | \
                                  (self.adaptor.bioentry.bioentry_id == other))
                else:
                    localquery = localquery | ((self.adaptor.bioentry.name == other) | \
                                              (self.adaptor.bioentry.accession == other) | \
                                              (self.adaptor.bioentry.description == other) | \
                                              (self.adaptor.bioentry.identifier == other) | \
                                              (self.adaptor.bioentry.bioentry_id == other))
        else:
            localquery = ((self.adaptor.bioentry.name == other) | \
                          (self.adaptor.bioentry.accession == other) | \
                          (self.adaptor.bioentry.description == other) | \
                          (self.adaptor.bioentry.identifier == other) | \
                          (self.adaptor.bioentry.bioentry_id == other))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    
    def contains(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry.name.lower().contains(element.lower())) | \
                                  (self.adaptor.bioentry.accession.lower().contains(element.lower())) | \
                                  (self.adaptor.bioentry.description.lower().contains(element.lower())) | \
                                  (self.adaptor.bioentry.identifier.lower().contains(element.lower()))  )
                else:
                    localquery = localquery | ((self.adaptor.bioentry.name.lower().contains(element.lower())) | \
                                              (self.adaptor.bioentry.accession.lower().contains(element.lower())) | \
                                              (self.adaptor.bioentry.description.lower().contains(element.lower())) | \
                                              (self.adaptor.bioentry.identifier.lower().contains(element.lower()))  )
        else:
            localquery = ((self.adaptor.bioentry.name.lower().contains(other.lower())) | \
                          (self.adaptor.bioentry.accession.lower().contains(other.lower())) | \
                          (self.adaptor.bioentry.description.lower().contains(other.lower())) | \
                          (self.adaptor.bioentry.identifier.lower().contains(other.lower())))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry.name.lower().startswith(element.lower())) | \
                                  (self.adaptor.bioentry.accession.lower().startswith(element.lower())) | \
                                  (self.adaptor.bioentry.description.lower().startswith(element.lower())) | \
                                  (self.adaptor.bioentry.identifier.lower().startswith(element.lower())) )
                else:
                    localquery = localquery | ((self.adaptor.bioentry.name.lower().startswith(element.lower())) | \
                                              (self.adaptor.bioentry.accession.lower().startswith(element.lower())) | \
                                              (self.adaptor.bioentry.description.lower().startswith(element.lower())) | \
                                              (self.adaptor.bioentry.identifier.lower().startswith(element.lower())))
        else:
            localquery = ((self.adaptor.bioentry.name.lower().startswith(other.lower())) | \
                          (self.adaptor.bioentry.accession.lower().startswith(other.lower())) | \
                          (self.adaptor.bioentry.description.lower().startswith(other.lower())) | \
                          (self.adaptor.bioentry.identifier.lower().startswith(other.lower())))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()    
   
class BioSQLQueryAccession(BioSQLQuery):
    '''Handles search by accession'''
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.accession == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.accession == element)
        else:
            localquery = (self.adaptor.bioentry.accession == other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.accession != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.accession != element)
        else:
            localquery = (self.adaptor.bioentry.accession != other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.accession.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.accession.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.accession.lower().contains(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.accession.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.accession.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.accession.lower().startswith(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()     
    
    
    
class BioSQLQueryName(BioSQLQuery):
    '''Handles search by name'''
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.name == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.name == element)
        else:
            localquery = (self.adaptor.bioentry.name == other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.name != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.name != element)
        else:
            localquery = (self.adaptor.bioentry.name != other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.name.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.name.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.name.lower().contains(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.name.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.name.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.name.lower().startswith(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()     

    
class BioSQLQueryID(BioSQLQuery):
    '''Handles search by identifier'''
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.identifier == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.identifier == element)
        else:
            localquery = (self.adaptor.bioentry.identifier == other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.identifier != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.identifier != element)
        else:
            localquery = (self.adaptor.bioentry.identifier != other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.identifier.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.identifier.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.identifier.lower().contains(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.identifier.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.identifier.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.identifier.lower().startswith(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()  
    
class BioSQLQueryDescription(BioSQLQuery):
    '''Handles search by description'''   
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.description == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.description == element)
        else:
            localquery = (self.adaptor.bioentry.description == other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.description != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.description != element)
        else:
            localquery = (self.adaptor.bioentry.description != other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.description.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.description.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.description.lower().contains(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.description.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry.description.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.bioentry.description.lower().startswith(other.lower()))
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()  
    
class BioSQLQueryBioentryID(BioSQLQuery):
    '''Handles search by bioentry_id'''
    def __eq__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id == element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id == other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id != element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id != other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    
    def __lt__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id < element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id < element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id < other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

    def __le__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id <= element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id <= element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id <= other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

    def __gt__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id > element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id > element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id > other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ge__(self, other):
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry.bioentry_id >= element)
                else:
                    localquery = localquery | (self.adaptor.bioentry.bioentry_id >= element)
        else:
            localquery = (self.adaptor.bioentry.bioentry_id >= other)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

class BioSQLQueryAnnotationType(BioSQLQuery):
    '''Handles search by annotation type'''
    def __eq__(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name == element)
                else:
                    localquery = localquery | (self.adaptor.term.name == element)
        else:
            localquery = (self.adaptor.term.name == other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
                    
        return self.count()
    
    def __ne__(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name != element)
                else:
                    localquery = localquery | (self.adaptor.term.name != element)
        else:
            localquery = (self.adaptor.term.name != other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery                   
        return self.count()
    
    def contains(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.term.name.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.term.name.lower().contains(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery 
        return self.count()
        
    def startswith(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.term.name.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.term.name.lower().startswith(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery 
        return self.count()
    
class BioSQLQueryAnnotationValue(BioSQLQuery):
    '''Handles search in values of every annotation field'''
    def __eq__(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry_qualifier_value.value == element)
                else:
                    localquery = localquery | (self.adaptor.bioentry_qualifier_value.value == element)
        else:
            localquery = (self.adaptor.bioentry_qualifier_value.value == other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry_qualifier_value.value != element)
                else:
                    localquery = localquery | (self.adaptor.bioentry_qualifier_value.value != element)
        else:
            localquery = (self.adaptor.bioentry_qualifier_value.value != other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry_qualifier_value.value.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry_qualifier_value.value.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.bioentry_qualifier_value.value.lower().contains(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        prequery = (self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.bioentry_qualifier_value.value.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.bioentry_qualifier_value.value.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.bioentry_qualifier_value.value.lower().startswith(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    
class BioSQLQueryAnnotationTypeValue(BioSQLQuery):
    '''Handles search in values a specific annotation field
    keys and values must be passed as a tuple of (key, value)
    all operators refers to the value, and the key must be provided
    exact (case-sensitive) since it is always set as =='''
    def __eq__(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                key, value = element
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value == value) & \
                                  (self.adaptor.term.name == key))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value == value) & \
                                               (self.adaptor.term.name == key))
        else:
            key, value = other
            localquery = ((self.adaptor.bioentry_qualifier_value.value == value) & \
                          (self.adaptor.term.name == key))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                key, value = element
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value != value) & \
                                  (self.adaptor.term.name == key))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value != value) & \
                                               (self.adaptor.term.name == key))
        else:
            key, value = other
            localquery = ((self.adaptor.bioentry_qualifier_value.value == value) & \
                          (self.adaptor.term.name == key))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                key, value = element
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().contains(value.lower())) & \
                                  (self.adaptor.term.name == key))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value.lower().contains(value.lower())) & \
                                               (self.adaptor.term.name == key))
        else:
            key, value = other
            localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().contains(value.lower())) & \
                          (self.adaptor.term.name == key))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                key, value = element
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(value.lower())) & \
                                  (self.adaptor.term.name == key))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(value.lower())) & \
                                               (self.adaptor.term.name == key))
        else:
            key, value = other
            localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(value.lower())) & \
                          (self.adaptor.term.name == key))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    
    
class BioSQLQueryFeatureType(BioSQLQuery):
    '''Handles search by feature type'''
    def __eq__(self, other):
        prequery = (self.adaptor.seqfeature.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name == element) 
                else:
                    localquery = localquery | (self.adaptor.term.name == element) 
        else:
            localquery = (self.adaptor.term.name == other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        prequery = (self.adaptor.seqfeature.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name != element) 
                else:
                    localquery = localquery | (self.adaptor.term.name != element) 
        else:
            localquery = (self.adaptor.term.name != other) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = (self.adaptor.seqfeature.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name.lower().contains(element.lower())) 
                else:
                    localquery = localquery | (self.adaptor.term.name.lower().contains(element.lower())) 
        else:
            localquery = (self.adaptor.term.name.lower().contains(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        prequery = (self.adaptor.seqfeature.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                   (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.term.name.lower().startswith(element.lower())) 
                else:
                    localquery = localquery | (self.adaptor.term.name.lower().startswith(element.lower())) 
        else:
            localquery = (self.adaptor.term.name.lower().startswith(other.lower())) 
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
class BioSQLQueryDBXref(BioSQLQuery):
    '''Handles search by dbxref'''
    def __eq__(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                   (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery =  (self.adaptor.dbxref.accession == element)
                else:
                    localquery = localquery | (self.adaptor.dbxref.accession == element)
        else:
            localquery = (self.adaptor.dbxref.accession == other)
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

    
    def __ne__(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                   (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery =  (self.adaptor.dbxref.accession != element)
                else:
                    localquery = localquery | (self.adaptor.dbxref.accession != element)
        else:
            localquery = (self.adaptor.dbxref.accession != other)
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                   (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery =  (self.adaptor.dbxref.accession.lower().contains(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.dbxref.accession.lower().contains(element.lower()))
        else:
            localquery = (self.adaptor.dbxref.accession.lower().contains(other.lower()))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()


        
    def startswith(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                   (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery =  (self.adaptor.dbxref.accession.lower().startswith(element.lower()))
                else:
                    localquery = localquery | (self.adaptor.dbxref.accession.lower().startswith(element.lower()))
        else:
            localquery = (self.adaptor.dbxref.accession.lower().startswith(other.lower()))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
class BioSQLQueryKeyword(BioSQLQuery):
    '''Handles search by keyword annotation value
    special type of type_value search, added for commodity
    searches for a value in all the qualifiers whose term starts with keyword'''

    def __eq__(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value == element) & \
                                  (self.adaptor.term.name.lower().startswith('keyword')))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value == element) & \
                                               (self.adaptor.term.name.lower().startswith('keyword')))
        else:
            localquery = ((self.adaptor.bioentry_qualifier_value.value == other) & \
                          (self.adaptor.term.name.lower().startswith('keyword')))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def __ne__(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value != element) & \
                                  (self.adaptor.term.name.lower().startswith('keyword')))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value != element) & \
                                               (self.adaptor.term.name.lower().startswith('keyword')))
        else:
            localquery = ((self.adaptor.bioentry_qualifier_value.value == other) & \
                          (self.adaptor.term.name.lower().startswith('keyword')))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().contains(element.lower())) & \
                                  (self.adaptor.term.name.lower().startswith('keyword')))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value.lower().contains(element.lower())) & \
                                               (self.adaptor.term.name.lower().startswith('keyword')))
        else:
            localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().contains(other.lower())) & \
                          (self.adaptor.term.name.lower().startswith('keyword')))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
        
    def startswith(self, other):
        prequery = ((self.adaptor.bioentry_qualifier_value.bioentry_id == self.adaptor.bioentry.bioentry_id) & \
                    (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id))
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(element.lower())) & \
                                  (self.adaptor.term.name.lower().startswith('keyword')))
                else:
                    localquery = localquery | ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(element.lower())) & \
                                               (self.adaptor.term.name.lower().startswith('keyword')))
        else:
            localquery = ((self.adaptor.bioentry_qualifier_value.value.lower().startswith(other.lower())) & \
                          (self.adaptor.term.name.lower().startswith('keyword')))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()


class BioSQLQuerySeq(BioSQLQuery):
    '''Handles simple text searches in sequence '''
    def __eq__(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.biosequence.bioentry_id)
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.biosequence.seq == element.upper())
                else:
                    localquery = localquery | (self.adaptor.biosequence.seq == element.upper())
        else:
            localquery = (self.adaptor.biosequence.seq == other.upper())
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

    
    def __ne__(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.biosequence.bioentry_id)
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.biosequence.seq != element.upper())
                else:
                    localquery = localquery | (self.adaptor.biosequence.seq != element.upper())
        else:
            localquery = (self.adaptor.biosequence.seq != other.upper())
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    def contains(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.biosequence.bioentry_id)
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.biosequence.seq.contains(element.upper()))
                else:
                    localquery = localquery | (self.adaptor.biosequence.seq.contains(element.upper()))
        else:
            localquery = (self.adaptor.biosequence.seq.contains(other.upper()))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()

        
    def startswith(self, other):
        prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.biosequence.bioentry_id)
        if isinstance(other, list):
            for i, element in enumerate(other):
                if not i:
                    localquery = (self.adaptor.biosequence.seq.startswith(element.upper()))
                else:
                    localquery = localquery | (self.adaptor.biosequence.seq.startswith(element.upper()))
        else:
            localquery = (self.adaptor.biosequence.seq.startswith(other.upper()))
        localquery = (prequery & localquery)
        if self.query:
            self.query = self.query | localquery    
        else:
            self.query = localquery
        return self.count()
    
    
class BioSQLQueryCreatedBy(BioSQLQuery):
    '''Handles search by creator user id '''
    def __eq__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_by == element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_by == element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_by == other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
        
    def __ne__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_by != element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_by != element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_by != other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
    
class BioSQLQueryModifiedBy(BioSQLQuery):
    '''Handles search by last modification user id '''
    def __eq__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_by == element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_by == element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_by == other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')

    def __ne__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_by != element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_by != element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_by != other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')


class BioSQLQueryCreatedOn(BioSQLQuery):
    ''' '''
    def __eq__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on == element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on == element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on == other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
            
    def __ne__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on != element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on != element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on != other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
            
    def __lt__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on < element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on < element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on < other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
    def __le__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on <= element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on <= element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on <= other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')    
    def __gt__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on > element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on > element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on > other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
    def __ge__(self, other):       
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.created_on >= element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.created_on >= element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.created_on >= other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
            
            
class BioSQLQueryModifiedOn(BioSQLQuery):
    ''' '''
    def __eq__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on == element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on == element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on == other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
            
    def __ne__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on != element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on != element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on != other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
            
    def __lt__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on < element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on < element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on < other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
    def __le__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on <= element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on <= element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on <= other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')    
    def __gt__(self, other):
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on > element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on > element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on > other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')
    def __ge__(self, other):       
        if self.handler.time_stamps_on_biosql:
            prequery = (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_timestamp.bioentry_id)
            if isinstance(other, list):
                for i, element in enumerate(other):
                    if not i:
                        localquery = (self.adaptor.bioentry_timestamp.modified_on >= element)
                    else:
                        localquery = localquery | (self.adaptor.bioentry_timestamp.modified_on >= element)
            else:
                localquery = (self.adaptor.bioentry_timestamp.modified_on >= other)
            localquery = (prequery & localquery)
            if self.query:
                self.query = self.query | localquery    
            else:
                self.query = localquery
            return self.count()
        else:
            raise NotImplementedError('''This method is not available for timestamps table in a different db. 
Please use a simple query over the corresponding "bioentry_timestamp" table''')

'''TODO
BioSQLQueryQuick ---> single query solo con contains on bioentry table in name, description, identifier and accession
BioSQLQueryFullText ---> query multiple su :
                                            - bioentry table
                                            - bioentry qualifiers values
                                            - seqfeature qualifeirs values
                                            - dbxrefs
                                            - sequence ??? potrebbe essere molto pesante
BioSQLQueryReferencePubmedID ---> just __eq__ and __ne__
BioSQLQueryReferenceTitle --->
BioSQLQueryReferenceAuthors --->
BioSQLQueryTaxon ---> ritornare tutte le bientry facendo un check dell'albero di taxon, o semplicemente sull'ultima foglia??

        elif request.vars.type=='pmid':
            results=abdb_handler.adaptor(biosequence_join &
                             biodatabase_query&\
                          (abdb_handler.adaptor.bioentry_reference.bioentry_id==abdb_handler.adaptor.bioentry.bioentry_id) & \
                          (abdb_handler.adaptor.bioentry_reference.reference_id==abdb_handler.adaptor.reference.reference_id) & \
                          (abdb_handler.adaptor.reference.dbxref_id==abdb_handler.adaptor.dbxref.dbxref_id) & \
                          (abdb_handler.adaptor.dbxref.dbname.lower()=='pubmed') & \
                          (abdb_handler.adaptor.dbxref.accession==request.vars.query)
                          ).select(abdb_handler.adaptor.bioentry.bioentry_id,
                                   abdb_handler.adaptor.biodatabase.name,
                                   abdb_handler.adaptor.bioentry.accession,
                                   abdb_handler.adaptor.bioentry.name,
                                   abdb_handler.adaptor.bioentry.description,
                                   abdb_handler.adaptor.biosequence.length,
                                   cache=(cache.ram, CacheTimes.search))


        elif request.vars.type=='titletext':
            results=abdb_handler.adaptor(biosequence_join &
                             biodatabase_query&\
                          (abdb_handler.adaptor.bioentry_reference.bioentry_id==abdb_handler.adaptor.bioentry.bioentry_id) & \
                          (abdb_handler.adaptor.bioentry_reference.reference_id==abdb_handler.adaptor.reference.reference_id) & \
                          (abdb_handler.adaptor.reference.title.lower().like('%%%s%%'%request.vars.query.lower()))                              
                          ).select(abdb_handler.adaptor.bioentry.bioentry_id,
                                   abdb_handler.adaptor.biodatabase.name,
                                   abdb_handler.adaptor.bioentry.accession,
                                   abdb_handler.adaptor.bioentry.name,
                                   abdb_handler.adaptor.bioentry.description,
                                   abdb_handler.adaptor.biosequence.length,
                                   cache=(cache.ram, CacheTimes.search))
'''

class BioSQLSearch(object):
    '''handles queries in a biosql db.
    use an adapter as input, and optionally a biodatabase name to restrict queries
    for each query  when the .select() method is called
    a list of matching bioentries is returned 
    
    >>> search = BioSQLSearch(BioSQL_handler, biodatabase = 'biodbname')
    >>> search.accession == 'AC123' 
    1
    >>> search.accession.select()
    set([1021])
    >>> search.accession.contains('AC123') #LIKE query, not case sensitive
    4
    >>> search.accession.select() 
    [1021,1023,1024,1025] # bioentry_id list
    >>> search.accession.startswith('AC') #LIKE query, not case sensitive
    10000
    >>> search.accession != 'AC123' 
    9999

    >>> search = BioSQLSearch(BioSQL_handler, biodatabase = 'biodbname')
    >>> search.accession == ['AC123', 'AC124'] # if a list is passed, multiple queries are run with an OR between them
    2
    
    currently available queries:
    - accession
    - id
    - description
    - bioentry_id --> methods: ==, !=, <, <=, >, >=
    - name
    - annotation_type  --> retrieves all the bioentries having a qualifier of the supplied type
    - annotation_value --> retrieves all the bioentries having the supplied qualifier value 
    - annotation_type_value --> retrieves all the bioentries having the supplied couple of key/value
    - feature_type --> retrieves all the bioentries having a feature of the supplied type
    - dbxref
    - seq  --> simple text search on the sequence, not case sensitive
    - keyword  --> retrieves all the bioentries having the supplied qualifier value in a qualifier field starting with 'keyword'
    and with a BioSQL db extended to handle audit trail
    - created_on --> methods: ==, !=
    - created_by --> methods: ==, !=, <, <=, >, >=
    - modified_on --> methods: ==, !=
    - modified_by --> methods: ==, !=, <, <=, >, >=
    
    
    '''
    def __init__(self, handler, biodatabase = None, limits = None):
        self.handler = handler
        self.adaptor = handler.adaptor
        self.limits = limits # limits to be passed in the query
        self.biodb = biodatabase
        
        '''query type mapping'''
        self.quick = BioSQLQueryQuick(self.handler, self.biodb)
        self.fulltext = BioSQLQuery(self.handler, self.biodb)
        self.accession = BioSQLQueryAccession(self.handler, self.biodb)
        self.id = BioSQLQueryID(self.handler, self.biodb)
        self.description = BioSQLQueryDescription(self.handler, self.biodb)
        self.bioentry_id = BioSQLQueryBioentryID(self.handler, self.biodb)
        self.name = BioSQLQueryName(self.handler, self.biodb)
        self.annotation_type = BioSQLQueryAnnotationType(self.handler, self.biodb)
        self.annotation_value = BioSQLQueryAnnotationValue(self.handler, self.biodb)
        self.annotation_type_value = BioSQLQueryAnnotationTypeValue(self.handler, self.biodb)
        self.feature_type = BioSQLQueryFeatureType(self.handler, self.biodb)
        self.dbxref = BioSQLQueryDBXref(self.handler, self.biodb)
        self.seq = BioSQLQuerySeq(self.handler, self.biodb)
        self.keyword = BioSQLQueryKeyword(self.handler, self.biodb)
        if self.handler.time_stamps:
            self.created_on = BioSQLQueryCreatedOn(self.handler, self.biodb)
            self.created_by = BioSQLQueryCreatedBy(self.handler, self.biodb)
            self.modified_on = BioSQLQueryModifiedOn(self.handler, self.biodb)
            self.modified_by = BioSQLQueryModifiedBy(self.handler, self.biodb)