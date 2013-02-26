#from gluon.sql import *
#from gluon.validators import *

if 0:
    from gluon.sql import Field,IS_IN_DB,DAL
    
    





'''This handles biosql dbs '''
# standard modules
from time import gmtime, strftime
import copy
import warnings

# biopython
from Bio import Alphabet
from Bio.SeqUtils.CheckSum import crc64
from Bio import Entrez
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio import SeqFeature




try:
    current_user = auth.user.id
except:
    current_user = 0


try:
    now = request.now
except:
    from datetime import datetime
    now = datetime.now

class ArchiveBioentry():
    '''this is used to make an archieved copy of a bioentry between the same biosql db or
    in an other biosql db used as archive 
    for the receiver db a biodatabase has to be specified
    
    use compatibility_mode = True for archive handler (receiver_handler) to avoid relationship errors. 
    relationships will be safely archived as qualifiers
    
    
    initialize as:
    >>> archiever = ArchiveBioentry(donor_handler, receiver_handler)
    to archive use:
    >>> archiever.archive(donor_bioentry_id)
    or 
    >>> archiever.archive(donor_bioentry_id, receiver_biodatabase)
    '''
    
    def __init__(self, donor_handler, receiver_handler):
        
        self.donor_handler = donor_handler
        self.receiver_handler = receiver_handler
        
    def archive(self, donor_bioentry_id, receiver_biodatabase = None):
        donor_bioentry_id = int(donor_bioentry_id)
        '''fetch seqrecord data'''
        seqrecord = self.donor_handler._retrieve_seqrecord(donor_bioentry_id)
        '''retrieve biodatabase if needed'''
        if not receiver_biodatabase:
            receiver_biodatabase = self.donor_handler._get_bioentry_biodatabase_name(donor_bioentry_id)
        if not receiver_biodatabase in self.receiver_handler.dbs:
            self.receiver_handler.make_new_db(dbname = receiver_biodatabase)
        self.receiver_handler.set_db(receiver_biodatabase)
        '''increase version '''
        accession, version =self.donor_handler._get_accession_from_seqrecord(seqrecord)
        version = int(version)
        last_archieved_version = self.receiver_handler._get_latest_version(accession = accession) # this will search for the latest version in the given biodatabase
        if version <= last_archieved_version:
            version = last_archieved_version +1
        seqrecord.id = '%s.%i'%(accession, version)
        '''load seqrecord to archieve'''
        archieved_bioentry_id = self.receiver_handler.load_seqrecord(seqrecord)
        '''populate timestamps '''
        if self.donor_handler.time_stamps:
            if self.donor_handler.time_stamps_on_biosql:
                timestamps = self.donor_handler.adaptor.bioentry_timestamp[donor_bioentry_id]
            else:
                timestamps = self.donor_handler.timestampsdb.bioentry_timestamp[donor_bioentry_id]
                
        if timestamps:
            self.receiver_handler._set_bioentry_time_stamp(archieved_bioentry_id, 
                                                           timestamps.created_by, 
                                                           timestamps.created_on, 
                                                           timestamps.modified_by, 
                                                           timestamps.modified_on)

            




class RelationRedundancyError(Exception):
     def __init__(self, object, subject, term_id, relation):
         self.object = object
         self.subject = subject
         self.term_id = term_id
         self.relation = relation
     def __str__(self):
         return 'This relation is already present in the db.\n\tobject: %i, subject: %i, term_id: %i, relation: %s'%(self.object,
                                                                                                                   self.subject,
                                                                                                                   self.term_id,
                                                                                                                   repr(self.relation))

class BioEntryDiff(object):
    '''make a diff between two bioentries 
    maybe it is simplier doing it at seqrecord level...
    check for name, accession, in bioentry table
    check for same sequence.
    check for missing dbxref using set in both directions
    check for annotations:
        check for missing anntoation keys using set in both directions
        for intersection check equalties in seqrecord annotation values
    check for features:
        build a feature equalilty object?
        check for same type, same position?
        check for additiona info?
        
    this function shouldn't be fast, will be rarely used, but should be exaustive'''
    def __init__(self,be1, be2):
        self.be1
        self.be2
        
    def diff(self):
         '''Calculate diff '''
         
         return 'NOT IMPLEMENTED'
   

class SeqFeatureEq(SeqFeature.SeqFeature):
    '''extend  the Biopython SeqFeature object to allow
    equalties (a == b) and inequalties (a != b). this will not affect indenty (a is b)  
    this will not check for sub_feature equalties'''
    def __init__(self, seqfeature):
        self.type = seqfeature.type
        self.strand = seqfeature.strand
        self.location = seqfeature.location
        self.id = seqfeature.id
        self.location_operator = seqfeature.location_operator
        self.qualifiers = seqfeature.qualifiers
        self.ref = seqfeature.ref
        self.ref_db = seqfeature.ref_db
        self.sub_features = seqfeature.sub_features
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.type == other.type) and \
                   (self.strand ==  other.strand)  and \
                   (self.location.start.position ==  other.location.start.position)  and \
                   (self.location.end.position ==  other.location.end.position) and \
                   (self.id ==  other.id)  and \
                   (self.location_operator ==  other.location_operator)  and \
                   (self.qualifiers ==  other.qualifiers)  and \
                   (self.ref ==  other.ref)  and \
                   (self.ref_db ==  other.ref_db) #stop at self.id  equalty for simple  behaviour
        else:
            return False
    def __ne__(self, other):
        return not self.__eq__(other)

     
       
class BaseBioSQL():
    '''Base class for to handlers for BioSQL data

    THIS CLASS AND ALL ITS SUBCLASSES ARE DEPRECATED
    SINCE THEY CANNOT BE USED IN A MULTIUSER ENVIROMENT
    FOR WRITE, JUST FOR READING.

    '''
    def __init__(self, handler,):
        warnings.warn("Call to deprecated class",
                      category=DeprecationWarning,
                      stacklevel=2)
        self.handler = handler
        self.adaptor = handler.adaptor
    
    def get(self,):
        raise NotImplementedError('This method is not available')
    def sync(self,):
        raise NotImplementedError('This method is not available')
    def delete(self,):
        raise NotImplementedError('This method is not available')

    def add(self,):
        raise NotImplementedError('This method is not available')

 

class BioSQLFeature(BaseBioSQL):
    '''WORKING ON'''
    
    def __init__(self, handler, bioentry_id, rank = 0, feature = None):
        '''handler is an initiated BioSQLHandler  object
        bioentry_id and  key (qualifier term name) are required to identify the annotation
        
        this handles a single feature add/update/delete
        
        if rank is specified an available annotation is searched,
        otherwise a new one will be created 
        
        usage:
        
        to download an existing feature use
        >>> feat = BioSQLFeature(handler, bioentry_id=123341, rank=1) 
        >>> feat.feature 
        SeqFeature object

        to create a new feature for a given  bioentry
        >>> feat = BioSQLFeature(handler, bioentry_id=123341) 
        >>> feat.get(SeqFeature object)
        or directly
        >>> feat = BioSQLFeature(handler, bioentry_id=123341, feature = SeqFeature object) 
        then
        >>> feat.feature 
        SeqFeature object
        
        modify the feat.feature
        if needed and then use 
        
        >>> feat.sync()
        
        for DB persistence 

        to delete a given feature form the db use 
        
        >>> feat.delete()
        

        '''
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.feature = feature
        self.rank = rank
        if rank != 0:
            self._get()
        
    def get(self, feature):
        ''' import or substitute a SeqFeature object to be loaded in the db
        you need to use sync() for db persistence'''
        self.feature = feature
        
    def sync(self):
        '''IN THE CURRENT IMPLEMENTAION:
        if the feature is not present in the db (rank 0) it is inserted in the db.
        if it is already present, the corresponding row in seqfeature table is updated to preserve seqfeature_id,
        and related location and qualifiers are first removed and then loaded again with updated values.
        TODO:
        granular updating only when needed'''

        if (not self.feature) and self.seqfeature_id:#the feature has to be deleted. AVOID THIS? there is an explicit public methods. check for consistency with other adapters
            self.delete()
        elif self.rank == 0:#check if it a new feature
            self._insert()
        else:#update the feature
            self._update()
        '''the feature data are always written in the db.
        to check only for needed update  download the current feature in the db, 
        and compare it with self.feature using SeqFeatureEq class'''
        self.handler._update_bioentry_time_stamp(self.bioentry_id)
        
        return
    
    def _update(self):
        '''do the same as insert except for _load_seqfeature_basic to maintain seqfeature_id'''
        
        self._update_seqfeature_basic(self.seqfeature_id,
                                      self.feature.type, 
                                      self.rank,
                                      self.bioentry_id)
        
        '''remove dbxrefs '''
        self.adaptor(self.adaptor.seqfeature_dbxref.seqfeature_id == self.seqfeature_id).delete()
        '''remove qualifiers '''
        self.adaptor(self.adaptor.seqfeature_qualifier_value.seqfeature_id == self.seqfeature_id).delete()
        '''remove locations '''
        self.adaptor(self.adaptor.location.seqfeature_id == self.seqfeature_id).delete()
        '''Add updated location and qualifiers '''
        self._load_seqfeature_locations(self.feature, self.seqfeature_id)
        self._load_seqfeature_qualifiers(self.feature.qualifiers, self.seqfeature_id)
    
    def _insert(self):
        """
        Load a biopython SeqFeature into the database (PRIVATE).
        """
        
 
    
        '''main load '''
        if self.rank == 0:#new seqfeature
            ranks = [row.rank for row in self.adaptor((self.adaptor.seqfeature.bioentry_id == self.bioentry_id) & \
                                                       (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) 
                                                       ).select(self.adaptor.seqfeature.rank)]
            if ranks:
                self.rank = max(ranks)+1
            else:
                self.rank = 1
        self.seqfeature_id = self._load_seqfeature_basic(self.feature.type, 
                                                    self.rank,
                                                    self.bioentry_id)
        self._load_seqfeature_locations(self.feature, self.seqfeature_id)
        self._load_seqfeature_qualifiers(self.feature.qualifiers, self.seqfeature_id)
    

    def delete(self):
        '''remove all the data related to a given seqfeature_id
        terms and dbxrefs will remain in the respective tables on the db, 
        since they can be used by other entities, but will no longer be linked to
        the deleted seqfeature '''
        
        '''remove dbxrefs '''
        self.adaptor(self.adaptor.seqfeature_dbxref.seqfeature_id == self.seqfeature_id).delete()
        
        '''remove qualifiers '''
        self.adaptor(self.adaptor.seqfeature_qualifier_value.seqfeature_id == self.seqfeature_id).delete()
        
        '''remove locations '''
        self.adaptor(self.adaptor.location.seqfeature_id == self.seqfeature_id).delete()
        
        '''remove seqfeature '''
        self.adaptor(self.adaptor.seqfeature.seqfeature_id == self.seqfeature_id).delete()
        
        self.handler._update_bioentry_time_stamp(self.bioentry_id)

        
        
    
    def _get(self,):

        '''get id and type'''        
        rows = self.adaptor((self.adaptor.seqfeature.bioentry_id == self.bioentry_id) & \
                           (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) & \
                           (self.adaptor.seqfeature.rank == self.rank)
                           ).select(self.adaptor.seqfeature.seqfeature_id,self.adaptor.term.name)
                           
        if rows:              
            self.seqfeature_id = rows[0].seqfeature.seqfeature_id
            self.type = rows[0].term.name
        else:
            raise IndexError('No feature available at rank %i for bioentry %i'%(self.rank, self.bioentry_id))
        
        '''get qualifiers'''
        qvs=((qualifier.term.name,qualifier.seqfeature_qualifier_value.value,) \
             for qualifier in self.adaptor((self.adaptor.seqfeature_qualifier_value.seqfeature_id == self.seqfeature_id) & \
                      (self.adaptor.seqfeature_qualifier_value.term_id == self.adaptor.term.term_id)\
                      ).select(self.adaptor.term.name,self.adaptor.seqfeature_qualifier_value.value,orderby=self.adaptor.seqfeature_qualifier_value.rank))

        self.qualifiers = {}
        for qv_name, qv_value in qvs:
            self.qualifiers.setdefault(qv_name, []).append(qv_value)
        # Get db_xrefs [special case of qualifiers]
        qvs = ((row.dbname,row.accession) \
             for row in self.adaptor((self.adaptor.seqfeature_dbxref.seqfeature_id == self.seqfeature_id) &  
                                     (self.adaptor.seqfeature_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
                                    ).select(self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,orderby=self.adaptor.seqfeature_dbxref.rank))
        for qv_name, qv_value in qvs:
            value = "%s:%s" % (qv_name, qv_value)
            self.qualifiers.setdefault("db_xref", []).append(value)




        ''' Get locations '''
        results= ((row.location_id,
                   row.start_pos,
                   row.end_pos,
                   row.strand) for row in self.adaptor(self.adaptor.location.seqfeature_id == self.seqfeature_id).select(self.adaptor.location.location_id,
                                                                                                                         self.adaptor.location.start_pos,
                                                                                                                         self.adaptor.location.end_pos,
                                                                                                                         self.adaptor.location.strand,
                                                                                                                         orderby=self.adaptor.location.rank))
        self.locations = []
        # convert to Python standard form
        # Convert strand = 0 to strand = None
        # re: comment in Loader.py:
        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        for location_id, start, end, strand in results:
            if start:
                start -= 1
            if strand == 0:
                strand = None
            if strand not in (+1, -1, None):
                raise ValueError("Invalid strand %s found in database for " \
                                 "seqfeature_id %s" % (strand, seqfeature_id))
            if end < start:
                import warnings
                warnings.warn("Inverted location start/end (%i and %i) for " \
                              "seqfeature_id %s" % (start, end, seqfeature_id))
            self.locations.append( (location_id, start, end, strand) )
            
            
        ''' Get possible remote reference information'''
        remote_results = ((row.location.location_id,
                           row.dbxref.dbname,
                           row.dbxref.accession,
                           row.dbxref.accession,) for row in self.adaptor((self.adaptor.location.seqfeature_id == self.seqfeature_id) &
                                                                      (self.adaptor.location.dbxref_id == self.adaptor.dbxref.dbxref_id)).select(self.adaptor.location.location_id,
                                                                                                                                         self.adaptor.dbxref.dbname,
                                                                                                                                         self.adaptor.dbxref.accession,
                                                                                                                                         self.adaptor.dbxref.accession,))
        self.ref_lookup = {}
        for location_id, dbname, accession, version in remote_results:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            # subfeature remote location db_ref are stored as a empty string when
            # not present
            if dbname == "":
                dbname = None
            self.ref_lookup[location_id] = (dbname, v)
            
        self.feature = SeqFeature.SeqFeature(type =self.type)
        self.feature._seqfeature_id = self.seqfeature_id #Store the key as a private property
        self.feature.qualifiers = self.qualifiers
        if len(self.locations) == 0:
            pass
        elif len(self.locations) == 1:
            location_id, start, end, strand = self.locations[0]
            #See Bug 2677, we currently don't record the location_operator
            #For consistency with older versions Biopython, default to "".
            self.feature.location_operator = \
                self.handler._retrieve_location_qualifier_value(location_id)
            dbname, version = self.ref_lookup.get(location_id, (None, None))
            self.feature.location = SeqFeature.FeatureLocation(start, end)
            self.feature.strand = strand
            self.feature.ref_db = dbname
            self.feature.ref = version
        else:
            assert feature.sub_features == []
            for location in locations:
                location_id, start, end, strand = location
                dbname, version = lookup.get(location_id, (None, None))
                subfeature = SeqFeature.SeqFeature()
                subfeature.type = seqfeature_type
                subfeature.location_operator = \
                    self.handler._retrieve_location_qualifier_value(location_id)
                #TODO - See Bug 2677 - we don't yet record location_operator,
                #so for consistency with older versions of Biopython default
                #to assuming its a join.
                if not subfeature.location_operator:
                    subfeature.location_operator="join"
                subfeature.location = SeqFeature.FeatureLocation(start, end)
                subfeature.strand = strand
                subfeature.ref_db = dbname
                subfeature.ref = version
                self.feature.sub_features.append(subfeature)
            # Assuming that the feature loc.op is the same as the sub_feature
            # loc.op:
            self.feature.location_operator = \
                feature.sub_features[0].location_operator
            # Locations are in order, but because of remote locations for
            # sub-features they are not necessarily in numerical order:
            start = self.locations[0][1]
            end = self.locations[-1][2]
            self.feature.location = SeqFeature.FeatureLocation(start, end)
            # To get the parent strand (as done when parsing GenBank files),
            # need to consider evil mixed strand examples like this,
            # join(complement(69611..69724),139856..140087,140625..140650)
            strands = set(sf.strand for sf in self.feature.sub_features)
            if len(strands)==1:
                self.feature.strand = self.feature.sub_features[0].strand
            else:
                self.feature.strand = None # i.e. mixed strands
    
        return 

    def _update_seqfeature_basic(self, seqfeature_id, feature_type, feature_rank, bioentry_id):
        """Load the first tables of a seqfeature and returns the id (PRIVATE).

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """
        ontology_id = self.handler._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self.handler._get_term_id(feature_type,
                                              ontology_id = ontology_id)
        # XXX source is always EMBL/GenBank/SwissProt here; it should depend on
        # the record (how?)
        source_cat_id = self.handler._get_ontology_id('SeqFeature Sources')
        source_term_id = self.handler._get_term_id('EMBL/GenBank/SwissProt',
                                      ontology_id = source_cat_id)
        
        self.adaptor.seqfeature[seqfeature_id] = dict (bioentry_id = bioentry_id, 
                                                       type_term_id = seqfeature_key_id, 
                                                       source_term_id = source_term_id, 
                                                       rank = feature_rank)

    def _load_seqfeature_basic(self, feature_type, feature_rank, bioentry_id):
        """Load the first tables of a seqfeature and returns the id (PRIVATE).

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """
        ontology_id = self.handler._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self.handler._get_term_id(feature_type,
                                              ontology_id = ontology_id)
        # XXX source is always EMBL/GenBank/SwissProt here; it should depend on
        # the record (how?)
        source_cat_id = self.handler._get_ontology_id('SeqFeature Sources')
        source_term_id = self.handler._get_term_id('EMBL/GenBank/SwissProt',
                                      ontology_id = source_cat_id)
        
        seqfeature_id = self.adaptor.seqfeature.insert(bioentry_id = bioentry_id, 
                                                       type_term_id = seqfeature_key_id, 
                                                       source_term_id = source_term_id, 
                                                       rank = feature_rank)
        

        return seqfeature_id

    def _load_seqfeature_locations(self, feature, seqfeature_id):
        """Load all of the locations for a SeqFeature into tables (PRIVATE).

        This adds the locations related to the SeqFeature into the
        seqfeature_location table. Fuzzies are not handled right now.
        For a simple location, ie (1..2), we have a single table row
        with seq_start = 1, seq_end = 2, location_rank = 1.

        For split locations, ie (1..2, 3..4, 5..6) we would have three
        row tables with:
            start = 1, end = 2, rank = 1
            start = 3, end = 4, rank = 2
            start = 5, end = 6, rank = 3
        """
        # TODO - Record an ontology for the locations (using location.term_id)
        # which for now as in BioPerl we leave defaulting to NULL.
        if feature.location_operator and feature.location_operator != "join":
            # e.g. order locations... we don't record "order" so it
            # will become a "join" on reloading. What does BioPerl do?
            import warnings
            warnings.warn("%s location operators are not fully supported" \
                          % feature.location_operator)
        
        # two cases, a simple location or a split location
        if not feature.sub_features:    # simple location
            self._insert_seqfeature_location(feature, 1, seqfeature_id)
        else: # split location
            for rank, cur_feature in enumerate(feature.sub_features):
                self._insert_seqfeature_location(cur_feature,
                                                 rank + 1,
                                                 seqfeature_id)

    def _insert_seqfeature_location(self, feature, rank, seqfeature_id):
        """Add a location of a SeqFeature to the seqfeature_location table (PRIVATE).

        TODO - Add location_operators to location_qualifier_value.
        """
        # convert biopython locations to the 1-based location system
        # used in bioSQL
        # XXX This could also handle fuzzies
        start = feature.location.nofuzzy_start + 1
        end = feature.location.nofuzzy_end

        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        strand = feature.strand or 0

        # TODO - Record an ontology term for the location (location.term_id)
        # which for now like BioPerl we'll leave as NULL.
        # This might allow us to record "between" positions properly, but I
        # doesn't really see how it could work for before/after fuzzy positions
        loc_term_id = None

        if feature.ref:
            # sub_feature remote locations when they are in the same db as the current
            # record do not have a value for ref_db, which the SeqFeature object
            # stores as None. BioSQL schema requires a varchar and is not NULL 
            dbxref_id = self.handler._get_dbxref_id(feature.ref_db or "", feature.ref)
        else:
            dbxref_id = None

        oid = self.adaptor.location.insert(seqfeature_id = seqfeature_id, 
                                           dbxref_id = dbxref_id, 
                                           term_id = loc_term_id,
                                           start_pos = start, 
                                           end_pos = end, 
                                           strand = strand, 
                                           rank = rank)
        

        """
        # See Bug 2677
        # TODO - Record the location_operator (e.g. "join" or "order")
        # using the location_qualifier_value table (which we and BioPerl
        # have historically left empty).
        # Note this will need an ontology term for the location qualifer
        # (location_qualifier_value.term_id) for which oddly the schema
        # does not allow NULL.
        if feature.location_operator:
            #e.g. "join" (common),
            #or "order" (see Tests/GenBank/protein_refseq2.gb)
            location_id = self.adaptor.last_id('location')
            loc_qual_term_id = None # Not allowed in BioSQL v1.0.1
            sql = r"INSERT INTO location_qualifier_value" \
                  r"(location_id, term_id, value)" \
                  r"VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (location_id, loc_qual_term_id,
                                       feature.location_operator))
        """

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature (PRIVATE).

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        tag_ontology_id = self.handler._get_ontology_id('Annotation Tags')
        for qualifier_key in qualifiers.keys():
            # Treat db_xref qualifiers differently to sequence annotation
            # qualifiers by populating the seqfeature_dbxref and dbxref
            # tables.  Other qualifiers go into the seqfeature_qualifier_value
            # and (if new) term tables.
            if qualifier_key != 'db_xref':
                qualifier_key_id = self.handler._get_term_id(qualifier_key,
                                                  ontology_id=tag_ontology_id)
                # now add all of the values to their table
                entries = qualifiers[qualifier_key]
                if not isinstance(entries, list):
                    # Could be a plain string, or an int or a float.
                    # However, we exect a list of strings here.
                    entries = [entries]
                for qual_value_rank in range(len(entries)):
                    qualifier_value = entries[qual_value_rank]
                    oid = self.adaptor.seqfeature_qualifier_value.insert(seqfeature_id = seqfeature_id, 
                                                                          term_id = qualifier_key_id, 
                                                                          rank = qual_value_rank + 1, 
                                                                          value = qualifier_value)

            else:
                # The dbxref_id qualifier/value sets go into the dbxref table
                # as dbname, accession, version tuples, with dbxref.dbxref_id
                # being automatically assigned, and into the seqfeature_dbxref
                # table as seqfeature_id, dbxref_id, and rank tuples
                self._load_seqfeature_dbxref(qualifiers[qualifier_key],
                                             seqfeature_id)


    def _load_seqfeature_dbxref(self, dbxrefs, seqfeature_id):
        """Add database crossreferences of a SeqFeature to the database (PRIVATE).

            o dbxrefs           List, dbxref data from the source file in the
                                format <database>:<accession>

            o seqfeature_id     Int, the identifier for the seqfeature in the
                                seqfeature table

            Insert dbxref qualifier data for a seqfeature into the
            seqfeature_dbxref and, if required, dbxref tables.
            The dbxref_id qualifier/value sets go into the dbxref table
            as dbname, accession, version tuples, with dbxref.dbxref_id
            being automatically assigned, and into the seqfeature_dbxref
            table as seqfeature_id, dbxref_id, and rank tuples
        """
        # NOTE - In older versions of Biopython, we would map the GenBank
        # db_xref "name", for example "GI" to "GeneIndex", and give a warning
        # for any unknown terms.  This was a long term maintainance problem,
        # and differed from BioPerl and BioJava's implementation.  See bug 2405
        for rank, value in enumerate(dbxrefs):
            # Split the DB:accession format string at colons.  We have to
            # account for multiple-line and multiple-accession entries
            try:
                dbxref_data = value.replace(' ','').replace('\n','').split(':')
                db = dbxref_data[0]
                accessions = dbxref_data[1:]
            except:
                raise ValueError("Parsing of db_xref failed: '%s'" % value)
            # Loop over all the grabbed accessions, and attempt to fill the
            # table
            for accession in accessions:
                # Get the dbxref_id value for the dbxref data
                dbxref_id = self.handler._get_dbxref_id(db, accession)
                # Insert the seqfeature_dbxref data
                self._get_seqfeature_dbxref(seqfeature_id, dbxref_id, rank+1)
        
    def _get_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Check for a pre-existing seqfeature_dbxref entry with the passed
            seqfeature_id and dbxref_id.  If one does not exist, insert new
            data

        """
        # Check for an existing record
        result = self.adaptor((self.adaptor.seqfeature_dbxref.seqfeature_id == seqfeature_id) &
                                 (self.adaptor.seqfeature_dbxref.dbxref_id == dbxref_id)).select(self.adaptor.seqfeature_dbxref.seqfeature_id,
                                                                                                self.adaptor.seqfeature_dbxref.dbxref_id)
        # If there was a record, return without executing anything, else create
        # the record and return
        if result:
            return (result[0].seqfeature_id, result[0].dbxref_id)
        return self._add_seqfeature_dbxref(seqfeature_id, dbxref_id, rank)

    def _add_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Insert a seqfeature_dbxref row and return the seqfeature_id and
            dbxref_id
        """
        
        self.adaptor.seqfeature_dbxref.insert(seqfeature_id = seqfeature_id, 
                                              dbxref_id = dbxref_id, 
                                              rank = rank)
        return (seqfeature_id, dbxref_id)


class BioSQLBioentryFeatures(BaseBioSQL):
    ''' TO DO
    handles feature comparisons and avoid duplications
    
        handler is an initiated BioSQLHandler  object
        bioentry_id and  key (qualifier term name) are required to identify the annotation
        
        this handles a single feature add/update/delete
        
        if rank is specified an available annotation is searched,
        otherwise a new one will be created 
        
        usage:
        
        to download an existing feature use
        >>> feat = BioSQLFeature(handler, bioentry_id=123341, rank=1) 
        >>> feat.feature 
        SeqFeature object

        to create a new feature for a given  bioentry
        >>> feat = BioSQLFeature(handler, bioentry_id=123341) 
        >>> feat.get(SeqFeature object)
        or directly
        >>> feat = BioSQLFeature(handler, bioentry_id=123341, feature = SeqFeature object) 
        then
        >>> feat.feature 
        SeqFeature object
        
        modify the feat.feature
        if needed and then use 
        
        >>> feat.sync()
        
        for DB persistence 

        to delete a given feature form the db use 
        
        >>> feat.delete()
        '''
        
    def __init__(self, handler, bioentry_id,):
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.features = dict() #rank:biosqlfeature dictionary
        self._get()
    
    def _get(self):
        '''get all vailable features '''
        rows = self.adaptor((self.adaptor.seqfeature.bioentry_id == self.bioentry_id)).select()
        if rows:   
            for row in rows:
                if row.rank in self.features:
                    raise ValueError('multiple seqfeatures present with the same rank')
                else:
                    self.features[row.rank] = BioSQLFeature(handler = self.handler, bioentry_id = self.bioentry_id, rank= row.rank)


    def get_ordered_seqfeatures(self):
        ordered = []
        for key in sorted(self.feature):
            ordered.append(self.feature[key].feature)
        return ordered


class BioSQLBioentryDBXrefs(BaseBioSQL):
    
    def __init__(self, handler, bioentry_id, ):
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.dbxrefs = []
        self._get()
        
    def _get(self,):
        dbxrefs=((dbxref.dbname, dbxref.accession, dbxref.version) for dbxref in self.adaptor((self.adaptor.bioentry.bioentry_id == self.bioentry_id) & \
                             (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                             (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id)).select( 
                                                        self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,self.adaptor.dbxref.version, orderby=self.adaptor.bioentry_dbxref.rank))
        for dbname, accession, version in dbxrefs:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            self.dbxrefs.append("%s:%s" % (dbname, v))
        
    
class BioSQLSeq(BaseBioSQL):
    def __init__(self, handler, bioentry_id):
        '''handler is an initiated BioSQLHandler  object
        bioentry_id is required to identify the bioentry sequence. only one 
        sequence per bioentry is allowed by the biosql schema. 
        
        usage:
        >>> seq = BioSQLSeq(handler, bioentry_id=123341, )
        >>> seq.seq 
        Seq('PRLNMLKAWHYTCVNGH', ProteinAlphabet)

        >>> seq.sync()
        
        if the sequence does not exist, a new one is created
        please note that by setting 
        
        >>> seq.seq = '' or [] or None
        

        
        >>> seq.delete()
        
        '''
        self.handler = handler
        self.adaptor = handler.adaptor
        self.bioentry_id = bioentry_id
        self.seq = None
        self.length = None
        self.alphabet = None
        self._get()
        
    def _get(self):
        '''modified from biopython, now returns a Seq object and not a DBSeq object '''
        #The database schema ensures there will be only one matching
        #row in the table.
    
        #If an UnknownSeq was recorded, seq will be NULL,
        #but length will be populated.  This means length(seq)
        #will return None.
        
        seqs = [(row.alphabet,row.length,len(row.seq),row.seq) for row in \
               self.adaptor(self.adaptor.biosequence.bioentry_id == self.bioentry_id \
               ).select(self.adaptor.biosequence.bioentry_id,
                        self.adaptor.biosequence.alphabet,
                        self.adaptor.biosequence.length,
                        self.adaptor.biosequence.seq,)]
        
        if not seqs : return
        assert len(seqs) == 1        
        moltype, given_length, length, seq = seqs[0]
    
        try:
            length = int(length)
            given_length = int(length)
            assert length == given_length
            have_seq = True
        except TypeError:
            assert length is None
            assert len(seqs) == 1
            assert seq is None or seq==""
            length = int(given_length)
            have_seq = False
            #del seq
        del given_length
            
        moltype = moltype.lower() #might be upper case in database
        #We have no way of knowing if these sequences will use IUPAC
        #alphabets, and we certainly can't assume they are unambiguous!
        if moltype == "dna":
            alphabet = Alphabet.generic_dna
        elif moltype == "rna":
            alphabet = Alphabet.generic_rna
        elif moltype == "protein":
            alphabet = Alphabet.generic_protein
        elif moltype == "unknown":
            #This is used in BioSQL/Loader.py and would happen
            #for any generic or nucleotide alphabets.
            alphabet = Alphabet.single_letter_alphabet
        else:
            raise AssertionError("Unknown moltype: %s" % moltype)
    
        if have_seq:
            self.seq = Seq(seq, alphabet,)
            self.length = length
            self.alphabet = moltype
        else:
            self.seq = UnknownSeq(length, alphabet) 
            self.length = length
            self.alphabet = moltype
    
    
    


    def _insert(self, value):
        '''PRIVATE
        load one annotation per time '''
                
        pass
        
        
    def get(self, annotation):
        pass
    
    def sync(self, clean = False):
        '''todo
        
        '''
        
        pass
        
    def delete(self):
        '''remove all the bioentry_qualifier_table rows with the specified key, 
        thus removing the annotation from the BioSQL db'''
        
        pass
  

class BioSQLRelation(BaseBioSQL):
    '''TO DO'''
    pass  

class BioSQLComment(BaseBioSQL):
    '''TO DO. needed?'''
    pass  
     
class BioSQLReference(BaseBioSQL):

    def __init__(self, handler, bioentry_id, ):
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.references = []
        self._get()
        
    def _get(self,):
        '''web2py has a bug with left joins putting left joins before joins in the generated SQL. this is allowed (even if incorrect) by any db-backend except postgres. 
        making a double query and merging the dbname when possible as a workaround here here
        single query working code should be:
        
        refs=self.adaptor((self.adaptor.bioentry_reference.bioentry_id==bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id==self.adaptor.reference.reference_id)
                                              )._select(orderby=self.adaptor.bioentry_reference.rank, left=self.adaptor.reference.on(self.adaptor.reference.dbxref_id==self.adaptor.dbxref.dbxref_id))
        '''

        refs=((row.bioentry_reference.start_pos,
               row.bioentry_reference.end_pos,
               row.reference.location,
               row.reference.title,
               row.reference.authors,
               row.reference.dbxref_id,) for row in self.adaptor((self.adaptor.bioentry_reference.bioentry_id == self.bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id == self.adaptor.reference.reference_id) 
                                              ).select(self.adaptor.reference.reference_id,
                                                       self.adaptor.reference.location,
                                                       self.adaptor.reference.title,
                                                       self.adaptor.reference.authors,
                                                       self.adaptor.reference.dbxref_id,
                                                       self.adaptor.bioentry_reference.start_pos,
                                                       self.adaptor.bioentry_reference.end_pos,
                                                       orderby=self.adaptor.bioentry_reference.rank,))
        refs_dbxrefs=dict(((row.reference.dbxref_id,dict(dbname=row.dbxref.dbname,accession=row.dbxref.accession)) for row in self.adaptor((self.adaptor.bioentry_reference.bioentry_id == self.bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id==self.adaptor.reference.reference_id) &
                                              (self.adaptor.reference.dbxref_id==self.adaptor.dbxref.dbxref_id)
                                              ).select(self.adaptor.reference.dbxref_id,
                                                       self.adaptor.dbxref.dbname,
                                                       self.adaptor.dbxref.accession,
                                                       orderby=self.adaptor.bioentry_reference.rank,)))
        
        for start, end, location, title, authors, dbxref_id in refs:
            if dbxref_id in refs_dbxrefs:
                dbname = refs_dbxrefs[dbxref_id]['dbname']
                accession = refs_dbxrefs[dbxref_id]['accession']
            else:
                dbname,accession=None,None
            reference = SeqFeature.Reference()
            #If the start/end are missing, reference.location is an empty list
            if (start is not None) or (end is not None):
                if start is not None: start -= 1 #python counting
                reference.location = [SeqFeature.FeatureLocation(start, end)]
            #Don't replace the default "" with None.
            if authors : reference.authors = authors
            if title : reference.title = title
            reference.journal = location
            if dbname == 'PUBMED':
                reference.pubmed_id = accession
            elif dbname == 'MEDLINE':
                reference.medline_id = accession
            self.references.append(reference)

    def _insert(self, reference):
    #def _load_reference(self, reference, rank, bioentry_id):
        """Record a SeqRecord's annotated references in the database (PRIVATE).
        record - a SeqRecord object with annotated references
        bioentry_id - corresponding database identifier
        """
        
        ranks = self.adaptor(self.adaptor.bioentry_reference.bioentry_id ==self.bioentry_id).select(self.adaptor.bioentry_reference.rank)
        
        if ranks:
            rank = max([row.rank for row in ranks]) +1
        else:
            rank = 1

        refs = None
        if reference.medline_id:
            refs = self.adaptor((self.adaptor.reference.dbxref_id == self.adaptor.dbxref.dbxref_id) &
                                (self.adaptor.dbxref.dbname == 'MEDLINE') &
                                (self.adaptor.dbxref.accession == reference.medline_id)).select(self.adaptor.reference.reference_id)

        if not refs and reference.pubmed_id:
            refs = self.adaptor((self.adaptor.reference.dbxref_id == self.adaptor.dbxref.dbxref_id) &
                                (self.adaptor.dbxref.dbname == 'PUBMED') &
                                (self.adaptor.dbxref.accession == reference.pubmed_id)).select(self.adaptor.reference.reference_id)
        if not refs:
            s = []
            for f in reference.authors, reference.title, reference.journal:
                s.append(f or "<undef>")
            crc = crc64("".join(s))
            refs = self.adaptor(self.adaptor.reference.crc == crc).select(self.adaptor.reference.reference_id)
        if not refs:
            if reference.medline_id:
                dbxref_id = self.handler._add_dbxref("MEDLINE",
                                             reference.medline_id, 0)
            elif reference.pubmed_id:
                dbxref_id = self.handler._add_dbxref("PUBMED",
                                             reference.pubmed_id, 0)
            else:
                dbxref_id = None
            authors = reference.authors or None
            title =  reference.title or None
            #The location/journal field cannot be Null, so default
            #to an empty string rather than None:
            journal = reference.journal or ""
            reference_id = self.adaptor.reference.insert(dbxref_id = dbxref_id, 
                                                         location = journal,
                                                         title = title, 
                                                         authors = authors, 
                                                         crc = crc)
        else:
            reference_id = refs[0].reference_id

        if reference.location:
            start = 1 + int(str(reference.location[0].start))
            end = int(str(reference.location[0].end))
        else:
            start = None
            end = None
            
        self.adaptor.bioentry_reference.insert(bioentry_id = self.bioentry_id, 
                                               reference_id = reference_id,
                                               start_pos = start, 
                                               end_pos = end, 
                                               rank = rank + 1)

class BioSQLTaxon(BaseBioSQL):
    '''TO DO. needed?'''
    pass  

     
class BioSQLQualifier(BaseBioSQL):
    '''only works for values stored in the bioentry qualifier table  '''
    def __init__(self, handler, bioentry_id, key):
        '''handler is an initiated BioSQLHandler  object
        bioentry_id and  key (qualifier term name) are required to identify the annotation
        annotation values MUST be a list
        
        usage:
        >>> ann = BioSQLQualifier(handler, bioentry_id=123341, key='ann key')
        >>> ann.key 
        'ann key'
        >>> ann.value 
        ['ann value 1', 'ann value 2']
        >>> ann.value = ['ann value 3']
        >>> ann.value.append('ann value 4')
        >>> ann.value 
        ['ann value 3', 'ann value 4']
        
        >>> ann.sync()
        
        if the annotation key do not exist, a new one is created
        please note that by setting 
        
        >>> ann.value = '' or [] or None
        
        will not delete the annotation, but leave the key term in the bioentry_qualifier_value table
        to delete an annotation key and all its values use
        
        >>> ann.delete()
        
        warning, if multiple users are using the BioSQL db, this class
        will not map the changes done to the DB since the object
        initialization. in the worst case you can have duplicated 
        annotations. in that case a new init and sync should fix the issue.
        a double query at the sync time is not performed for performance.
        
        '''
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.key = key
        self.value = []
        self.rank2value = {}
        self._get()
        
    def _get(self):
        qvs=[qualifier for qualifier in self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == self.bioentry_id) & \
                                                      (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) & \
                                                      (self.adaptor.term.name == self.key) \
                                                      ).select(self.adaptor.bioentry_qualifier_value.value, 
                                                               self.adaptor.bioentry_qualifier_value.rank, 
                                                               orderby=self.adaptor.bioentry_qualifier_value.rank)]
        if qvs:
            for qv in qvs:
                self.value.append(qv.value)
                self.rank2value[qv.rank] = qv.value

        '''save original data to avoid non needed overwrites'''
        self._original_values = copy.deepcopy(self.value)
        self._original_rank2value = copy.deepcopy(self.rank2value)

        if qvs:
            self._lower_rank = max(self.rank2value) + 1
        else:
            self._lower_rank = 0

    def _insert(self, value, rank = None):
        '''PRIVATE
        load one annotation per time'''
        if rank == None:
            self.adaptor.bioentry_qualifier_value.insert(bioentry_id = self.bioentry_id,
                                                         term_id =self._term_id,
                                                         value = str(value),
                                                         rank = self._lower_rank )
            self._lower_rank  += 1
        else:
            self.adaptor.bioentry_qualifier_value.insert(bioentry_id = self.bioentry_id,
                                                         term_id =self._term_id,
                                                         value = str(value),
                                                         rank = rank )
    def update_single_value(self, value, rank ):
        '''updates a single qualifier value given it's rank and new value'''
        key_term = self.adaptor(self.adaptor.term.name == self.key).select(self.adaptor.term.term_id).first()

        self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == self.bioentry_id) & \
                     (self.adaptor.bioentry_qualifier_value.term_id == key_term.term_id) & \
                     (self.adaptor.bioentry_qualifier_value.rank == rank)).update(value = str(value))
        self.handler._update_bioentry_time_stamp(self.bioentry_id)
                                                         
        
    def get(self, annotation):
        if isinstance(annotation, str):
            self.value = [annotation]
        elif isinstance(annotation, list) or isinstance(annotation, tuple):
            self.value = annotation
        else:
            raise TypeError('Type "%s" is not supported for upload to database. Only strings, lists and tuples supported'%type(self.value))
    
    def sync(self, clean = False):
        '''
        qualifiers rank will be normalized ar each sync
        
        todo
        since there is no id in bioentry_qualifier_value, each missing entry must be deleted with a query, is it feasible???
        
        #============== Web2py DAL ===================#
        count, delete, update
        
        
        You can count records in a set:
        
        count
        >>> print db(db.person.id > 0).count()

        
        
        You can delete records in a set:
        
        delete
        >>> db(db.person.id > 3).delete()
        
        delete works with single table sets, no joins allowed!
        
        
        
        And you can update all records in a set by passing named arguments corresponding to the fields that need to be updated:
        
        update
        >>> db(db.person.id > 3).update(name='Ken')
        #=================================#
        
        additional values can be added easily
        
        if clean is called, all previously present values will be deleted and all the values will be added from scratch.
        is this useful????
        
        '''
        
        if isinstance(self.value, str):
            self.value = [self.value]
        elif isinstance(self.value, list) or isinstance(self.value, tuple):
            pass
        else:
            raise TypeError('Type "%s" is not supported for upload to database. Only strings, lists and tuples supported'%type(self.value))
        
        '''remove previous values no longer present'''
        self._tag_ontology_id = self.handler._get_ontology_id('Annotation Tags')
        self._term_id = self.handler._get_term_id(self.key, ontology_id=self._tag_ontology_id)
        
        if clean:
            values_to_remove = self._original_values
        else:
            values_to_remove = set(self._original_values) - set(self.value)
        for entry in values_to_remove:
            self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == self.bioentry_id) & \
                          (self.adaptor.bioentry_qualifier_value.term_id == self._term_id) & \
                          (self.adaptor.bioentry_qualifier_value.value == entry) \
                          ).delete() #delete works with single table sets, no joins allowed
                              
        '''add additional values if needed.
        not using sets to mantain list order'''
        for entry in self.value:
            if entry not in self._original_values:
                if isinstance(entry, str) or isinstance(entry, int):
                    self._insert(entry)
                else:
                    raise TypeError('Type "%s" value is not supported for upload to database'%type(entry))
                
        self.handler._update_bioentry_time_stamp(self.bioentry_id)
        return 

    def delete(self):
        '''remove all the bioentry_qualifier_table rows with the specified key, 
        thus removing the annotation from the BioSQL db'''
        
        tag_ontology_id = self.handler._get_ontology_id('Annotation Tags')
        term_id = self.handler._get_term_id(self.key, ontology_id=tag_ontology_id)
        
        self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == self.bioentry_id) & \
                     (self.adaptor.bioentry_qualifier_value.term_id == term_id)).delete() 
                     
        self.handler._update_bioentry_time_stamp(self.bioentry_id)
                
class BioSQLBioentryQualifiers(BaseBioSQL):

    def __init__(self, handler, bioentry_id, ):
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.taxon_id = self.adaptor.bioentry[bioentry_id].taxon_id
        self.simple_qualifiers = {}
        self.qualifiers = []#not needed? each qualifiers can be referred as a single entry with its name when needed without making a query for each qualifier here
        self._get()
        
    def _get(self):
        annotations = {}
        qvs=((qualifier.term.name,qualifier.bioentry_qualifier_value.value,) \
                     for qualifier in self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == self.bioentry_id) & \
                              (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id)\
                              ).select(self.adaptor.term.name,self.adaptor.bioentry_qualifier_value.value,orderby=self.adaptor.bioentry_qualifier_value.rank))
        added_names = []
        for name, value in qvs:
            annotations.setdefault(name, []).append(value)
            if not name in added_names:
                self.qualifiers.append(BioSQLQualifier(self.handler, bioentry_id=self.bioentry_id, key=name))
                added_names.append(name)

        for key, val in annotations.items():
            if isinstance(val, list):
                val = [self.handler._make_unicode_into_string(x) for x in val]
            elif isinstance(val, unicode):
                val = str(val)
            self.simple_qualifiers[key] = val

               
                             
class BioSQLBioentryAnnotations(BaseBioSQL):
    ''' this is used to handle all the BioSQL inforamtion
    that will be reported as a seqrecord annotation of a given
    bioentry_id '''
         
    def __init__(self, handler, bioentry_id, taxon_id = None):
        
        BaseBioSQL.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.taxon_id = self.adaptor.bioentry[bioentry_id].taxon_id
        self.value = self._get()
        
    def _get(self):
        annotations = {}
        annotations.update(self.handler._retrieve_qualifier_value(self.bioentry_id))
        annotations.update(self.handler._retrieve_reference(self.bioentry_id))
        annotations.update(self.handler._retrieve_taxon(self.bioentry_id, self.taxon_id))
        annotations.update(self.handler._retrieve_comment(self.bioentry_id))
        annotations.update(self.handler._retrieve_bioentry_relationship(self.bioentry_id))
        str_anns = {}
        
        for key, val in annotations.items():
            if isinstance(val, list):
                val = [self.handler._make_unicode_into_string(x) for x in val]
            elif isinstance(val, unicode):
                val = str(val)
            str_anns[key] = val
        return str_anns
    
class BioSQLBioentry(BaseBioSQL):
    
    '''
    USAGE:
    
    loads data from db
    >>> sr = BioSQLBioentry(bioentry_id = 1234)
    
    add a seqrecord on db:
    >>> new_sr = SeqRecord('MKVSLAMSPRL')
    >>> sr = BioSQLBioentry(seqrecord = new_sr)
    >>> sr.sync()
    
    modify seqrecord on db:
    >>> sr_mod = SeqRecord('MKVSLAMSPRL')
    >>> sr = BioSQLBioentry(bioentry_id = 1234, seqrecord = sr_mod)
    >>> sr.sync()
    
    or
    >>> sr_mod = SeqRecord('MKVSLAMSPRL')
    >>> sr = BioSQLBioentry(bioentry_id = 1234)
    >>> sr.get(sr_mod)
    >>> sr.sync()
    
    delete seqrecord from db:
    >>> BioSQLBioentry(bioentry_id = 1234).delete()
    
    '''
    
    
    def __init__(self, handler, bioentry_id = None, seqrecord = None):
        self.bioentry_id = bioentry_id
        if bioentry_id and not seqrecord: #load data from db
            self._get(bioentry_id)
        elif seqrecord: #load data from supplied seqrecord
            self.bioentry_id = 0
            self.get(seqrecord)

    
    def _update(self):
        pass
    def _insert(self):
        pass
    def _get(self, bioentry_id):
        
        rows = self.adaptor(self.adaptor.bioentry_id == bioentry_id).select()
        if rows:
            self.bioentry_id == rows[0].bioentry_id
            self.biodatabase_id == rows[0].biodatabase_id
            self.biodatabase = self.adaptor.biodatabase[self.biodatabase_id].name
            self.taxon_id == rows[0].taxon_id
            self.name == rows[0].name
            self.accession == rows[0].accession
            self.identifier == rows[0].identifier
            self.division == rows[0].division
            self.description == rows[0].description
            self.version == rows[0].version
            self.seqrecord = self._retrieve_seqrecord(bioentry_id)
        else:
            raise IndexError('bioentry_id %s not available in db'%bioentry_id)

        pass

    
    def _retrieve_seqrecord(self,bioentry_id,):
        """create a SeqRecord object from BioSQL DB and populate other proprierties in self."""

        if self.version and self.version != "0":
            self.seqrec_id = "%s.%s" % (self.accession, self.version)
        else:
            self.seqrec_id = self.accession

        self.seq = BioSQLBioentrySeq(bioentry_id)
        try:
            length = len(self.seqrec_seq)
        except:#Could be no sequence in the database!
            length = 0
        self.per_letter_annotations = _RestrictedDict(length=length) # not handled in biopython
        self.dbxrefs = BioSQLBioentryDBXrefs(bioentry_id)
        self.features = BioSQLBioentryFeatures(bioentry_id)
        self.annotations = BioSQLBioentryAnnotations(bioentry_id, taxon_id)

        self.seqrecord=SeqRecord(self.seq.export(), id = self.seqrec_id, name = self.name, description = self.description)
        self.seqrecord.dbxrefs = self.dbxrefs.export()
        self.seqrecord.features = self.features.export()
        self.seqrecord.annotations = self.annotations.export()

    def export(self):
        return self.seqrecord
    
    

    
    
class BioSQLHandler(object):
    '''handles connection for read and (partial) writes to a biosql db
    
    postgres://postgres:bioseq@localhost/humangenome
    
    db = BioSQLHandler(connection_string)
    
    query example:
    
    for row in db.adaptor(db.adaptor.bioentry.id>0).select():
        print row.bioentry_id, row.accession, row.identifier, row.name
    
    add relation between bioentry:
    
    relation_oid = db._add_bioentry_relation(obj_id,sbj_id,'relation type')
    
    differences between creating tables from web2py and loading official sql code:
    - tables with missing unique incremental id are added with an id field. this will not be used 
    for compatibility in the queries or anywhere in the code. web2py needs it to be there.
    - multicolon dababase constrains are not encoded in the db backend, but are used in the web2py software layer.
    web2py cannot add multicolumn constrains
    - a bioentry_timestamps table is added for usage logging. this can be created in a different db for retrocompatibility.
    - ondelete CASCADE is currently not working under postgresql (bug?) and sqlite (triggers have to be used)
    
    '''
    
    def __init__(self, 
                 conn_string,
                 fetch_NCBI_taxonomy=False, 
                 compatibility_mode = False,
                 time_stamps = False,
                 pool_size = 5):
        '''compatibility_mode makes seqrecord upload similar to biopython, thus relationships are uploaded as annotation (useful for archive purpose)
        time_stamps is set to a connection string to enable timestamps'''
        
        self._build_error = False
        self.comp_mode = compatibility_mode #TODO, if True follow the biopython schema to save data, else use improved shema filling
        self.dbtype = conn_string.split(':')[0]
        if time_stamps:
            self.time_stamps = True
            if time_stamps == conn_string or time_stamps == True:
                self.time_stamps_on_biosql = True
            else:
                self.time_stamps_on_biosql = False
        else:
            self.time_stamps = False
            self.time_stamps_on_biosql = False
        
        try:
            self.adaptor= DAL(conn_string, pool_size=pool_size)
            self._build_model()#builds the web2py model
            self.dbs = dict([(row.name, row.biodatabase_id) for row in self.adaptor(self.adaptor.biodatabase.biodatabase_id > 0).select()])
        except:
            self.adaptor= DAL(conn_string, pool_size=pool_size)
            self._build_model(create_tables=True)#create new db from scratch
            self.adaptor= DAL(conn_string, pool_size=pool_size)
            self._build_model()#builds the web2py model with tables set at migrate = False
            try:
                self.dbs = dict([(row.name, row.biodatabase_id) for row in self.adaptor(self.adaptor.biodatabase.biodatabase_id > 0).select()])
            except:
                self._build_error = True
        self.fetch_NCBI_taxonomy = fetch_NCBI_taxonomy
        if self.time_stamps:
            if not self.time_stamps_on_biosql:
                self._build_external_time_stamps_table(time_stamps)

        
    def set_db(self,dbname):
        if dbname in self.dbs:
            self.dbid = self.dbs[dbname]
            return
        else:
            raise Exception('biodatabase %s not present in BioSQL db.'%dbname)
        
    def make_new_db(self,dbname, authority='', description = ''):
         if dbname in self.dbs:
             raise Exception('biodatabase %s is already present in BioSQL db.'%dbname)
         else:
             dbid = self.adaptor.biodatabase.insert(name = dbname,
                                                    authority = authority,
                                                    description = description)
             self.dbs[dbname]=dbid
             return dbid
             
        
    def load_seqrecord(self, record, db = None):
        """Load a Biopython SeqRecord into the database.
        added db value to enable usage like this:
        
        >>> biodb_handler.load_seqrecord(seqrecord, db = 'DB Name')
        1
        
        """
        if not db:
            if 'dbid' not in self.__dict__:
                raise AttributeError('You need to set a biodatabase to load a seqrecord')
        else:
            if 'dbid' not in self.__dict__:
                self.set_db(db)
                original_dbid = self.dbid
            else:
                original_dbid = self.dbid 
                self.set_db(db)
        bioentry_id = self._load_bioentry_table(record)
        self._load_bioentry_date(record, bioentry_id)
        self._load_biosequence(record, bioentry_id)
        self._load_comment(record, bioentry_id)
        self._load_dbxrefs(record, bioentry_id)
        references = record.annotations.get('references', ())
        for reference, rank in zip(references, range(len(references))):
            self._load_reference(reference, rank, bioentry_id)
        self._load_annotations(record, bioentry_id)
        for seq_feature_num in range(len(record.features)):
            seq_feature = record.features[seq_feature_num]
            self._load_seqfeature(seq_feature, seq_feature_num, bioentry_id)
        if self.time_stamps:
            self._load_bioentry_time_stamp(bioentry_id)
        self.dbid = original_dbid
        return bioentry_id #added return instead of None in biopython
        #make a single commit at this point!!!

    def load_seqrecord_as_new_version(self, record):
        """Load a Biopython SeqRecord into the database.
        if a bioentry exixts with the same accession in the same biodb
        it create a newer version. if no other accession exists it 
        will initialize the version to 1.
        """
        accession,version = self._get_accession_from_seqrecord(record)
        record.id = '%s.%i'%(accession,version+1)
        return self.load_seqrecord(record)
    
    

    def _build_external_time_stamps_table(self, time_stamps):
        
        self.timestampsdb = DAL(time_stamps, pool_size=5)
        table_name = 'bioentry_timestamp'
        self.timestampsdb.define_table(table_name,
            Field('bioentry_id','integer', requires = IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name')),
            Field('created_by', default = current_user,),# writable = False),
            Field('created_on','datetime', default = now,),# writable = False) ,
            Field('modified_by', default = current_user, ),# writable = False),
            Field('modified_on','datetime', default = now, ))# writable = False),)

    def _build_model(self, create_tables = False, FAKE_MIGRATE  = False ):
        '''BioSQL db model
        WARNING: addional primarykeys are created by web2py since it requires 'id' fields,
        thus there will also be additional corresponding seqeunces.
        multicolum primary key are not supported, and the uniqueness must be checked
        at application level.'''
                
        def get_sequence_name(tablename):
            
            table_seq = None #use web2py standard sequence
            
            if self.dbtype == 'postgres':#biosql schema uses a non standard sequence name in postgresql
                '''default, if the tables were created by web2py.
                warning this will create different sequence names in postgresql with respect to the official biosql schema'''
                table_seq = '%s_pk_seq'%table_name
                '''check for already available sequences, if the tables were created by injecting sql from the official biosql schema '''
                try:
                    col_defaults = self.adaptor.executesql('''SELECT column_default, column_name  FROM information_schema.columns WHERE table_name = '%s'; '''%tablename, as_dict = True)
                    for col in col_defaults:
                        if col['column_default']:
                            table_seq = col['column_default'].split("'")[1]
                            break
                except: #table do not exists
                    pass
            return table_seq
            
        
        alter_db_list = ['postgres']# list of supported dbs to apply schema changing
        
        table_name = 'biodatabase'
        self.adaptor.define_table(table_name,
            Field('biodatabase_id','id', unique = True),
            Field('name','string', length=128, notnull = True, unique = True),
            Field('authority','string', length=128),
            Field('description','text'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX db_auth ON biodatabase (authority);''')
        

        table_name = 'taxon'
        self.adaptor.define_table(table_name,
            Field('taxon_id','id', unique = True),
            Field('ncbi_taxon_id','integer', unique = True),
            Field('parent_taxon_id','integer'),
            Field('node_rank','string', length=32),
            Field('genetic_code','integer'),
            Field('mito_genetic_code','integer'),
            Field('left_value','integer', unique = True),
            Field('right_value','integer',unique = True),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX taxparent ON taxon (parent_taxon_id);''')
            
        table_name = 'taxon_name'
        self.adaptor.define_table(table_name,
            Field('taxon_id',
                  'reference taxon', 
                  unique = True, 
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'taxon.taxon_id',), 
                  ondelete='CASCADE'),
            Field('name','string',length=255, notnull = True,),
            Field('name_class','string',length=32, notnull = True,),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            primarykey=['taxon_id', 'name','name_class'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE taxon_name DROP COLUMN id;''')
                self.adaptor.executesql('''ALTER TABLE taxon_name ADD UNIQUE ( taxon_id,name,name_class );''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
            self.adaptor.executesql('''CREATE INDEX taxnametaxonid ON taxon_name(taxon_id);''')
            self.adaptor.executesql('''CREATE INDEX taxnamename    ON taxon_name(name);''')
            


            
        table_name = 'ontology'
        self.adaptor.define_table(table_name,
            Field('ontology_id','id',  unique = True),
            Field('name','string', length=32, notnull = True, unique = True),
            Field('definition','text'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        
        table_name = 'term'
        self.adaptor.define_table(table_name,
            Field('term_id','id',  unique = True),
            Field('name','string', length=255, notnull = True, unique = True),
            Field('definition','text'),
            Field('identifier','string', length=40, unique = True),
            Field('is_obsolete','string', length=1),
            Field('ontology_id','reference ontology', 
                                notnull = True, 
                                requires=IS_IN_DB(self.adaptor,'ontology.ontology_id','ontology.name'), 
                                ondelete='CASCADE'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX term_ont ON term (ontology_id);''')
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE term ADD UNIQUE (name, ontology_id, is_obsolete);''')

        table_name = 'term_synonym'
        self.adaptor.define_table(table_name,
            Field('synonym','string',length=255, notnull = True),
            Field('term_id','reference term',requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), ondelete='CASCADE'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['term_id', 'synonym'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE term_synonym DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE term_synonym ADD PRIMARY KEY (term_id, synonym);''')

        table_name = 'dbxref'
        self.adaptor.define_table(table_name,
            Field('dbxref_id','id',  unique = True),
            Field('dbname','string', notnull = True),
            Field('accession','string', length=128, notnull = True),
            Field('version','integer', notnull = True),
            migrate = create_tables,
            #missing multicolumn constraint :'CONSTRAINT dbxref_accession_key UNIQUE (accession, dbname, version)'
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX dbxref_db ON dbxref (dbname);''')
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE dbxref ADD UNIQUE (accession, dbname, version);''')

    
        table_name = 'term_dbxref'
        self.adaptor.define_table(table_name,
            Field('term_id', 'reference term',requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), ondelete='CASCADE'),
            Field('dbxref_id','reference dbxref',requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id','dbxref.accession'), ondelete='CASCADE'),    
            Field('rank','integer',),# default = 0),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['term_id', 'dbxref_id'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE term_dbxref DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE term_dbxref ADD PRIMARY KEY (term_id, dbxref_id);''')
            self.adaptor.executesql('''CREATE INDEX trmdbxref_dbxrefid ON term_dbxref(dbxref_id);''')


        
        table_name = 'term_relationship'
        self.adaptor.define_table(table_name,
            Field('term_relationship_id','id', unique = True),
            Field('subject_term_id','reference term',requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), ondelete='CASCADE'),
            Field('predicate_term_id','reference term',requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), ondelete='CASCADE'),
            Field('object_term_id','reference term',requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), ondelete='CASCADE'),
            Field('ontology_id','reference ontology',requires=IS_IN_DB(self.adaptor,'ontology.ontology_id','ontology.name'), ondelete='CASCADE'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE term_relationship ADD UNIQUE (subject_term_id,predicate_term_id,object_term_id,ontology_id);''')
            self.adaptor.executesql('''CREATE INDEX trmrel_predicateid ON term_relationship (predicate_term_id);''')
            self.adaptor.executesql('''CREATE INDEX trmrel_objectid ON term_relationship (object_term_id);''')
            self.adaptor.executesql('''CREATE INDEX trmrel_ontid ON term_relationship (ontology_id);''')
        
        table_name = 'term_relationship_term'
        self.adaptor.define_table(table_name,
            Field('term_relationship_id','id', 
                  unique = True, 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term', 
                  unique = True, 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='CASCADE'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
            
            
        table_name = 'term_path'
        self.adaptor.define_table(table_name,
            Field('term_path_id','id', unique = True),
            Field('subject_term_id',
                  'reference term', 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='CASCADE'),
            Field('predicate_term_id',
                 'reference term', 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='CASCADE'),
            Field('object_term_id',
                  'reference term', 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='CASCADE'),
            Field('ontology_id',
                  'reference ontology',
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'ontology.ontology_id','ontology.name'), 
                  ondelete='CASCADE'),
            Field('distance','integer', ),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE term_path ADD UNIQUE (subject_term_id, predicate_term_id, object_term_id, ontology_id,distance);''')
            self.adaptor.executesql('''CREATE INDEX trmpath_predicateid ON term_path ( predicate_term_id );''')
            self.adaptor.executesql('''CREATE INDEX trmpath_objectid ON term_path ( object_term_id );''')
            self.adaptor.executesql('''CREATE INDEX trmpath_ontid ON term_path ( ontology_id );''')

        table_name = 'bioentry'
        self.adaptor.define_table(table_name,
            Field('bioentry_id','id', unique = True),
            Field('biodatabase_id','reference biodatabase',  
                  requires = IS_IN_DB(self.adaptor,'biodatabase.biodatabase_id','biodatabase.name'),  
                  ondelete = 'NO ACTION',
                  notnull = True),
            Field('taxon_id','reference taxon',  
                  requires = IS_IN_DB(self.adaptor,'taxon.taxon_id',),  
                  ondelete = 'NO ACTION'),
            Field('name','string', length=40, notnull = True),
            Field('accession','string', length=128, notnull = True, unique = True),
            Field('identifier','string', length=40, unique = True),
            Field('division','string', length=6),
            Field('description','text'),
            Field('version','integer', notnull = True, default = 0),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry ADD UNIQUE (identifier, biodatabase_id);''')
                self.adaptor.executesql('''ALTER TABLE bioentry ADD UNIQUE (accession, biodatabase_id, version);''')
            self.adaptor.executesql('''CREATE INDEX bioentry_name ON bioentry (name);''')
            self.adaptor.executesql('''CREATE INDEX bioentry_acc ON bioentry (accession);''')#added not present in original biosql
            self.adaptor.executesql('''CREATE INDEX bioentry_tax ON bioentry (taxon_id);''')
            self.adaptor.executesql('''CREATE INDEX bioentry_db ON bioentry (biodatabase_id);''')

            
        if self.time_stamps_on_biosql:
            table_name = 'bioentry_timestamp'
            self.adaptor.define_table(table_name,
                Field('bioentry_id','reference bioentry',  unique = True, requires = IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name')),
                Field('created_by', default = current_user,),# writable = False),
                Field('created_on','datetime', default = now,),# writable = False) ,
                Field('modified_by', default = current_user, ),# writable = False),
                Field('modified_on','datetime', default = now, ),# writable = False),
                migrate = create_tables,# you need write access to use timestamps
                sequence_name = get_sequence_name(table_name))

        table_name = 'bioentry_relationship'
        self.adaptor.define_table(table_name,
            Field('bioentry_relationship_id','id',  unique = True, ),
            Field('object_bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('subject_bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('rank','integer', ),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry_relationship ADD UNIQUE (object_bioentry_id,subject_bioentry_id,term_id);''')
            self.adaptor.executesql('''CREATE INDEX bioentryrel_child ON bioentry_relationship (subject_bioentry_id);''')
            self.adaptor.executesql('''CREATE INDEX bioentryrel_trm ON bioentry_relationship (term_id);''')

        table_name = 'bioentry_path'
        self.adaptor.define_table(table_name,
            Field('object_bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('subject_bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('distance','integer', ),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            primarykey=['object_bioentry_id', 'subject_bioentry_id', 'term_id', 'distance'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry_path DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE bioentry_path ADD UNIQUE (object_bioentry_id, subject_bioentry_id, term_id, distance);''')
            self.adaptor.executesql('''CREATE INDEX bioentrypath_child ON bioentry_path (subject_bioentry_id);''')
            self.adaptor.executesql('''CREATE INDEX bioentrypath_trm ON bioentry_path (term_id);''')

 

        table_name = 'biosequence'
        self.adaptor.define_table(table_name,
            Field('bioentry_id','reference bioentry', unique = True, requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), ondelete='CASCADE'), 
            Field('length','integer'),
            Field('alphabet','string', length=10),
            Field('seq','text'),
            Field('version','integer'),#not used in biopython
            ##primarykey=['bioentry_id'],#cannot be used since in web2py it is not possible to make a reference form a keyedtable to a normal table
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['bioentry_id',],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE biosequence DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE biosequence ADD PRIMARY KEY (bioentry_id);''')   
        
        table_name = 'dbxref_qualifier_value'
        self.adaptor.define_table(table_name,
            Field('dbxref_id',
                  'reference dbxref', 
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id',), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('rank','integer', notnull = True,  default = 0),
            Field('value','text'),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['dbxref_id','term_id', 'rank'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE dbxref_qualifier_value DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE dbxref_qualifier_value ADD PRIMARY KEY (dbxref_id, term_id, rank);''')  
            self.adaptor.executesql('''CREATE INDEX dbxrefqual_dbx ON dbxref_qualifier_value (dbxref_id);''')
            self.adaptor.executesql('''CREATE INDEX dbxrefqual_trm ON dbxref_qualifier_value (term_id);''')

        table_name = 'bioentry_dbxref'
        self.adaptor.define_table(table_name,
            Field('bioentry_id',
                 'reference bioentry',
                 notnull = True,
                 requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                 ondelete='CASCADE'),
            Field('dbxref_id',
                  'reference dbxref',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id','dbxref.name'), 
                  ondelete='CASCADE'),
            Field('rank','integer',),# default = 0),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['bioentry_id','dbxref_id',],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry_dbxref DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE bioentry_dbxref ADD PRIMARY KEY (bioentry_id, dbxref_id);''')  
            self.adaptor.executesql('''CREATE INDEX dblink_dbx ON bioentry_dbxref (dbxref_id);''')
            
            
        table_name = 'reference'
        self.adaptor.define_table(table_name,
            Field('reference_id','id', unique = True),
            Field('dbxref_id',
                  'reference dbxref',
                  unique = True,
                  requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id','dbxref.accession'), 
                  ondelete='NO ACTION'),    
            Field('location','text', notnull= True),
            Field('title','text'),
            Field('authors','text'),
            Field('crc','string', length=32, unique = True),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        
        table_name = 'bioentry_reference'
        self.adaptor.define_table(table_name,
            Field('bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('reference_id',
                  'reference reference',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'reference.reference_id','reference.title'), 
                  ondelete='CASCADE'),
            Field('start_pos','integer'),
            Field('end_pos','integer'),
            Field('rank','integer', default = 0, notnull = True),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['bioentry_id','reference_id', 'rank'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry_reference DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE bioentry_reference ADD PRIMARY KEY (bioentry_id, reference_id, rank);''')  
            self.adaptor.executesql('''CREATE INDEX bioentryref_ref ON bioentry_reference (reference_id);''')
            
        table_name = 'comment'  
        self.adaptor.define_table(table_name,
            Field('comment_id','id', unique = True),
            Field('bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),   
            Field('rank','integer', default = 0, notnull = True),
            Field('comment_text','text', notnull = True),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE comment ADD UNIQUE (bioentry_id, rank);''')  

            
        table_name = 'bioentry_qualifier_value'
        self.adaptor.define_table(table_name,
            Field('bioentry_id',
                  'reference bioentry', 
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True, 
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('value','text'),
            Field('rank','integer', notnull = True, default = 0),#cannot use default = 0 in postgresql, web2py replicates the given field, report an issue!!!
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            primarykey=[],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE bioentry_qualifier_value DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE bioentry_qualifier_value ADD UNIQUE (bioentry_id, term_id, rank);''')  
            self.adaptor.executesql('''CREATE INDEX bioentryqual_trm ON bioentry_qualifier_value (term_id);''')



        table_name = 'seqfeature'
        self.adaptor.define_table(table_name,
            Field('seqfeature_id','id',  unique = True),
            Field('bioentry_id',
                  'reference bioentry',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'bioentry.bioentry_id','bioentry.name'), 
                  ondelete='CASCADE'),    
            Field('type_term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('source_term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('display_name','string', length=64),
            Field('rank','integer', notnull = True, default = 0),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX seqfeature_trm  ON seqfeature(type_term_id);''')
            self.adaptor.executesql('''CREATE INDEX seqfeature_fsrc ON seqfeature(source_term_id);''')
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE seqfeature ADD UNIQUE (bioentry_id,type_term_id,source_term_id,rank);''')  
            
            
        table_name = 'seqfeature_relationship'
        self.adaptor.define_table(table_name,
            Field('seqfeature_relationship_id','id',  unique = True),
            Field('object_seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.name'), 
                  ondelete='CASCADE'),
            Field('subject_seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.seqfeature_id'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('rank','integer',),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            self.adaptor.executesql('''CREATE INDEX seqfeaturerel_trm   ON seqfeature_relationship(term_id);''')
            self.adaptor.executesql('''CREATE INDEX seqfeaturerel_child ON seqfeature_relationship(subject_seqfeature_id);''')
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE seqfeature_relationship ADD UNIQUE (object_seqfeature_id,subject_seqfeature_id,term_id);''')  
            
        table_name = 'seqfeature_path'
        self.adaptor.define_table(table_name,
            Field('object_seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.name'), 
                  ondelete='CASCADE'),
            Field('subject_seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.seqfeature_id'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('distance','integer', ),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=[],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE seqfeature_path DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE seqfeature_path ADD UNIQUE (object_seqfeature_id, subject_seqfeature_id, term_id,distance);''')  
            self.adaptor.executesql('''CREATE INDEX seqfeaturepath_trm   ON seqfeature_path(term_id);''')
            self.adaptor.executesql('''CREATE INDEX seqfeaturepath_child ON seqfeature_path(subject_seqfeature_id);''')
            
        
        table_name = 'seqfeature_qualifier_value'
        self.adaptor.define_table(table_name,
            Field('seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.seqfeature_id'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),    
            Field('rank','integer', default = 0, notnull = True,),
            Field('value','text', notnull = True,),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['seqfeature_id','term_id', 'rank'],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list: 
                self.adaptor.executesql('''ALTER TABLE seqfeature_qualifier_value DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE seqfeature_qualifier_value ADD PRIMARY KEY (seqfeature_id, term_id, rank);''')
            self.adaptor.executesql('''CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value (term_id);''')
                    
        table_name = 'seqfeature_dbxref'
        self.adaptor.define_table(table_name,
            Field('seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.seqfeature_id'), 
                  ondelete='CASCADE'),
            Field('dbxref_id',
                  'reference dbxref',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id','dbxref.accession'), 
                  ondelete='CASCADE'),    
            Field('rank','integer',),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['seqfeature_id','dbxref_id', ],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list:
                self.adaptor.executesql('''ALTER TABLE seqfeature_dbxref DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE seqfeature_dbxref ADD PRIMARY KEY (seqfeature_id, dbxref_id);''')
            self.adaptor.executesql('''CREATE INDEX feadblink_dbx  ON seqfeature_dbxref(dbxref_id);''')
           






        table_name = 'location'
        self.adaptor.define_table(table_name,
            Field('location_id','id',  unique = True),
            Field('seqfeature_id',
                  'reference seqfeature',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'seqfeature.seqfeature_id','seqfeature.seqfeature_id'), 
                  ondelete='CASCADE'),
            Field('dbxref_id',
                  'reference dbxref',
                  requires=IS_IN_DB(self.adaptor,'dbxref.dbxref_id','dbxref.accession'), 
                  ondelete='NO ACTION'),
            Field('term_id',
                  'reference term',
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),
            Field('start_pos','integer',),
            Field('end_pos','integer',),
            Field('rank','integer', default = 0, notnull = True,),    
            Field('strand','integer', default = 0, notnull = True,),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            sequence_name = get_sequence_name(table_name))
        if create_tables:
            if self.dbtype in alter_db_list:
                self.adaptor.executesql('''ALTER TABLE location ADD UNIQUE (seqfeature_id, rank);''')
            self.adaptor.executesql('''CREATE INDEX seqfeatureloc_start ON location(start_pos, end_pos);''')
            self.adaptor.executesql('''CREATE INDEX seqfeatureloc_dbx   ON location(dbxref_id);''')
            self.adaptor.executesql('''CREATE INDEX seqfeatureloc_trm   ON location(term_id);''')
            
        table_name = 'location_qualifier_value'
        self.adaptor.define_table(table_name,
            Field('location_id',
                  'reference location',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'location.location_id','location.location_id'), 
                  ondelete='CASCADE'),
            Field('term_id',
                  'reference term',
                  notnull = True,
                  requires=IS_IN_DB(self.adaptor,'term.term_id','term.name'), 
                  ondelete='NO ACTION'),    
            Field('int_value','integer'),
            Field('value','string', length=255, notnull = True,),
            migrate = create_tables,
            fake_migrate = FAKE_MIGRATE,
            #primarykey=['location_id','dbxref_id', ],
            sequence_name = None)
        if create_tables:
            if self.dbtype in alter_db_list:
                self.adaptor.executesql('''ALTER TABLE location_qualifier_value DROP COLUMN id;''')
                #self.adaptor.executesql('''DROP SEQUENCE %s;''' % self.adaptor[table_name]._sequence_name)
                self.adaptor.executesql('''ALTER TABLE location_qualifier_value ADD PRIMARY KEY (location_id, term_id);''')
            self.adaptor.executesql('''CREATE INDEX locationqual_trm ON location_qualifier_value(term_id);''')


        
        

        
 
        
        
        
    #========================== LOADER METHODS =====================================
        
    def _get_ontology_id(self, name, definition=None):
        """Returns the identifier for the named ontology (PRIVATE).

        This looks through the onotology table for a the given entry name.
        If it is not found, a row is added for this ontology (using the
        definition if supplied).  In either case, the id corresponding to
        the provided name is returned, so that you can reference it in
        another table.
        """
        oids = self.adaptor(self.adaptor.ontology.name == name).select(self.adaptor.ontology.ontology_id)
        if oids:
            return oids[0].ontology_id
        oid = self.adaptor.ontology.insert(name = name, definition = definition)
        self.adaptor.commit()
        return oid

    def _get_term_id(self,
                     name,
                     ontology_id=None,
                     definition=None,
                     identifier=None):
        """Get the id that corresponds to a term (PRIVATE).

        This looks through the term table for a the given term. If it
        is not found, a new id corresponding to this term is created.
        In either case, the id corresponding to that term is returned, so
        that you can reference it in another table.

        The ontology_id should be used to disambiguate the term.
        """

        # try to get the term id
        query = (self.adaptor.term.name == name)
        if ontology_id:
            query=(query & (self.adaptor.term.ontology_id == ontology_id))
        
        results = self.adaptor(query).select(self.adaptor.term.term_id)
        # something is wrong
        if len(results) > 1:
            raise ValueError("Multiple term ids for %s: %r" % 
                             (name, [row.term_id for row in results]))
        elif len(results) == 1:
            return results[0].term_id
        else:
            oid = self.adaptor.term.insert(name = name, 
                                           definition = definition, 
                                           identifier = identifier,
                                           ontology_id = ontology_id) 
            self.adaptor.commit()
            return oid


        
    def _add_dbxref(self, dbname, accession, version):
        """Insert a dbxref and return its id."""
        
        oid = self.adaptor.dbxref.insert(dbname = dbname,
                                         accession = accession,
                                         version = version)
        return oid

           
    def _get_taxon_id(self, record):
        """Get the taxon id for this record (PRIVATE).

        record - a SeqRecord object

        This searches the taxon/taxon_name tables using the
        NCBI taxon ID, scientific name and common name to find
        the matching taxon table entry's id.
        
        If the species isn't in the taxon table, and we have at
        least the NCBI taxon ID, scientific name or common name,
        at least a minimal stub entry is created in the table.

        Returns the taxon id (database key for the taxon table,
        not an NCBI taxon ID), or None if the taxonomy information
        is missing.

        See also the BioSQL script load_ncbi_taxonomy.pl which
        will populate and update the taxon/taxon_name tables
        with the latest information from the NCBI.
        """
        
        # To find the NCBI taxid, first check for a top level annotation
        ncbi_taxon_id = None
        if "ncbi_taxid" in record.annotations:
            #Could be a list of IDs.
            if isinstance(record.annotations["ncbi_taxid"],list):
                if len(record.annotations["ncbi_taxid"])==1:
                    ncbi_taxon_id = record.annotations["ncbi_taxid"][0]
            else:
                ncbi_taxon_id = record.annotations["ncbi_taxid"]
        if not ncbi_taxon_id:
            # Secondly, look for a source feature
            for f in record.features:
                if f.type == 'source':
                    quals = getattr(f, 'qualifiers', {})
                    if "db_xref" in quals:
                        for db_xref in f.qualifiers["db_xref"]:
                            if db_xref.startswith("taxon:"):
                                ncbi_taxon_id = int(db_xref[6:])
                                break
                if ncbi_taxon_id: break

        try:
            scientific_name = record.annotations["organism"][:255]
        except KeyError:
            scientific_name = None
        try:
            common_name = record.annotations["source"][:255]
        except KeyError:
            common_name = None
        # Note: The maximum length for taxon names in the schema is 255.
        # Cropping it now should help in getting a match when searching,
        # and avoids an error if we try and add these to the database.


        if ncbi_taxon_id:
            #Good, we have the NCBI taxon to go on - this is unambiguous :)
            #Note that the scientific name and common name will only be
            #used if we have to record a stub entry.
            return self._get_taxon_id_from_ncbi_taxon_id(ncbi_taxon_id,
                                                         scientific_name,
                                                         common_name)
        
        if not common_name and not scientific_name:
            # Nothing to go on... and there is no point adding
            # a new entry to the database.  We'll just leave this
            # sequence's taxon as a NULL in the database.
            return None

        # Next, we'll try to find a match based on the species name
        # (stored in GenBank files as the organism and/or the source).
        if scientific_name:
            taxa = self.adaptor(self.adaptor.taxon_name.name == scientific_name).select()
            if taxa:
                #Good, mapped the scientific name to a taxon table entry
                return taxa[0].taxon_id

        # Last chance...
        if common_name:
            taxa = self.adaptor(self.adaptor.taxon_name.name == common_name).select(groupby = self.adaptor.taxon_id)

            #Its natural that several distinct taxa will have the same common
            #name - in which case we can't resolve the taxon uniquely.
            if len(taxa) > 1:
                raise ValueError("Taxa: %d species have name %r" % (
                    len(taxa),
                    common_name))
            if taxa:
                #Good, mapped the common name to a taxon table entry
                return taxa[0].taxon_id

        # At this point, as far as we can tell, this species isn't
        # in the taxon table already.  So we'll have to add it.
        # We don't have an NCBI taxonomy ID, so if we do record just
        # a stub entry, there is no simple way to fix this later.
        #
        # TODO - Should we try searching the NCBI taxonomy using the
        # species name?
        #
        # OK, let's try inserting the species.
        # Chances are we don't have enough information ...
        # Furthermore, it won't be in the hierarchy.

        lineage = []
        for c in record.annotations.get("taxonomy", []):
            lineage.append([None, None, c])
        if lineage:
            lineage[-1][1] = "genus"
        lineage.append([None, "species", record.annotations["organism"]])
        # XXX do we have them?
        if "subspecies" in record.annotations:
            lineage.append([None, "subspecies",
                            record.annotations["subspecies"]])
        if "variant" in record.annotations:
            lineage.append([None, "varietas",
                            record.annotations["variant"]])
        lineage[-1][0] = ncbi_taxon_id
        
        left_value = self.adaptor(self.adaptor.taxon.taxon_id > 0).select(self.adaptor.taxon.left_value, 
                                                                          orderby = ~self.adaptor.taxon.left_value)
        if left_value:
            left_value =left_value[0].left_value + 1
        else:
            left_value = 0
        
        # XXX -- Brad: Fixing this for now in an ugly way because
        # I am getting overlaps for right_values. I need to dig into this
        # more to actually understand how it works. I'm not sure it is
        # actually working right anyhow.
        right_start_value = self.adaptor(self.adaptor.taxon.taxon_id > 0).select(self.adaptor.taxon.right_value, 
                                                                          orderby = ~self.adaptor.taxon.right_value)
        if right_start_value:
            right_start_value = right_start_value[0].right_value
        else:
            right_start_value = 0
        right_value = right_start_value + 2 * len(lineage) - 1

        parent_taxon_id = None
        for taxon in lineage:
            taxon_id = self.adaptor.taxon.insert(parent_taxon_id = parent_taxon_id,
                                                 ncbi_taxon_id = taxon[0], 
                                                 node_rank = taxon[1],
                                                 left_value = left_value, 
                                                 right_value = right_value)

            taxon_name_id = self.adaptor.taxon_name.insert(taxon_id = taxon_id,
                                                           name = taxon[2][:255],
                                                           name_class = 'scientific name')
            #Note the name field is limited to 255, some SwissProt files
            #have a multi-species name which can be longer.  So truncate this.
            left_value += 1
            right_value -= 1
            parent_taxon_id = taxon_id
        if common_name:
            taxon_name_id = self.adaptor.taxon_name.insert(taxon_id = taxon_id,
                                                           name = common_name,
                                                           name_class = 'common name')

        return taxon_id

    def _fix_name_class(self, entrez_name):
        """Map Entrez name terms to those used in taxdump (PRIVATE).

        We need to make this conversion to match the taxon_name.name_class
        values used by the BioSQL load_ncbi_taxonomy.pl script.
        
        e.g.
        "ScientificName" -> "scientific name",
        "EquivalentName" -> "equivalent name",
        "Synonym" -> "synonym",
        """
        #Add any special cases here:
        #
        #known = {}
        #try:
        #    return known[entrez_name]
        #except KeyError:
        #    pass

        #Try automatically by adding spaces before each capital
        def add_space(letter):
            if letter.isupper():
                return " "+letter.lower()
            else:
                return letter
        answer = "".join([add_space(letter) for letter in entrez_name]).strip()
        assert answer == answer.lower()
        return answer

    def _get_taxon_id_from_ncbi_taxon_id(self, ncbi_taxon_id,
                                         scientific_name = None,
                                         common_name = None):
        """Get the taxon id for this record from the NCBI taxon ID (PRIVATE).

        ncbi_taxon_id - string containing an NCBI taxon id
        scientific_name - string, used if a stub entry is recorded
        common_name - string, used if a stub entry is recorded
        
        This searches the taxon table using ONLY the NCBI taxon ID
        to find the matching taxon table entry's ID (database key).
        
        If the species isn't in the taxon table, and the fetch_NCBI_taxonomy
        flag is true, Biopython will attempt to go online using Bio.Entrez
        to fetch the official NCBI lineage, recursing up the tree until an
        existing entry is found in the database or the full lineage has been
        fetched.

        Otherwise the NCBI taxon ID, scientific name and common name are
        recorded as a minimal stub entry in the taxon and taxon_name tables.
        Any partial information about the lineage from the SeqRecord is NOT
        recorded.  This should mean that (re)running the BioSQL script
        load_ncbi_taxonomy.pl can fill in the taxonomy lineage.

        Returns the taxon id (database key for the taxon table, not
        an NCBI taxon ID).
        """
        assert ncbi_taxon_id

        taxon_id = self.adaptor(self.adaptor.taxon.taxon_id == ncbi_taxon_id).select()
        if taxon_id:
            #Good, we have mapped the NCBI taxid to a taxon table entry
            return taxon_id[0].taxon_id

        # At this point, as far as we can tell, this species isn't
        # in the taxon table already.  So we'll have to add it.

        parent_taxon_id = None
        rank = "species"
        genetic_code = None
        mito_genetic_code = None
        species_names = []
        if scientific_name:
            species_names.append(("scientific name", scientific_name))
        if common_name:
            species_names.append(("common name", common_name))
        
        if self.fetch_NCBI_taxonomy:
            #Go online to get the parent taxon ID!
            handle = Entrez.efetch(db="taxonomy",id=ncbi_taxon_id,retmode="XML")
            taxonomic_record = Entrez.read(handle)
            if len(taxonomic_record) == 1:
                assert taxonomic_record[0]["TaxId"] == str(ncbi_taxon_id), \
                       "%s versus %s" % (taxonomic_record[0]["TaxId"],
                                         ncbi_taxon_id)
                parent_taxon_id = self._get_taxon_id_from_ncbi_lineage( \
                                            taxonomic_record[0]["LineageEx"])
                rank = taxonomic_record[0]["Rank"]
                genetic_code = taxonomic_record[0]["GeneticCode"]["GCId"]
                mito_genetic_code = taxonomic_record[0]["MitoGeneticCode"]["MGCId"]
                species_names = [("scientific name",
                                  taxonomic_record[0]["ScientificName"])]
                try:
                    for name_class, names in taxonomic_record[0]["OtherNames"].iteritems():
                        name_class = self._fix_name_class(name_class)
                        if not isinstance(names, list):
                            #The Entrez parser seems to return single entry
                            #lists as just a string which is annoying.
                            names = [names]
                        for name in names:
                            #Want to ignore complex things like ClassCDE entries
                            if isinstance(name, basestring):
                                species_names.append((name_class, name))
                except KeyError:
                    #OtherNames isn't always present,
                    #e.g. NCBI taxon 41205, Bromheadia finlaysoniana
                    pass
        else:
            pass
            # If we are not allowed to go online, we will record the bare minimum;
            # as long as the NCBI taxon id is present, then (re)running
            # load_ncbi_taxonomy.pl should fill in the taxonomomy lineage
            # (and update the species names).
            #
            # I am NOT going to try and record the lineage, even if it
            # is in the record annotation as a list of names, as we won't
            # know the NCBI taxon IDs for these parent nodes.
        
        taxon_id = self.adaptor.taxon.insert(parent_taxon_id = parent_taxon_id, 
                                             ncbi_taxon_id = ncbi_taxon_id, 
                                             node_rank = rank,
                                             genetic_code = genetic_code, 
                                             mito_genetic_code = mito_genetic_code, 
                                             left_value = None, 
                                             right_value = None)


        #Record the scientific name, common name, etc
        for name_class, name in species_names:
            taxon_name_id = self.adaptor.taxon_name.insert(taxon_id = taxon_id,
                                                           name = name[:255],
                                                           name_class = name_class)

        return taxon_id

    def _get_taxon_id_from_ncbi_lineage(self, taxonomic_lineage):# TODO
        """This is recursive! (PRIVATE).

        taxonomic_lineage - list of taxonomy dictionaries from Bio.Entrez

        First dictionary in list is the taxonomy root, highest would be the species.
        Each dictionary includes:
        - TaxID (string, NCBI taxon id)
        - Rank (string, e.g. "species", "genus", ..., "phylum", ...)
        - ScientificName (string)
        (and that is all at the time of writing)

        This method will record all the lineage given, returning the the taxon id
        (database key, not NCBI taxon id) of the final entry (the species).
        """
        ncbi_taxon_id = taxonomic_lineage[-1]["TaxId"]
        
        
        #Is this in the database already?  Check the taxon table...
        taxon_id = self.adaptor(self.adaptor.taxon.taxon_id == ncbi_taxon_id).select()
        if taxon_id:
            # we could verify that the Scientific Name etc in the database
            # is the same and update it or print a warning if not...
            return taxon_id[0].taxon_id

        #We have to record this.
        if len(taxonomic_lineage) > 1:
            #Use recursion to find out the taxon id (database key) of the parent.
            parent_taxon_id = self._get_taxon_id_from_ncbi_lineage(taxonomic_lineage[:-1])
            assert isinstance(parent_taxon_id, int) or isinstance(parent_taxon_id, long), repr(parent_taxon_id)
        else:
            parent_taxon_id = None

        # INSERT new taxon
        rank = taxonomic_lineage[-1].get("Rank", None)
        taxon_id = self.adaptor.taxon.insert(parent_taxon_id = parent_taxon_id, 
                                             ncbi_taxon_id = ncbi_taxon_id, 
                                             node_rank = rank,)
        assert isinstance(taxon_id, int) or isinstance(taxon_id, long), repr(taxon_id)
        # ... and its name in taxon_name
        scientific_name = taxonomic_lineage[-1].get("ScientificName", None)
        if scientific_name:
            taxon_name_id = self.adaptor.taxon_name.insert(taxon_id = taxon_id,
                                                           name = scientific_name[:255],
                                                           name_class = 'scientific name')
        return taxon_id

    def _get_accession_from_seqrecord(self, record):
        '''defines accession and version given a seqrecord '''
        if record.id.count(".") == 1: # try to get a version from the id
            #This assumes the string is something like "XXXXXXXX.123"
            accession, version = record.id.split('.')
            try:
                version = int(version)
            except ValueError:
                accession = record.id
                version = 0
        else: # otherwise just use a version of 0
            accession = record.id
            version = 0

        if "accessions" in record.annotations \
        and isinstance(record.annotations["accessions"],list) \
        and record.annotations["accessions"]:
            #Take the first accession (one if there is more than one)
            accession = record.annotations["accessions"][0]
            
        return accession, version

    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information (PRIVATE).

        record - SeqRecord object to add to the database.
        """
        # get the pertinent info and insert it
            
        accession,version = self._get_accession_from_seqrecord(record)

        #Find the taxon id (this is not just the NCBI Taxon ID)
        #NOTE - If the species isn't defined in the taxon table,
        #a new minimal entry is created.
        taxon_id = self._get_taxon_id(record)#TO DO

        if "gi" in record.annotations:
            identifier = record.annotations["gi"]
        else:
            if len(record.id) <= 40:
                identifier = record.id
            else:
                identifier = None

        #Allow description and division to default to NULL as in BioPerl.
        description = getattr(record, 'description', None)
        division = record.annotations.get("data_file_division", None)
        
        
        bioentry_id = self.adaptor.bioentry.insert(biodatabase_id = self.dbid,
                                                   taxon_id = taxon_id,
                                                   name = record.name, 
                                                   accession = accession,
                                                   identifier = identifier,
                                                   division = division,
                                                   description = description,
                                                   version = version)

        return bioentry_id

    def _load_bioentry_date(self, record, bioentry_id):
        """Add the effective date of the entry into the database.

        record - a SeqRecord object with an annotated date
        bioentry_id - corresponding database identifier
        """
        # dates are GenBank style, like:
        # 14-SEP-2000
        date = record.annotations.get("date",
                                      strftime("%d-%b-%Y", gmtime()).upper())
        if isinstance(date, list) : date = date[0]
        annotation_tags_id = self._get_ontology_id("Annotation Tags")
        date_id = self._get_term_id("date_changed", annotation_tags_id)
        
        date_oid = self.adaptor.bioentry_qualifier_value.insert(bioentry_id = bioentry_id, 
                                                                term_id = date_id, 
                                                                value = date, 
                                                                rank = 1)
        


    def _load_biosequence(self, record, bioentry_id):
        """Record a SeqRecord's sequence and alphabet in the database (PRIVATE).

        record - a SeqRecord object with a seq property
        bioentry_id - corresponding database identifier
        """
        if record.seq is None:
            #The biosequence table entry is optional, so if we haven't
            #got a sequence, we don't need to write to the table.
            return
        
        # determine the string representation of the alphabet
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            alphabet = "dna"
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            alphabet = "rna"
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            alphabet = "protein"
        else:
            alphabet = "unknown"

        if isinstance(record.seq, UnknownSeq):
            seq_str = None
        else:
            seq_str = str(record.seq)

        sequence_oid = self.adaptor.biosequence.insert(bioentry_id = bioentry_id, 
                                                       version = 0, 
                                                       length = len(record.seq), 
                                                       seq = seq_str,
                                                       alphabet = alphabet,)


    def _load_comment(self, record, bioentry_id):
        """Record a SeqRecord's annotated comment in the database (PRIVATE).

        record - a SeqRecord object with an annotated comment
        bioentry_id - corresponding database identifier
        """
        comments = record.annotations.get('comment')
        if not comments:
            return
        if not isinstance(comments, list):
            #It should be a string then...
            comments = [comments]

        for index, comment in enumerate(comments):
            comment = comment.replace('\n', ' ')
            #TODO - Store each line as a separate entry?  This would preserve
            #the newlines, but we should check BioPerl etc to be consistent.
            
            sequence_oid = self.adaptor.comment.insert(bioentry_id = bioentry_id, 
                                                       comment_text = comment,
                                                       rank = index+1)
            
            
        
    def _load_annotations(self, record, bioentry_id):
        """Record a SeqRecord's misc annotations in the database (PRIVATE).

        The annotation strings are recorded in the bioentry_qualifier_value
        table, except for special cases like the reference, comment and
        taxonomy which are handled with their own tables.

        record - a SeqRecord object with an annotations dictionary
        bioentry_id - corresponding database identifier
        
        #===============================#
        WARNING:
        in biopython if a string is passed no rank is added (leave to db default). 
        this is changed in BioSQLAnnotation
        the rank is needed to discriminate between two duplicated couples o term and pair,
        and better cope with the absence of a primary key, even if this will not resolve the issue.
        also no DB default will be set when creating the db with web2py dal (possible bug) so check for consistence.
        in that case however an id field is available to discriminate between duplicated couples.
        if possible always pass an annotation as a list. or use the BioSQLAnnotation for adding annotations.
        
        UPDATE:
        for consistence a rank=0 is set even for string annotations.
        use this function just for the fist upload of a seqrecord.
        #===============================#
        """
        tag_ontology_id = self._get_ontology_id('Annotation Tags')
        for key, value in record.annotations.iteritems():
            if key in ["references", "comment", "ncbi_taxid"]:
                #Handled separately
                continue
            elif key.startswith('relationship_') and not self.comp_mode:
                self._load_relation_from_annotation(bioentry_id, key, value)
            term_id = self._get_term_id(key, ontology_id=tag_ontology_id)
            if isinstance(value, list) or isinstance(value, tuple):
                rank = 0
                for entry in value:
                    if isinstance(entry, str) or isinstance(entry, int):
                        #Easy case
                        rank += 1
                        self.adaptor.bioentry_qualifier_value.insert(bioentry_id = bioentry_id,
                                                                     term_id =term_id,
                                                                     value = str(entry),
                                                                     rank = rank)
                    else:
                        pass
                        #print "Ignoring annotation '%s' sub-entry of type '%s'" \
                        #      % (key, str(type(entry)))
            elif isinstance(value, str) or isinstance(value, int):
                #Have a simple single entry, leave rank as the DB default# modified: set rank = 0, differs from biopython
                rank = 0
                self.adaptor.bioentry_qualifier_value.insert(bioentry_id = bioentry_id,
                                                             term_id =term_id,
                                                             value = str(value),
                                                             rank = rank)
            else:
                pass
                #print "Ignoring annotation '%s' entry of type '%s'" \
                #      % (key, type(value))


    def _load_reference(self, reference, rank, bioentry_id):
        """Record a SeqRecord's annotated references in the database (PRIVATE).
        record - a SeqRecord object with annotated references
        bioentry_id - corresponding database identifier
        """

        refs = None
        if reference.medline_id:
            refs = self.adaptor((self.adaptor.reference.dbxref_id == self.adaptor.dbxref.dbxref_id) &
                                (self.adaptor.dbxref.dbname == 'MEDLINE') &
                                (self.adaptor.dbxref.accession == reference.medline_id)).select(self.adaptor.reference.reference_id)

        if not refs and reference.pubmed_id:
            refs = self.adaptor((self.adaptor.reference.dbxref_id == self.adaptor.dbxref.dbxref_id) &
                                (self.adaptor.dbxref.dbname == 'PUBMED') &
                                (self.adaptor.dbxref.accession == reference.pubmed_id)).select(self.adaptor.reference.reference_id)
        if not refs:
            s = []
            for f in reference.authors, reference.title, reference.journal:
                s.append(f or "<undef>")
            crc = crc64("".join(s))
            refs = self.adaptor(self.adaptor.reference.crc == crc).select(self.adaptor.reference.reference_id)
        if not refs:
            if reference.medline_id:
                dbxref_id = self._add_dbxref("MEDLINE",
                                             reference.medline_id, 0)
            elif reference.pubmed_id:
                dbxref_id = self._add_dbxref("PUBMED",
                                             reference.pubmed_id, 0)
            else:
                dbxref_id = None
            authors = reference.authors or None
            title =  reference.title or None
            #The location/journal field cannot be Null, so default
            #to an empty string rather than None:
            journal = reference.journal or ""
            reference_id = self.adaptor.reference.insert(dbxref_id = dbxref_id, 
                                                         location = journal,
                                                         title = title, 
                                                         authors = authors, 
                                                         crc = crc)
        else:
            reference_id = refs[0].reference_id

        if reference.location:
            start = 1 + int(str(reference.location[0].start))
            end = int(str(reference.location[0].end))
        else:
            start = None
            end = None
            
        self.adaptor.bioentry_reference.insert(bioentry_id = bioentry_id, 
                                               reference_id = reference_id,
                                               start_pos = start, 
                                               end_pos = end, 
                                               rank = rank + 1)
        
        
    def _load_seqfeature(self, feature, feature_rank, bioentry_id): 
        """Load a biopython SeqFeature into the database (PRIVATE).
        """
        seqfeature_id = self._load_seqfeature_basic(feature.type, feature_rank,
                                                    bioentry_id)
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

    def _load_seqfeature_basic(self, feature_type, feature_rank, bioentry_id):
        """Load the first tables of a seqfeature and returns the id (PRIVATE).

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """
        ontology_id = self._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self._get_term_id(feature_type,
                                              ontology_id = ontology_id)
        # XXX source is always EMBL/GenBank/SwissProt here; it should depend on
        # the record (how?)
        source_cat_id = self._get_ontology_id('SeqFeature Sources')
        source_term_id = self._get_term_id('EMBL/GenBank/SwissProt',
                                      ontology_id = source_cat_id)
        
        seqfeature_id = self.adaptor.seqfeature.insert(bioentry_id = bioentry_id, 
                                                       type_term_id = seqfeature_key_id, 
                                                       source_term_id = source_term_id, 
                                                       rank = feature_rank + 1)
        

        return seqfeature_id

    def _load_seqfeature_locations(self, feature, seqfeature_id):
        """Load all of the locations for a SeqFeature into tables (PRIVATE).

        This adds the locations related to the SeqFeature into the
        seqfeature_location table. Fuzzies are not handled right now.
        For a simple location, ie (1..2), we have a single table row
        with seq_start = 1, seq_end = 2, location_rank = 1.

        For split locations, ie (1..2, 3..4, 5..6) we would have three
        row tables with:
            start = 1, end = 2, rank = 1
            start = 3, end = 4, rank = 2
            start = 5, end = 6, rank = 3
        """
        # TODO - Record an ontology for the locations (using location.term_id)
        # which for now as in BioPerl we leave defaulting to NULL.
        if feature.location_operator and feature.location_operator != "join":
            # e.g. order locations... we don't record "order" so it
            # will become a "join" on reloading. What does BioPerl do?
            import warnings
            warnings.warn("%s location operators are not fully supported" \
                          % feature.location_operator)
        
        # two cases, a simple location or a split location
        if not feature.sub_features:    # simple location
            self._insert_seqfeature_location(feature, 1, seqfeature_id)
        else: # split location
            for rank, cur_feature in enumerate(feature.sub_features):
                self._insert_seqfeature_location(cur_feature,
                                                 rank + 1,
                                                 seqfeature_id)

    def _insert_seqfeature_location(self, feature, rank, seqfeature_id):
        """Add a location of a SeqFeature to the seqfeature_location table (PRIVATE).

        TODO - Add location_operators to location_qualifier_value.
        """
        # convert biopython locations to the 1-based location system
        # used in bioSQL
        # XXX This could also handle fuzzies
        start = feature.location.nofuzzy_start + 1
        end = feature.location.nofuzzy_end

        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        strand = feature.strand or 0

        # TODO - Record an ontology term for the location (location.term_id)
        # which for now like BioPerl we'll leave as NULL.
        # This might allow us to record "between" positions properly, but I
        # doesn't really see how it could work for before/after fuzzy positions
        loc_term_id = None

        if feature.ref:
            # sub_feature remote locations when they are in the same db as the current
            # record do not have a value for ref_db, which the SeqFeature object
            # stores as None. BioSQL schema requires a varchar and is not NULL 
            dbxref_id = self._get_dbxref_id(feature.ref_db or "", feature.ref)
        else:
            dbxref_id = None

        oid = self.adaptor.location.insert(seqfeature_id = seqfeature_id, 
                                           dbxref_id = dbxref_id, 
                                           term_id = loc_term_id,
                                           start_pos = start, 
                                           end_pos = end, 
                                           strand = strand, 
                                           rank = rank)
        

        """
        # See Bug 2677
        # TODO - Record the location_operator (e.g. "join" or "order")
        # using the location_qualifier_value table (which we and BioPerl
        # have historically left empty).
        # Note this will need an ontology term for the location qualifer
        # (location_qualifier_value.term_id) for which oddly the schema
        # does not allow NULL.
        if feature.location_operator:
            #e.g. "join" (common),
            #or "order" (see Tests/GenBank/protein_refseq2.gb)
            location_id = self.adaptor.last_id('location')
            loc_qual_term_id = None # Not allowed in BioSQL v1.0.1
            sql = r"INSERT INTO location_qualifier_value" \
                  r"(location_id, term_id, value)" \
                  r"VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (location_id, loc_qual_term_id,
                                       feature.location_operator))
        """

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature (PRIVATE).

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        tag_ontology_id = self._get_ontology_id('Annotation Tags')
        for qualifier_key in qualifiers.keys():
            # Treat db_xref qualifiers differently to sequence annotation
            # qualifiers by populating the seqfeature_dbxref and dbxref
            # tables.  Other qualifiers go into the seqfeature_qualifier_value
            # and (if new) term tables.
            if qualifier_key != 'db_xref':
                qualifier_key_id = self._get_term_id(qualifier_key,
                                                  ontology_id=tag_ontology_id)
                # now add all of the values to their table
                entries = qualifiers[qualifier_key]
                if not isinstance(entries, list):
                    # Could be a plain string, or an int or a float.
                    # However, we exect a list of strings here.
                    entries = [entries]
                for qual_value_rank in range(len(entries)):
                    qualifier_value = entries[qual_value_rank]
                    oid = self.adaptor.seqfeature_qualifier_value.insert(seqfeature_id = seqfeature_id, 
                                                                          term_id = qualifier_key_id, 
                                                                          rank = qual_value_rank + 1, 
                                                                          value = qualifier_value)

            else:
                # The dbxref_id qualifier/value sets go into the dbxref table
                # as dbname, accession, version tuples, with dbxref.dbxref_id
                # being automatically assigned, and into the seqfeature_dbxref
                # table as seqfeature_id, dbxref_id, and rank tuples
                self._load_seqfeature_dbxref(qualifiers[qualifier_key],
                                             seqfeature_id)


    def _load_seqfeature_dbxref(self, dbxrefs, seqfeature_id):
        """Add database crossreferences of a SeqFeature to the database (PRIVATE).

            o dbxrefs           List, dbxref data from the source file in the
                                format <database>:<accession>

            o seqfeature_id     Int, the identifier for the seqfeature in the
                                seqfeature table

            Insert dbxref qualifier data for a seqfeature into the
            seqfeature_dbxref and, if required, dbxref tables.
            The dbxref_id qualifier/value sets go into the dbxref table
            as dbname, accession, version tuples, with dbxref.dbxref_id
            being automatically assigned, and into the seqfeature_dbxref
            table as seqfeature_id, dbxref_id, and rank tuples
        """
        # NOTE - In older versions of Biopython, we would map the GenBank
        # db_xref "name", for example "GI" to "GeneIndex", and give a warning
        # for any unknown terms.  This was a long term maintainance problem,
        # and differed from BioPerl and BioJava's implementation.  See bug 2405
        for rank, value in enumerate(dbxrefs):
            # Split the DB:accession format string at colons.  We have to
            # account for multiple-line and multiple-accession entries
            try:
                dbxref_data = value.replace(' ','').replace('\n','').split(':')
                db = dbxref_data[0]
                accessions = dbxref_data[1:]
            except:
                raise ValueError("Parsing of db_xref failed: '%s'" % value)
            # Loop over all the grabbed accessions, and attempt to fill the
            # table
            for accession in accessions:
                # Get the dbxref_id value for the dbxref data
                dbxref_id = self._get_dbxref_id(db, accession)
                # Insert the seqfeature_dbxref data
                self._get_seqfeature_dbxref(seqfeature_id, dbxref_id, rank+1)
        
    def _get_dbxref_id(self, db, accession):
        """ _get_dbxref_id(self, db, accession) -> Int

            o db          String, the name of the external database containing
                          the accession number

            o accession   String, the accession of the dbxref data

            Finds and returns the dbxref_id for the passed data.  The method
            attempts to find an existing record first, and inserts the data
            if there is no record.
        """
        # Check for an existing record
        dbxref_id = self.adaptor((self.adaptor.dbxref.dbname == db) &
                                 (self.adaptor.dbxref.accession == accession)).select(self.adaptor.dbxref.dbxref_id)
        # If there was a record, return the dbxref_id, else create the
        # record and return the created dbxref_id
        if dbxref_id:
            return dbxref_id[0].dbxref_id
        return self._add_dbxref(db, accession, 0)

    def _get_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Check for a pre-existing seqfeature_dbxref entry with the passed
            seqfeature_id and dbxref_id.  If one does not exist, insert new
            data

        """
        # Check for an existing record
        result = self.adaptor((self.adaptor.seqfeature_dbxref.seqfeature_id == seqfeature_id) &
                                 (self.adaptor.seqfeature_dbxref.dbxref_id == dbxref_id)).select(self.adaptor.seqfeature_dbxref.seqfeature_id,
                                                                                                self.adaptor.seqfeature_dbxref.dbxref_id)
        # If there was a record, return without executing anything, else create
        # the record and return
        if result:
            return (result[0].seqfeature_id, result[0].dbxref_id)
        return self._add_seqfeature_dbxref(seqfeature_id, dbxref_id, rank)

    def _add_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Insert a seqfeature_dbxref row and return the seqfeature_id and
            dbxref_id
        """
        
        self.adaptor.seqfeature_dbxref.insert(seqfeature_id = seqfeature_id, 
                                              dbxref_id = dbxref_id, 
                                              rank = rank)
        return (seqfeature_id, dbxref_id)

    def _load_dbxrefs(self, record, bioentry_id):
        """Load any sequence level cross references into the database (PRIVATE).

        See table bioentry_dbxref."""
        for rank, value in enumerate(record.dbxrefs):
            # Split the DB:accession string at first colon.
            # We have to cope with things like:
            # "MGD:MGI:892" (db="MGD", accession="MGI:892")
            # "GO:GO:123" (db="GO", accession="GO:123")
            #
            # Annoyingly I have seen the NCBI use both the style
            # "GO:GO:123" and "GO:123" in different vintages.
            assert value.count("\n")==0
            try:
                db, accession = value.split(':',1)
                db = db.strip()
                accession = accession.strip()
            except:
                raise ValueError("Parsing of dbxrefs list failed: '%s'" % value)
            # Get the dbxref_id value for the dbxref data
            dbxref_id = self._get_dbxref_id(db, accession)
            # Insert the bioentry_dbxref  data
            self._get_bioentry_dbxref(bioentry_id, dbxref_id, rank+1)

    def _get_bioentry_dbxref(self, bioentry_id, dbxref_id, rank):
        """ Check for a pre-existing bioentry_dbxref entry with the passed
            seqfeature_id and dbxref_id.  If one does not exist, insert new
            data

        """
        # Check for an existing record
        result = self.adaptor((self.adaptor.bioentry_dbxref.bioentry_id == bioentry_id) &
                                 (self.adaptor.bioentry_dbxref.dbxref_id == dbxref_id)).select(self.adaptor.bioentry_dbxref.bioentry_id,
                                                                                               self.adaptor.bioentry_dbxref.dbxref_id)

        # If there was a record, return without executing anything, else create
        # the record and return
        if result:
            return (result[0].bioentry_id, result[0].dbxref_id)
        return self._add_bioentry_dbxref(bioentry_id, dbxref_id, rank)

    def _add_bioentry_dbxref(self, bioentry_id, dbxref_id, rank):
        """ Insert a bioentry_dbxref row and return the seqfeature_id and
            dbxref_id
        """
        
        self.adaptor.bioentry_dbxref.insert(bioentry_id = bioentry_id, 
                                            dbxref_id = dbxref_id, 
                                            rank = rank)
        return (bioentry_id, dbxref_id)

    #========================== READER METHODS =====================================

    def _retrieve_taxon(self,bioentry_id, taxon_id):
        a = {}
        common_names = [row.name for row in self.adaptor((self.adaptor.taxon_name.taxon_id==taxon_id) & \
                                                   (self.adaptor.taxon_name.name_class=='genbank common name')).select(self.adaptor.taxon_name.name)]
        if common_names:
            a['source'] = common_names[0]

        scientific_names = [row.name for row in self.adaptor((self.adaptor.taxon_name.taxon_id==taxon_id) & \
                                                   (self.adaptor.taxon_name.name_class=='scientific name')).select(self.adaptor.taxon_name.name)]
        if scientific_names:
            a['organism'] = scientific_names[0]
        
        ncbi_taxids = [row.ncbi_taxon_id for row in self.adaptor(self.adaptor.taxon.taxon_id==taxon_id).select(self.adaptor.taxon.ncbi_taxon_id)]
        if ncbi_taxids and ncbi_taxids[0] and ncbi_taxids[0] != "0":
            a['ncbi_taxid'] = ncbi_taxids[0]

        taxonomy = []
        while taxon_id:
            name, rank, parent_taxon_id = [(row.taxon_name.name,row.taxon.node_rank, row.taxon.parent_taxon_id)  for \
                                            row in self.adaptor((self.adaptor.taxon.taxon_id==taxon_id) & \
                                             (self.adaptor.taxon.taxon_id==self.adaptor.taxon_name.taxon_id) & \
                                             (self.adaptor.taxon_name.name_class=='scientific name')).select(
                                                self.adaptor.taxon_name.name, self.adaptor.taxon.node_rank, self.adaptor.taxon.parent_taxon_id)][0]
            if taxon_id == parent_taxon_id:
                break
            if rank != "no rank":
                taxonomy.insert(0, name)
            taxon_id = parent_taxon_id
    
        if taxonomy:
            a['taxonomy'] = taxonomy
        return a
    
    def _retrieve_dbxrefs(self,bioentry_id):
        _dbxrefs = []
        dbxrefs=((dbxref.dbname, dbxref.accession, dbxref.version) for dbxref in self.adaptor((self.adaptor.bioentry.bioentry_id == bioentry_id) & \
                             (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) & \
                             (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id)).select( 
                                                        self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,self.adaptor.dbxref.version, orderby=self.adaptor.bioentry_dbxref.rank))
        for dbname, accession, version in dbxrefs:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            _dbxrefs.append("%s:%s" % (dbname, v))
        return _dbxrefs

    def _retrieve_comment(self,bioentry_id):
        comments = [row.comment_text for row in self.adaptor(self.adaptor.comment.bioentry_id == bioentry_id).select(orderby=self.adaptor.comment.rank)]
        if comments:
            return {"comment": comments}
        else:
            return {}

    def _retrieve_reference(self,bioentry_id):
        # XXX dbxref_qualifier_value

        '''web2py has a bug with left joins putting left joins before joins in the generated SQL. this is allowed (even if incorrect) by any db-backend except postgres. 
        making a double query and merging the dbname when possible as a workaround here here
        single query working code should be:
        
        refs=self.adaptor((self.adaptor.bioentry_reference.bioentry_id==bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id==self.adaptor.reference.reference_id)
                                              )._select(orderby=self.adaptor.bioentry_reference.rank, left=self.adaptor.reference.on(self.adaptor.reference.dbxref_id==self.adaptor.dbxref.dbxref_id))
        '''

        refs=((row.bioentry_reference.start_pos,
               row.bioentry_reference.end_pos,
               row.reference.location,
               row.reference.title,
               row.reference.authors,
               row.reference.dbxref_id,) for row in self.adaptor((self.adaptor.bioentry_reference.bioentry_id==bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id==self.adaptor.reference.reference_id) 
                                              ).select(self.adaptor.reference.reference_id,
                                                       self.adaptor.reference.location,
                                                       self.adaptor.reference.title,
                                                       self.adaptor.reference.authors,
                                                       self.adaptor.reference.dbxref_id,
                                                       self.adaptor.bioentry_reference.start_pos,
                                                       self.adaptor.bioentry_reference.end_pos,
                                                       orderby=self.adaptor.bioentry_reference.rank,))
        refs_dbxrefs=dict(((row.reference.dbxref_id,dict(dbname=row.dbxref.dbname,accession=row.dbxref.accession)) for row in self.adaptor((self.adaptor.bioentry_reference.bioentry_id==bioentry_id) & 
                                              (self.adaptor.bioentry_reference.reference_id==self.adaptor.reference.reference_id) &
                                              (self.adaptor.reference.dbxref_id==self.adaptor.dbxref.dbxref_id)
                                              ).select(self.adaptor.reference.dbxref_id,
                                                       self.adaptor.dbxref.dbname,
                                                       self.adaptor.dbxref.accession,
                                                       orderby=self.adaptor.bioentry_reference.rank,)))
        references = []
        for start, end, location, title, authors, dbxref_id in refs:
            if dbxref_id in refs_dbxrefs:
                dbname = refs_dbxrefs[dbxref_id]['dbname']
                accession = refs_dbxrefs[dbxref_id]['accession']
            else:
                dbname,accession=None,None
            reference = SeqFeature.Reference()
            #If the start/end are missing, reference.location is an empty list
            if (start is not None) or (end is not None):
                if start is not None: start -= 1 #python counting
                reference.location = [SeqFeature.FeatureLocation(start, end)]
            #Don't replace the default "" with None.
            if authors : reference.authors = authors
            if title : reference.title = title
            reference.journal = location
            if dbname == 'PUBMED':
                reference.pubmed_id = accession
            elif dbname == 'MEDLINE':
                reference.medline_id = accession
            references.append(reference)
        if references:
            return {'references': references}
        else:
            return {}   

    def _make_unicode_into_string(self,text):
        if isinstance(text, unicode):
            return str(text)
        else :
            return text
        
    def _retrieve_qualifier_value(self,bioentry_id):

        qvs=((qualifier.term.name,qualifier.bioentry_qualifier_value.value,) \
                     for qualifier in self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == bioentry_id) & \
                              (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id)\
                              ).select(self.adaptor.term.name,self.adaptor.bioentry_qualifier_value.value,orderby=self.adaptor.bioentry_qualifier_value.rank))
        qualifiers = {}
        for name, value in qvs:
            if self.comp_mode:#do we really need this stuff?
                if name == "keyword": name = "keywords"
                elif name == "date_changed": name = "dates"
                elif name == "secondary_accession": name = "accessions"
            qualifiers.setdefault(name, []).append(value)
        return qualifiers
    
    def _retrieve_qualifier_term(self,bioentry_id, qualifier_term_name = '%'):
        'given a like query returns all qualifiers terms for a given bioentry'
        
        terms = [qualifier.name for qualifier in self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == bioentry_id) & \
                              (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) & \
                              (self.adaptor.term.name.like(qualifier_term_name)) \
                              ).select(self.adaptor.term.name,orderby=self.adaptor.bioentry_qualifier_value.rank)]
        return terms
        
        
    
    def _retrieve_single_qualifier_value(self,bioentry_id,qualifier_term_name):
        '''accepts like query'''
        
        qvs=((qualifier.term.name,qualifier.bioentry_qualifier_value.value,) \
                     for qualifier in self.adaptor((self.adaptor.bioentry_qualifier_value.bioentry_id == bioentry_id) & \
                              (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) & \
                              (self.adaptor.term.name.like(qualifier_term_name)) \
                              ).select(self.adaptor.term.name,self.adaptor.bioentry_qualifier_value.value,orderby=self.adaptor.bioentry_qualifier_value.rank))
        qualifiers = {}
        for name, value in qvs:
            qualifiers.setdefault(name, []).append(value)
        return qualifiers
        
    def _retrieve_annotations(self,bioentry_id, taxon_id):
        annotations = {}
        annotations.update(self._retrieve_qualifier_value(bioentry_id))
        annotations.update(self._retrieve_reference(bioentry_id))
        annotations.update(self._retrieve_taxon(bioentry_id, taxon_id))
        annotations.update(self._retrieve_comment(bioentry_id))
        annotations.update(self._retrieve_bioentry_relationship(bioentry_id))
        str_anns = {}
        
        for key, val in annotations.items():
            if isinstance(val, list):
                val = [self._make_unicode_into_string(x) for x in val]
            elif isinstance(val, unicode):
                val = str(val)
            str_anns[key] = val
        return str_anns
    
    def _retrieve_seq(self,bioentry_id):
        '''modified from biopython, now returns a Seq object and not a DBSeq object '''
        #The database schema ensures there will be only one matching
        #row in the table.
    
        #If an UnknownSeq was recorded, seq will be NULL,
        #but length will be populated.  This means length(seq)
        #will return None.
        
        seqs = [(row.alphabet,row.length,len(row.seq),row.seq) for row in \
               self.adaptor(self.adaptor.biosequence.bioentry_id == bioentry_id \
               ).select(self.adaptor.biosequence.bioentry_id,
                        self.adaptor.biosequence.alphabet,
                        self.adaptor.biosequence.length,
                        self.adaptor.biosequence.seq,)]
        
        if not seqs : return
        assert len(seqs) == 1        
        moltype, given_length, length, seq = seqs[0]
    
        try:
            length = int(length)
            given_length = int(length)
            assert length == given_length
            have_seq = True
        except TypeError:
            assert length is None
            assert len(seqs) == 1
            assert seq is None or seq==""
            length = int(given_length)
            have_seq = False
            #del seq
        del given_length
            
        moltype = moltype.lower() #might be upper case in database
        #We have no way of knowing if these sequences will use IUPAC
        #alphabets, and we certainly can't assume they are unambiguous!
        if moltype == "dna":
            alphabet = Alphabet.generic_dna
        elif moltype == "rna":
            alphabet = Alphabet.generic_rna
        elif moltype == "protein":
            alphabet = Alphabet.generic_protein
        elif moltype == "unknown":
            #This is used in BioSQL/Loader.py and would happen
            #for any generic or nucleotide alphabets.
            alphabet = Alphabet.single_letter_alphabet
        else:
            raise AssertionError("Unknown moltype: %s" % moltype)
    
        if have_seq:
            return Seq(seq, alphabet,)
        else:
            return UnknownSeq(length, alphabet) 
            
    def _retrieve_location_qualifier_value(self,location_id):
        value= [row.value for row in self.adaptor(self.adaptor.location_qualifier_value.location_id==location_id).select(self.adaptor.location_qualifier_value.value)]
        try:
            return value[0] 
        except IndexError:
            return ""

   
    def _retrieve_features(self,bioentry_id):
        results = ((row.seqfeature.seqfeature_id,row.term.name,row.seqfeature.rank) \
                   for row in self.adaptor((self.adaptor.seqfeature.bioentry_id == bioentry_id) & \
                                      (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id)
                                      ).select(self.adaptor.seqfeature.seqfeature_id,self.adaptor.term.name,self.adaptor.seqfeature.rank,orderby=self.adaptor.seqfeature.rank))

        seq_feature_list = []
        for seqfeature_id, seqfeature_type, seqfeature_rank in results:
            # Get qualifiers [except for db_xref which is stored separately]
            qvs=((qualifier.term.name,qualifier.seqfeature_qualifier_value.value,) \
                         for qualifier in self.adaptor((self.adaptor.seqfeature_qualifier_value.seqfeature_id == seqfeature_id) & \
                                  (self.adaptor.seqfeature_qualifier_value.term_id == self.adaptor.term.term_id)\
                                  ).select(self.adaptor.term.name,self.adaptor.seqfeature_qualifier_value.value,orderby=self.adaptor.seqfeature_qualifier_value.rank))

            qualifiers = {}
            for qv_name, qv_value in qvs:
                qualifiers.setdefault(qv_name, []).append(qv_value)
            # Get db_xrefs [special case of qualifiers]
            qvs = ((row.dbname,row.accession) \
                     for row in self.adaptor((self.adaptor.seqfeature_dbxref.seqfeature_id==seqfeature_id) &  
                                          (self.adaptor.seqfeature_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id) 
                                         ).select(self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,orderby=self.adaptor.seqfeature_dbxref.rank))
            for qv_name, qv_value in qvs:
                value = "%s:%s" % (qv_name, qv_value)
                qualifiers.setdefault("db_xref", []).append(value)
            # Get locations
            results= ((row.location_id,
                       row.start_pos,
                       row.end_pos,
                       row.strand) for row in self.adaptor(self.adaptor.location.seqfeature_id==seqfeature_id).select(self.adaptor.location.location_id,
                                                                                                              self.adaptor.location.start_pos,
                                                                                                              self.adaptor.location.end_pos,
                                                                                                              self.adaptor.location.strand,
                                                                                                              orderby=self.adaptor.location.rank))
            locations = []
            # convert to Python standard form
            # Convert strand = 0 to strand = None
            # re: comment in Loader.py:
            # Biopython uses None when we don't know strand information but
            # BioSQL requires something (non null) and sets this as zero
            # So we'll use the strand or 0 if Biopython spits out None
            for location_id, start, end, strand in results:
                if start:
                    start -= 1
                if strand == 0:
                    strand = None
                if strand not in (+1, -1, None):
                    raise ValueError("Invalid strand %s found in database for " \
                                     "seqfeature_id %s" % (strand, seqfeature_id))
                if end < start:
                    import warnings
                    warnings.warn("Inverted location start/end (%i and %i) for " \
                                  "seqfeature_id %s" % (start, end, seqfeature_id))
                locations.append( (location_id, start, end, strand) )
            # Get possible remote reference information
            remote_results = ((row.location.location_id,
                               row.dbxref.dbname,
                               row.dbxref.accession,
                               row.dbxref.accession,) for row in self.adaptor((self.adaptor.location.seqfeature_id == seqfeature_id) &
                                                                          (self.adaptor.location.dbxref_id == self.adaptor.dbxref.dbxref_id)).select(self.adaptor.location.location_id,
                                                                                                                                             self.adaptor.dbxref.dbname,
                                                                                                                                             self.adaptor.dbxref.accession,
                                                                                                                                             self.adaptor.dbxref.accession,))
            lookup = {}
            for location_id, dbname, accession, version in remote_results:
                if version and version != "0":
                    v = "%s.%s" % (accession, version)
                else:
                    v = accession
                # subfeature remote location db_ref are stored as a empty string when
                # not present
                if dbname == "":
                    dbname = None
                lookup[location_id] = (dbname, v)
            
            feature = SeqFeature.SeqFeature(type = seqfeature_type)
            feature._seqfeature_id = seqfeature_id #Store the key as a private property
            feature.qualifiers = qualifiers
            if len(locations) == 0:
                pass
            elif len(locations) == 1:
                location_id, start, end, strand = locations[0]
                #See Bug 2677, we currently don't record the location_operator
                #For consistency with older versions Biopython, default to "".
                feature.location_operator = \
                    self._retrieve_location_qualifier_value( location_id)
                dbname, version = lookup.get(location_id, (None, None))
                feature.location = SeqFeature.FeatureLocation(start, end)
                feature.strand = strand
                feature.ref_db = dbname
                feature.ref = version
            else:
                assert feature.sub_features == []
                for location in locations:
                    location_id, start, end, strand = location
                    dbname, version = lookup.get(location_id, (None, None))
                    subfeature = SeqFeature.SeqFeature()
                    subfeature.type = seqfeature_type
                    subfeature.location_operator = \
                        self._retrieve_location_qualifier_value(location_id)
                    #TODO - See Bug 2677 - we don't yet record location_operator,
                    #so for consistency with older versions of Biopython default
                    #to assuming its a join.
                    if not subfeature.location_operator:
                        subfeature.location_operator="join"
                    subfeature.location = SeqFeature.FeatureLocation(start, end)
                    subfeature.strand = strand
                    subfeature.ref_db = dbname
                    subfeature.ref = version
                    feature.sub_features.append(subfeature)
                # Assuming that the feature loc.op is the same as the sub_feature
                # loc.op:
                feature.location_operator = \
                    feature.sub_features[0].location_operator
                # Locations are in order, but because of remote locations for
                # sub-features they are not necessarily in numerical order:
                start = locations[0][1]
                end = locations[-1][2]
                feature.location = SeqFeature.FeatureLocation(start, end)
                # To get the parent strand (as done when parsing GenBank files),
                # need to consider evil mixed strand examples like this,
                # join(complement(69611..69724),139856..140087,140625..140650)
                strands = set(sf.strand for sf in feature.sub_features)
                if len(strands)==1:
                    feature.strand = feature.sub_features[0].strand
                else:
                    feature.strand = None # i.e. mixed strands
    
            seq_feature_list.append(feature)
    
        return seq_feature_list

    
    def _retrieve_bioentry_relationship(self,bioentry_id):
        """ retrieve the relationships that involve the bioentry_id as a subject"""
        
        objs=((relation.term.name,relation.bioentry_relationship.object_bioentry_id,) \
                     for relation in self.adaptor((self.adaptor.bioentry_relationship.subject_bioentry_id == bioentry_id) & \
                              (self.adaptor.bioentry_relationship.term_id == self.adaptor.term.term_id)\
                              ).select(self.adaptor.term.name,
                                       self.adaptor.bioentry_relationship.object_bioentry_id, 
                                       orderby=self.adaptor.bioentry_relationship.rank))
        relationships = {}
        for relation, obj in objs:
            relation = 'relationship_'+relation.replace(' ','_',)
            if relation not in relationships:
                relationships[relation] = []
            obj_data = self.adaptor(self.adaptor.bioentry.bioentry_id == obj).select()
            if not obj_data:
                raise IndexError('Object bioentry does not exist: %i'%obj)
            obj_biodb = self.adaptor(self.adaptor.biodatabase.biodatabase_id == obj_data[0].biodatabase_id).select(self.adaptor.biodatabase.name)[0].name
            if obj_data[0].accession:
                relationships[relation].append('%s:%s'%(obj_biodb, obj_data[0].accession))
            else:
                relationships[relation].append(str(obj))
        return relationships
    
    def _retrieve_bioentry_relationship_obj(self,bioentry_id):
        """ retrieve the relationships that involve the bioentry_id as an object"""
        
        subjs=((relation.term.name,relation.bioentry_relationship.subject_bioentry_id,) \
                     for relation in self.adaptor((self.adaptor.bioentry_relationship.object_bioentry_id == bioentry_id) & \
                              (self.adaptor.bioentry_relationship.term_id == self.adaptor.term.term_id)\
                              ).select(self.adaptor.term.name,
                                       self.adaptor.bioentry_relationship.subject_bioentry_id, 
                                       orderby=self.adaptor.bioentry_relationship.rank))
        relationships = {}
        for relation, sbj in subjs:
            relation = 'relationship_'+relation.replace(' ','_',)
            if relation not in relationships:
                relationships[relation] = []
            sbj_data = self.adaptor(self.adaptor.bioentry.bioentry_id == sbj).select()
            if not sbj_data:
                raise IndexError('Object bioentry does not exist: %i'%sbj)
            sbj_biodb = self.adaptor(self.adaptor.biodatabase.biodatabase_id == sbj_data[0].biodatabase_id).select(self.adaptor.biodatabase.name)[0].name
            if sbj_data[0].accession:
                relationships[relation].append('%s:%s'%(sbj_biodb, sbj_data[0].accession))
            else:
                relationships[relation].append( str(sbj))
        return relationships
    
    def _retrieve_seqrecord(self, bioentry_id = None, accession = None):
        """create a SeqRecord object from BioSQL.
        use _retrieve_seqrecord(None, accession_code) for accession-based retrieving,
        needed to mantain biopython retrocompatibility"""
        if  accession:
            bioentry_id = self._get_bioentry_id(accession = accession)
        elif not bioentry_id:
            raise IndexError('Specify a bioentry_id or an accession')
        
        bioentry_data= self.adaptor(self.adaptor.bioentry.bioentry_id == bioentry_id).select()
        if not bioentry_data:
            raise IndexError('Provided bioentry does not exist')
        else:
            bioentry_data=bioentry_data[0]
        biodatabase_id = bioentry_data.biodatabase_id
        taxon_id = bioentry_data.taxon_id
        name = bioentry_data.name
        accession = bioentry_data.accession
        version = bioentry_data.version
        identifier = bioentry_data.identifier
        division = bioentry_data.division
        description = bioentry_data.description

        if version and version != "0":
            seqrec_id = "%s.%s" % (accession, version)
        else:
            seqrec_id = accession

        seqrec_seq = self._retrieve_seq(bioentry_id)
        #We don't yet record any per-letter-annotations in the
        #BioSQL database, but we should set this property up
        #for completeness (and the __str__ method).
        try:
            length = len(seqrec_seq)
        except:
            #Could be no sequence in the database!
            length = 0
        seqrec_per_letter_annotations = _RestrictedDict(length=length)
        seqrec_dbxrefs = self._retrieve_dbxrefs(bioentry_id)
        seqrec_features = self._retrieve_features(bioentry_id)
        seqrec_annotations = self._retrieve_annotations(bioentry_id, taxon_id)

        seqrec=SeqRecord(seqrec_seq, id = seqrec_id, name = name, description = description)
        seqrec.dbxrefs = seqrec_dbxrefs
        seqrec.features = seqrec_features
        seqrec.annotations = seqrec_annotations

        return seqrec

    #====================== NEW METHODS =====================================================
    
    def _get_latest_version(self, bioentry_id = None, accession = None):
        '''returns -1 if not found '''
        
        if bioentry_id:
            rows = self.adaptor(self.adaptor.bioentry.bioentry_id == bioentry_id).select()
            if rows:
                accession = rows[0].accession
                biodb = rows[0].biodatabase_id
            else:
                raise IndexError('BIoentry_id %s do not exist'%bioentry_id)
        else:
            biodb = self.dbid
            rows = self.adaptor((self.adaptor.bioentry.biodatabase_id == biodb) &
                                        (self.adaptor.bioentry.accession == accession) 
                                        ).select(self.adaptor.bioentry.version)
            if rows:
                return max([row.version for row in rows])
            else:
                return -1


    def _load_relation_from_annotation (self, bioentry_id, key, values):
        '''create a relationship starting froma a relationship_ annotation'''
        
        relation = ' '.join(key.split('_')[1:])
        for value in values:
            value = value.split(':')
            biodb_name = value[0]
            accession = ':'.join(value[1:]) 
            object_bioentry_id =self._get_bioentry_id(accession = accession, dbid = self.dbs[biodb_name])
            self._add_bioentry_relation(object_bioentry_id, bioentry_id, relation)
        return
    
    def _add_bioentry_relation(self, object_bioentry_id, subject_bioentry_id, relation, rank = 0):
        '''bioentry_id,bioentry_id,'is child of' 
        if the same relation exist no insert is made. because of web2py or because of postgresql?
        '''
        tag_ontology_id = self._get_ontology_id('Bioentry relation')
        term_id = self._get_term_id(relation, ontology_id=tag_ontology_id)
        if rank ==0:
            last_rank = self._get_last_relation_rank(object_bioentry_id, subject_bioentry_id)
            rank = last_rank + 1
        
        if not self._check_bioentry_relation(object_bioentry_id, subject_bioentry_id, relation):
            oid = self.adaptor.bioentry_relationship.insert(object_bioentry_id = object_bioentry_id, 
                                                            subject_bioentry_id = subject_bioentry_id,
                                                            term_id = term_id,
                                                            rank = rank) 
            self.adaptor.commit()
            return oid
        else:    
            raise RelationRedundancyError(object_bioentry_id, subject_bioentry_id, term_id, relation)
    
    def _check_bioentry_relation(self, object_bioentry_id, subject_bioentry_id, relation):
        ''' checks  if a relations is already present in the bioentry_relationship table
        to be used to avoid duplicate relations'''
        tag_ontology_id = self._get_ontology_id('Bioentry relation')
        term_id = self._get_term_id(relation, ontology_id=tag_ontology_id)
        
        if self.adaptor((self.adaptor.bioentry_relationship.subject_bioentry_id == subject_bioentry_id) &
                        (self.adaptor.bioentry_relationship.object_bioentry_id == object_bioentry_id) &
                        (self.adaptor.bioentry_relationship.term_id == term_id)).count():
            return True
        else:
            return False
        
    def _get_last_relation_rank(self,object_bioentry_id, subject_bioentry_id):
        '''returns -1 if not found '''
        rows = self.adaptor((self.adaptor.bioentry_relationship.subject_bioentry_id == subject_bioentry_id) &
                            (self.adaptor.bioentry_relationship.object_bioentry_id == object_bioentry_id)
                           ).select(self.adaptor.bioentry_relationship.rank)
        if rows:
            return max([int(row.rank) for row in rows])
        else:
            return -1

    def _load_bioentry_time_stamp(self, bioentry_id):
        if self.time_stamps_on_biosql:
            self.adaptor.bioentry_timestamp.insert(bioentry_id = bioentry_id)
        else:
            self.timestampsdb.bioentry_timestamp.insert(bioentry_id = bioentry_id)
            
  
    def _set_bioentry_time_stamp(self, bioentry_id, created_by, created_on, modified_by, modified_on):
        '''this is used to override default timestamps,
        needed for archive purpose '''
        if self.time_stamps_on_biosql:
            self.adaptor(self.adaptor.bioentry_timestamp.bioentry_id == bioentry_id).update(created_by = created_by, 
                                                                                            created_on = created_on, 
                                                                                            modified_by = modified_by,
                                                                                            modified_on = modified_on)
#            self.adaptor.bioentry_timestamp[bioentry_id] = dict(created_by = created_by, 
#                                                                created_on = created_on, 
#                                                                modified_by = modified_by,
#                                                                modified_on = modified_on)
        else:
            self.timestampsdb(self.timestampsdb.bioentry_timestamp.bioentry_id == bioentry_id).update(created_by = created_by, 
                                                                                                        created_on = created_on, 
                                                                                                        modified_by = modified_by,
                                                                                                        modified_on = modified_on)
#            self.timestampsdb.bioentry_timestamp[bioentry_id] = dict(created_by = created_by, 
#                                                                     created_on = created_on, 
#                                                                     modified_by = modified_by,
#                                                                     modified_on = modified_on)

    def _update_bioentry_time_stamp(self, bioentry_id, modified_by = current_user, modified_on = now):
        '''updates timestamps to current user and current time.
        cab update to custom user and date if passed '''
        if self.time_stamps_on_biosql:
            self.adaptor(self.adaptor.bioentry_timestamp.bioentry_id == bioentry_id).update(modified_by = modified_by, 
                                                                                                   modified_on = modified_on)
        else:
            self.timestampsdb(self.timestampsdb.bioentry_timestamp.bioentry_id == bioentry_id).update(modified_by = modified_by, 
                                                                                                             modified_on = modified_on)


    def _get_bioentry_time_stamp(self, bioentry_id):
        if self.time_stamps_on_biosql:
            timestamps =  self.adaptor(self.adaptor.bioentry_timestamp.bioentry_id == bioentry_id).select().first()
        else:
            timestamps =  self.timestampsdb(self.timestampsdb.bioentry_timestamp.bioentry_id == bioentry_id).select().first()
        if timestamps:
            return dict(created_by = timestamps.created_by,
                        created_on = timestamps.created_on, 
                        modified_by = timestamps.modified_by, 
                        modified_on = timestamps.modified_on)
                                                            
    def _get_bioentry_biodatabase_name(self, bioentry_id):
        rows = self.adaptor((self.adaptor.bioentry.bioentry_id == bioentry_id) &
                            (self.adaptor.bioentry.biodatabase_id == self.adaptor.biodatabase.biodatabase_id)
                           ).select(self.adaptor.biodatabase.name)
        if rows:
            return rows[0].name
        else:
            raise IndexError('bioentry_id %s not present in db'%bioentry_id)

    def get_bioentries_with_feature(self, type, start=0,end=-1):
        '''Returns the list of bioentries in the current biodatabase 
        having the the requested feature type '''

        
        bioentries = [row.bioentry_id for row in self.adaptor((self.adaptor.bioentry.biodatabase_id == self.dbid) &
                                                                (self.adaptor.bioentry.bioentry_id == self.adaptor.seqfeature.bioentry_id) &
                                                                (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id) &
                                                                (self.adaptor.term.name == type)).select(self.adaptor.bioentry.bioentry_id)]
                            
        return bioentries
        
    
    
    def get_bioentries_with_annotation(self, type, value=''):
        '''Returns the list of bioentries in the current biodatabase 
        having the the requested annotation type
        if a value is passed only entries matching that value are returned
        the value matching is done with LIKE and is case INsensitive'''

        if value:
            value_query = (self.adaptor.bioentry_qualifier_value.value.lower().like('%%%s%%'%value.lower()))
        else:
            value_query = (self.adaptor.bioentry_qualifier_value.bioentry_id >0)
        bioentries = [row.bioentry_id for row in self.adaptor((self.adaptor.bioentry.biodatabase_id == self.dbid) &
                                                                (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_qualifier_value.bioentry_id) &
                                                                (self.adaptor.bioentry_qualifier_value.term_id == self.adaptor.term.term_id) &
                                                                value_query &
                                                                (self.adaptor.term.name == type)).select(self.adaptor.bioentry.bioentry_id)]
                            
        return bioentries
    
    def _get_bioentry_id(self, accession = None, id = None, dbid = None):
        '''returns bioentry id, given an accession or an identifier
        if multiple versions are available for the query the last one is returned'''
        if not dbid:
            dbid = self.dbid
        
        if accession:
            if accession.count(".") == 1: # try to get a version from the id
                #This assumes the string is something like "XXXXXXXX.123"
                original_accession=accession
                accession, version = accession.split('.')
                try:
                    version = int(version)
                except ValueError:
                    accession = original_accession
                    version = 0
            else: # otherwise just use a version of 0
                accession = accession
                version = 0
            if version:
                bioentry = self.adaptor((self.adaptor.bioentry.biodatabase_id == dbid) &
                                        (self.adaptor.bioentry.accession == accession) &
                                        (self.adaptor.bioentry.version == version) 
                                        ).select(self.adaptor.bioentry.bioentry_id)
            else:#return the highest available version
                 bioentry = self.adaptor((self.adaptor.bioentry.biodatabase_id == dbid) &
                                        (self.adaptor.bioentry.accession == accession)
                                        ).select(self.adaptor.bioentry.bioentry_id, orderby = ~self.adaptor.bioentry.version)
                                        
            if bioentry:
                return bioentry[0].bioentry_id
            else:
                raise ValueError('No entry found with accession: %s'%(accession))
            
        elif id:
            bioentry = self.adaptor((self.adaptor.bioentry.biodatabase_id == self.dbid) &
                                    (self.adaptor.bioentry.identifier == id) 
                                    ).select(self.adaptor.bioentry.bioentry_id)
            if bioentry:
                return bioentry[0].bioentry_id
            else:
                raise ValueError('No entry found with identifier: %s'%(accession))
            
    def get_bioentry_names(self):
        '''returns  a dictionary of bioentries id:name for the current set db'''
        return  dict([(row.bioentry_id, row.name)  for row in  self.adaptor(self.adaptor.bioentry.biodatabase_id == self.dbid).select(self.adaptor.bioentry.bioentry_id,
                                                                                                                    self.adaptor.bioentry.name,)])

    def get_annotation(self, bioentry_id, qualifier_term_name):
        
        
        return
    
    




    
    