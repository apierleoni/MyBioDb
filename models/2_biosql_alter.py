__author__ = 'pierleonia'



# biopython
from Bio import Alphabet
from Bio.SeqUtils.CheckSum import crc64
from Bio import Entrez
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio import SeqFeature


class BaseBioSQLAlter():
    '''Base class for to handlers for BioSQL data '''
    def __init__(self, handler,):
        self.handler = handler
        self.adaptor = handler.adaptor

    def create(self,):
        raise NotImplementedError('This method is not available')

    def read(self,):
        raise NotImplementedError('This method is not available')

    def update(self,):
        raise NotImplementedError('This method is not available')

    def delete(self,):
        raise NotImplementedError('This method is not available')





class BioSQLAlterFeature(BaseBioSQLAlter):

#    def __init__(self, handler):
#        '''handler is an initiated BioSQLHandler  object
#        bioentry_id required
#        seqfeature_id if 0 means new, else is the id in the database
#        feature  pass to upload a new one
#
#        '''
#        BaseBioSQLAlter.__init__(self, handler = handler)




    def create(self,feature, bioentry_id):
        """
        Load a new biopython SeqFeature into the database (PRIVATE).
        """
        ranks = [row.rank for row in self.adaptor((self.adaptor.seqfeature.bioentry_id == bioentry_id)
                                                 ).select(self.adaptor.seqfeature.rank)]
        if ranks:
            rank = max(ranks)+1
        else:
            rank = 1
        seqfeature_id = self._load_seqfeature_basic(feature.type, rank, bioentry_id)
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

        return seqfeature_id

    def read(self, seqfeature_id):

        '''get id and type'''
        rows = self.adaptor((self.adaptor.seqfeature.seqfeature_id == seqfeature_id) &\
                            (self.adaptor.seqfeature.type_term_id == self.adaptor.term.term_id)
                           ).select(self.adaptor.seqfeature.seqfeature_id,self.adaptor.term.name)

        if rows:
            feature_type = rows[0].term.name
        else:
            raise IndexError('No feature available for seqfeature_id %i'%(seqfeature_id))

        '''get qualifiers'''
        qvs=((qualifier.term.name,qualifier.seqfeature_qualifier_value.value,)\
            for qualifier in self.adaptor((self.adaptor.seqfeature_qualifier_value.seqfeature_id == seqfeature_id) &\
                                          (self.adaptor.seqfeature_qualifier_value.term_id == self.adaptor.term.term_id)\
        ).select(self.adaptor.term.name,self.adaptor.seqfeature_qualifier_value.value,orderby=self.adaptor.seqfeature_qualifier_value.rank))

        qualifiers = {}
        for qv_name, qv_value in qvs:
            qualifiers.setdefault(qv_name, []).append(qv_value)
            # Get db_xrefs [special case of qualifiers]
        qvs = ((row.dbname,row.accession)\
            for row in self.adaptor((self.adaptor.seqfeature_dbxref.seqfeature_id == seqfeature_id) &
                                    (self.adaptor.seqfeature_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id)
        ).select(self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,orderby=self.adaptor.seqfeature_dbxref.rank))
        for qv_name, qv_value in qvs:
            value = "%s:%s" % (qv_name, qv_value)
            qualifiers.setdefault("db_xref", []).append(value)




        ''' Get locations '''
        results= ((row.location_id,
                   row.start_pos,
                   row.end_pos,
                   row.strand) for row in self.adaptor(self.adaptor.location.seqfeature_id == seqfeature_id).select(
                                                        self.adaptor.location.location_id,
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
                raise ValueError("Invalid strand %s found in database for "\
                                 "seqfeature_id %s" % (strand, seqfeature_id))
            if end < start:
                import warnings
                warnings.warn("Inverted location start/end (%i and %i) for "\
                              "seqfeature_id %s" % (start, end, seqfeature_id))
            locations.append( (location_id, start, end, strand) )


        ''' Get possible remote reference information'''
        remote_results = ((row.location.location_id,
                           row.dbxref.dbname,
                           row.dbxref.accession,
                           row.dbxref.accession,) for row in self.adaptor((self.adaptor.location.seqfeature_id == seqfeature_id) &
                                                                          (self.adaptor.location.dbxref_id == self.adaptor.dbxref.dbxref_id)).select(self.adaptor.location.location_id,
            self.adaptor.dbxref.dbname,
            self.adaptor.dbxref.accession,
            self.adaptor.dbxref.accession,))
        ref_lookup = {}
        for location_id, dbname, accession, version in remote_results:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
                # subfeature remote location db_ref are stored as a empty string when
            # not present
            if dbname == "":
                dbname = None
            ref_lookup[location_id] = (dbname, v)

        feature = SeqFeature.SeqFeature(type =feature_type)
        feature._seqfeature_id = seqfeature_id #Store the key as a private property
        feature.qualifiers = qualifiers
        if len(locations) == 0:
            pass
        elif len(locations) == 1:
            location_id, start, end, strand = locations[0]
            #See Bug 2677, we currently don't record the location_operator
            #For consistency with older versions Biopython, default to "".
            feature.location_operator =\
            self.handler._retrieve_location_qualifier_value(location_id)
            dbname, version = ref_lookup.get(location_id, (None, None))
            feature.location = SeqFeature.FeatureLocation(start, end)
            feature.strand = strand
            feature.ref_db = dbname
            feature.ref = version
        else:
            assert feature.sub_features == []
            for location in locations:
                location_id, start, end, strand = location
                dbname, version = ref_lookup.get(location_id, (None, None))
                subfeature = SeqFeature.SeqFeature()
                subfeature.type = feature_type
                subfeature.location_operator =\
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
                feature.sub_features.append(subfeature)
                # Assuming that the feature loc.op is the same as the sub_feature
            # loc.op:
            feature.location_operator =\
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

        return feature

    def update(self, feature, seqfeature_id):
        '''do the same as insert except for _load_seqfeature_basic to maintain seqfeature_id'''

        seqfeature_row = self.adaptor((self.adaptor.seqfeature.seqfeature_id == seqfeature_id)
                     ).select(self.adaptor.seqfeature.rank, self.adaptor.seqfeature.bioentry_id)
        self._update_seqfeature_basic(seqfeature_id,
            feature.type,
            seqfeature_row.rank,
            seqfeature_row.bioentry_id)

        '''remove dbxrefs '''
        self.adaptor(self.adaptor.seqfeature_dbxref.seqfeature_id == seqfeature_id).delete()
        '''remove qualifiers '''
        self.adaptor(self.adaptor.seqfeature_qualifier_value.seqfeature_id == seqfeature_id).delete()
        '''remove locations '''
        self.adaptor(self.adaptor.location.seqfeature_id == seqfeature_id).delete()
        '''Add updated location and qualifiers '''
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

        self.handler._update_bioentry_time_stamp(seqfeature_row.bioentry_id)


    def delete(self, seqfeature_id):
        '''remove all the data related to a given seqfeature_id
        terms and dbxrefs will remain in the respective tables on the db,
        since they can be used by other entities, but will no longer be linked to
        the deleted seqfeature '''

        bioentry_id = self.adaptor.seqfeature[seqfeature_id].bioentry_id
        '''remove dbxrefs '''
        self.adaptor(self.adaptor.seqfeature_dbxref.seqfeature_id == seqfeature_id).delete()

        '''remove qualifiers '''
        self.adaptor(self.adaptor.seqfeature_qualifier_value.seqfeature_id == seqfeature_id).delete()

        '''remove locations '''
        self.adaptor(self.adaptor.location.seqfeature_id == seqfeature_id).delete()

        '''remove seqfeature '''
        self.adaptor(self.adaptor.seqfeature.seqfeature_id == seqfeature_id).delete()

        self.handler._update_bioentry_time_stamp(bioentry_id)


    def update_seqfeature_type(self, seqfeature_id, new_value):
        ontology_id = self.handler._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self.handler._get_term_id(new_value,
            ontology_id = ontology_id)

        self.adaptor.seqfeature[seqfeature_id] = dict (type_term_id = seqfeature_key_id)

    def update_seqfeature_start(self, seqfeature_id, new_value):
        raise NotImplementedError('TODO')
    def update_seqfeature_end(self, seqfeature_id, new_value):
        raise NotImplementedError('TODO')
    def update_seqfeature_identifier(self, seqfeature_id, new_value):
        self._update_qualifier(seqfeature_id, 'id', new_value)
    def update_seqfeature_description(self, seqfeature_id, new_value):
        self._update_qualifier(seqfeature_id, 'description', new_value)


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
            warnings.warn("%s location operators are not fully supported"\
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

    def _update_qualifier(self, seqfeature_id, qualifier_key, qualifier_values):
        tag_ontology_id = self.handler._get_ontology_id('Annotation Tags')
        if qualifier_key != 'db_xref':
            qualifier_key_id = self.handler._get_term_id(qualifier_key,
                ontology_id=tag_ontology_id)
            # now add all of the values to their table
            entries = qualifier_values
            if not isinstance(entries, list):
                # Could be a plain string, or an int or a float.
                # However, we exect a list of strings here.
                entries = [entries]
            for qual_value_rank in range(len(entries)):
                qualifier_value = entries[qual_value_rank]
                stored_value = self.adaptor((self.adaptor.seqfeature_qualifier_value.seqfeature_id == seqfeature_id) &\
                                            (self.adaptor.seqfeature_qualifier_value.term_id == qualifier_key_id) &\
                                            (self.adaptor.seqfeature_qualifier_value.rank == qual_value_rank + 1))
                if qualifier_value:
                    if stored_value.count():
                        stored_value.update(value = qualifier_value)
                    else:
                        self.adaptor.seqfeature_qualifier_value.insert(seqfeature_id = seqfeature_id,
                            term_id = qualifier_key_id,
                            rank = qual_value_rank + 1,
                            value = qualifier_value)
                else:
                    stored_value.delete()



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


class BioSQLAlterBioentryFeatures(BaseBioSQLAlter):
    ''' TO DO
        '''

    def __init__(self, handler, bioentry_id,):
        BaseBioSQLAlter.__init__(self, handler = handler)
        self.bioentry_id = bioentry_id
        self.features = dict() #{rank:{seqfeature_id, feature}} dictionary
        self._get()

    def _get(self):
        '''get all avilable features '''
        rows = self.adaptor((self.adaptor.seqfeature.bioentry_id == self.bioentry_id)).select()
        featureHandler = BioSQLAlterFeature(handler = self.handler)
        if rows:
            for row in rows:
                if row.rank in self.features:
                    raise ValueError('multiple seqfeatures present with the same rank')
                else:
                    self.features[row.rank] = dict(seqfeature_id = row.seqfeature_id,
                                                   seqfeature = featureHandler.read(seqfeature_id = row.seqfeature_id))
        self.features

    def get_seqfeatures(self):
        ordered = []
        for key in sorted(self.features):
            ordered.append((self.features[key]['seqfeature_id'], self.features[key]['seqfeature']))
        return ordered


class BioSQLAlterBioentryDBXref(BaseBioSQLAlter):

    def create(self):
        raise NotImplementedError('This method is not available')

    def read(self,):
        raise NotImplementedError('This method is not available')

    def update(self,dbxref_id, dbname, accession):
        self.adaptor.dbxref[int(dbxref_id)] = dict (dbname = dbname, accession = accession)
    def update_dbname(self,dbxref_id, dbname):
        self.adaptor.dbxref[int(dbxref_id)] = dict (dbname = dbname)
    def update_accession(self,dbxref_id, accession):
        self.adaptor.dbxref[int(dbxref_id)] = dict ( accession = accession)

    def delete(self,bioentry_id, dbxref_id):
        self.adaptor((self.adaptor.bioentry_dbxref.bioentry_id == bioentry_id) &\
                     (self.adaptor.bioentry_dbxref.dbxref_id == dbxref_id)).delete()


    def _get(self,):
        dbxrefs=((dbxref.dbname, dbxref.accession, dbxref.version) for dbxref in self.adaptor((self.adaptor.bioentry.bioentry_id == self.bioentry_id) &\
                                                                                              (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) &\
                                                                                              (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id)).select(
            self.adaptor.dbxref.dbname,self.adaptor.dbxref.accession,self.adaptor.dbxref.version, orderby=self.adaptor.bioentry_dbxref.rank))
        for dbname, accession, version in dbxrefs:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            self.dbxrefs.append("%s:%s" % (dbname, v))


class BioSQLAlterBioentryDBXrefs(BaseBioSQLAlter):

    def create(self):
        raise NotImplementedError('This method is not available')

    def read(self,bioentry_id):
        dbxrefs=((dbxref.dbname, dbxref.accession, dbxref.version, dbxref.dbxref_id) for dbxref in self.adaptor((self.adaptor.bioentry.bioentry_id == bioentry_id) &\
                                                                                              (self.adaptor.bioentry.bioentry_id == self.adaptor.bioentry_dbxref.bioentry_id) &\
                                                                                              (self.adaptor.bioentry_dbxref.dbxref_id == self.adaptor.dbxref.dbxref_id)).select(
                                                                                                            self.adaptor.dbxref.dbname,
                                                                                                            self.adaptor.dbxref.dbxref_id,
                                                                                                            self.adaptor.dbxref.accession,
                                                                                                            self.adaptor.dbxref.version,
                                                                                                            orderby=self.adaptor.bioentry_dbxref.rank))
        returnlist = []
        for dbname, accession, version, dbxref_id in dbxrefs:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            returnlist.append(("%s:%s" % (dbname, v), dbxref_id))
        return returnlist

    def update(self,):
        raise NotImplementedError('This method is not available')

    def delete(self,):
        raise NotImplementedError('This method is not available')